!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cumulus Parameterization
!!
!! @par Description
!!          Cumulus convection by Kain-Fritsch parameterization
!!          Reference: Kain and Fritsch(1990)
!!                     Kain (2004)
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2016-06-27 (S.Matsugishi) [new]
!!
!<
#include "inc_openmp.h"
module scale_atmos_phy_cp_kf
  !------------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index

  use scale_tracer, only: QA
  !------------------------------------------------------------------------------
  implicit none
  private
  !------------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_CP_kf_setup
  public :: ATMOS_PHY_CP_kf

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private, allocatable :: lifetime  (:,:) ! convectime lifetime [s]
  integer , private, allocatable :: I_convflag(:,:) ! convection type 0:deep convection 1:shallow convection 2: no convection

  ! kf time controll
  integer,  private :: TIME_RES_KF   ! time step for kf
  integer,  private :: TIME_DSTEP_KF ! time interval
  logical,  private :: TIME_DOKF     ! exclude kf trigger

  ! tuning parameter
  logical,  private :: PARAM_ATMOS_PHY_CP_kf_wadapt = .true.
  integer,  private :: PARAM_ATMOS_PHY_CP_kf_w_time = 16

  !------------------------------------------------------------------------------
contains
  !------------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_CP_kf_setup (CP_TYPE)
    use scale_process, only: &
       PRC_MPIstop
    use scale_time , only :&
       TIME_DTSEC,             &
       TIME_DTSEC_ATMOS_PHY_CP
    use scale_atmos_phy_cp_kf_sub,only: &
       kf_lutab, &
       CP_kf_param
    implicit none

    character(len=*), intent(in) :: CP_TYPE

    ! tunning parameters, original parameter set is from KF2004(WRF) and JMA-NHM
    real(RP) :: PARAM_ATMOS_PHY_CP_kf_rate      =   0.03_RP ! ratio of cloud water and precipitation (Ogura and Cho 1973)
    integer  :: PARAM_ATMOS_PHY_CP_kf_trigger   = 1         ! trigger function type 1:WRF 3:JMA
    logical  :: PARAM_ATMOS_PHY_CP_kf_qs        = .true.    ! qs is allowed?
    logical  :: PARAM_ATMOS_PHY_CP_kf_qi        = .true.    ! qi is allowed?
    real(DP) :: PARAM_ATMOS_PHY_CP_kf_dt        =    5.0_DP ! KF convection check time interval [min]
    real(RP) :: PARAM_ATMOS_PHY_CP_kf_dlcape    =    0.1_RP ! cape decleace rate
    real(RP) :: PARAM_ATMOS_PHY_CP_kf_dlifetime = 1800.0_RP ! minimum lifetime scale of deep convection [sec]
    real(RP) :: PARAM_ATMOS_PHY_CP_kf_slifetime = 2400.0_RP ! lifetime of shallow convection [sec]
    real(RP) :: PARAM_ATMOS_PHY_CP_kf_DEPTH_USL =  300.0_RP ! depth of updraft source layer [hPa]
    integer  :: PARAM_ATMOS_PHY_CP_kf_cond      = 1         ! condload type 1:Ogura-Cho(1973) 2:Kessler
    real(RP) :: PARAM_ATMOS_PHY_CP_kf_thres     = 1.E-3_RP  ! autoconversion rate for Kessler
    logical  :: PARAM_ATMOS_PHY_CP_kf_warmrain  = .false.   ! QQA is less equal to 3?
    logical  :: PARAM_ATMOS_PHY_CP_kf_LOG       = .false.   ! output log?

    NAMELIST / PARAM_ATMOS_PHY_CP_KF / &
       PARAM_ATMOS_PHY_CP_kf_rate,      &
       PARAM_ATMOS_PHY_CP_kf_trigger,   &
       PARAM_ATMOS_PHY_CP_kf_qs,        &
       PARAM_ATMOS_PHY_CP_kf_qi,        &
       PARAM_ATMOS_PHY_CP_kf_dt,        &
       PARAM_ATMOS_PHY_CP_kf_dlcape,    &
       PARAM_ATMOS_PHY_CP_kf_dlifetime, &
       PARAM_ATMOS_PHY_CP_kf_slifetime, &
       PARAM_ATMOS_PHY_CP_kf_DEPTH_USL, &
       PARAM_ATMOS_PHY_CP_kf_cond,      &
       PARAM_ATMOS_PHY_CP_kf_thres,     &
       PARAM_ATMOS_PHY_CP_kf_warmrain,  &
       PARAM_ATMOS_PHY_CP_kf_LOG,       &
       PARAM_ATMOS_PHY_CP_kf_wadapt, &
       PARAM_ATMOS_PHY_CP_kf_w_time

    integer :: ierr
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[CUMULUS] / Categ[ATMOS PHYSICS] / Origin[SCALElib]'
    if( IO_L ) write(IO_FID_LOG,*) '+++ Kain-Fritsch scheme'

    if ( CP_TYPE /= 'KF' ) then
       write(*,*) 'xxx ATMOS_PHY_CP_TYPE is not KF. Check!'
       call PRC_MPIstop
    endif

    if ( abs(TIME_DTSEC_ATMOS_PHY_CP-TIME_DTSEC) > 0.0_DP ) then
       write(*,*) 'xxx TIME_DTSEC_ATMOS_PHY_CP should be same as TIME_DTSEC for KF scheme. STOP.'
       call PRC_MPIstop
    endif

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_CP_KF,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_CP_KF. Check!'
       call PRC_MPIstop
    endif
    if( IO_L ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_CP_KF)

    call kf_lutab ! set kf look up table

    ! set kf convection check step
    TIME_DSTEP_KF = nint( PARAM_ATMOS_PHY_CP_kf_dt * 60.0_DP / TIME_DTSEC_ATMOS_PHY_CP )
    TIME_DSTEP_KF = max(TIME_DSTEP_KF,1) ! kf time interval step
    TIME_RES_KF   = -1                   ! initialize to keep consistent for below step
    TIME_DOKF     = .true.               ! initialize

    call CP_kf_param( PARAM_ATMOS_PHY_CP_kf_rate,      &
                      PARAM_ATMOS_PHY_CP_kf_trigger,   &
                      PARAM_ATMOS_PHY_CP_kf_qs,        &
                      PARAM_ATMOS_PHY_CP_kf_qi,        &
                      PARAM_ATMOS_PHY_CP_kf_dt,        &
                      PARAM_ATMOS_PHY_CP_kf_dlcape,    &
                      PARAM_ATMOS_PHY_CP_kf_dlifetime, &
                      PARAM_ATMOS_PHY_CP_kf_slifetime, &
                      PARAM_ATMOS_PHY_CP_kf_DEPTH_USL, &
                      PARAM_ATMOS_PHY_CP_kf_cond,      &
                      PARAM_ATMOS_PHY_CP_kf_thres,     &
                      PARAM_ATMOS_PHY_CP_kf_warmrain,  &
                      PARAM_ATMOS_PHY_CP_kf_LOG ,      &
                      TIME_DSTEP_KF                    )

    ! output parameter lists
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) "*** Interval for check [step]                       : ", TIME_DSTEP_KF
    if( IO_L ) write(IO_FID_LOG,*) "*** Ogura-Cho condense material convert rate        : ", PARAM_ATMOS_PHY_CP_kf_rate
    if( IO_L ) write(IO_FID_LOG,*) "*** Trigger function type, 1:WRF 3:JMA              : ", PARAM_ATMOS_PHY_CP_kf_trigger
    if( IO_L ) write(IO_FID_LOG,*) "*** Exist qi?                                       : ", PARAM_ATMOS_PHY_CP_kf_qi
    if( IO_L ) write(IO_FID_LOG,*) "*** Exist qs?                                       : ", PARAM_ATMOS_PHY_CP_kf_qs
    if( IO_L ) write(IO_FID_LOG,*) "*** CAPE decrease rate                              : ", PARAM_ATMOS_PHY_CP_kf_dlcape
    if( IO_L ) write(IO_FID_LOG,*) "*** Minimum lifetime scale of deep convection [sec] : ", PARAM_ATMOS_PHY_CP_kf_dlifetime
    if( IO_L ) write(IO_FID_LOG,*) "*** Lifetime of shallow convection            [sec] : ", PARAM_ATMOS_PHY_CP_kf_slifetime
    if( IO_L ) write(IO_FID_LOG,*) "*** Updraft souce layer depth                 [hPa] : ", PARAM_ATMOS_PHY_CP_kf_DEPTH_USL
    if( IO_L ) write(IO_FID_LOG,*) "*** Condload type 1:Ogura-Cho(1973) 2:Kessler       : ", PARAM_ATMOS_PHY_CP_kf_cond
    if( IO_L ) write(IO_FID_LOG,*) "*** Kessler type condload's threshold               : ", PARAM_ATMOS_PHY_CP_kf_thres
    if( IO_L ) write(IO_FID_LOG,*) "*** Warm rain?                                      : ", PARAM_ATMOS_PHY_CP_kf_warmrain
    if( IO_L ) write(IO_FID_LOG,*) "*** Use running mean of w in adaptive timestep?     : ", PARAM_ATMOS_PHY_CP_kf_wadapt
    if( IO_L ) write(IO_FID_LOG,*) "*** Fixed time scale for running mean of w          : ", PARAM_ATMOS_PHY_CP_kf_w_time
    if( IO_L ) write(IO_FID_LOG,*) "*** Output log?                                     : ", PARAM_ATMOS_PHY_CP_kf_LOG

    ! output variables
    allocate( lifetime  (IA,JA) )
    allocate( I_convflag(IA,JA) )
    lifetime  (:,:) = 0.0_RP
    I_convflag(:,:) = 2

    return
  end subroutine ATMOS_PHY_CP_kf_setup

  !------------------------------------------------------------------------------
  subroutine ATMOS_PHY_CP_kf( &
       DENS,           &
       MOMZ,           &
       MOMX,           &
       MOMY,           &
       RHOT,           &
       QTRC,           &
       DENS_t_CP,      &
       MOMZ_t_CP,      &
       MOMX_t_CP,      &
       MOMY_t_CP,      &
       RHOT_t_CP,      &
       RHOQ_t_CP,      &
       MFLX_cloudbase, &
       SFLX_convrain,  &
       cloudtop,       &
       cloudbase,      &
       cldfrac_dp,     &
       cldfrac_sh,     &
       nca,            &
       w0avg           )
    use scale_grid_index
    use scale_tracer, only: &
       MP_QA
    use scale_time, only: &
       TIME_DTSEC_ATMOS_PHY_CP
    use scale_atmos_phy_cp_kf_sub, only: &
       CP_kf_main
    use scale_history, only: &
       HIST_in
    implicit none

    real(RP), intent(in)    :: DENS          (KA,IA,JA)
    real(RP), intent(in)    :: MOMX          (KA,IA,JA)
    real(RP), intent(in)    :: MOMY          (KA,IA,JA)
    real(RP), intent(in)    :: MOMZ          (KA,IA,JA)
    real(RP), intent(in)    :: RHOT          (KA,IA,JA)
    real(RP), intent(in)    :: QTRC          (KA,IA,JA,QA)
    real(RP), intent(inout) :: DENS_t_CP     (KA,IA,JA)
    real(RP), intent(inout) :: MOMZ_t_CP     (KA,IA,JA)    ! not used
    real(RP), intent(inout) :: MOMX_t_CP     (KA,IA,JA)    ! not used
    real(RP), intent(inout) :: MOMY_t_CP     (KA,IA,JA)    ! not used
    real(RP), intent(inout) :: RHOT_t_CP     (KA,IA,JA)
    real(RP), intent(inout) :: RHOQ_t_CP     (KA,IA,JA,QA)
    real(RP), intent(inout) :: MFLX_cloudbase(IA,JA)       ! not used
    real(RP), intent(inout) :: SFLX_convrain (IA,JA)       ! convective rain rate [kg/m2/s]
    real(RP), intent(inout) :: cloudtop      (IA,JA)       ! cloud top height  [m]
    real(RP), intent(inout) :: cloudbase     (IA,JA)       ! cloud base height [m]
    real(RP), intent(inout) :: cldfrac_dp    (KA,IA,JA)    ! cloud fraction (deep convection)
    real(RP), intent(inout) :: cldfrac_sh    (KA,IA,JA)    ! cloud fraction (shallow convection)
    real(RP), intent(inout) :: nca           (IA,JA)       ! convection active time [sec]
    real(RP), intent(inout) :: w0avg         (KA,IA,JA)    ! running mean of vertical velocity [m/s]

    real(RP) :: U      (KA,IA,JA)       ! x-direction velocity [m/s]
    real(RP) :: V      (KA,IA,JA)       ! y-direction velocity [m/s]
    real(RP) :: TEMP   (KA,IA,JA)       ! temperature [K]
    real(RP) :: PRES   (KA,IA,JA)       ! pressure [Pa]
    real(RP) :: QV     (KA,IA,JA)       ! water vaper mixing ratio [kg/kg]
    real(RP) :: QHYD   (KA,IA,JA,MP_QA) ! water mixing ratio [kg/kg]
    real(RP) :: QSAT   (KA,IA,JA)       ! saturate water vaper mixing ratio [kg/kg]
    real(RP) :: rh     (KA,IA,JA)       ! relative humidity [%]
    real(RP) :: deltap (KA,IA,JA)       ! pressure interval [hPa]
    real(RP) :: deltaz (KA,IA,JA)       ! height interval (center level) [m]
    real(RP) :: Z      (KA,IA,JA)       ! centerlevel real height [m]
    real(RP) :: deltax                  ! delta x [m]

    real(RP) :: cldfrac(KA,2)           ! 1 shallow , 2 deep

    integer  :: i, j
    !------------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Cumulus Parameterization(KF)'

    call KF_wmean( w0avg(:,:,:), & ! [OUT]
                   DENS (:,:,:), & ! [IN]
                   MOMZ (:,:,:)  ) ! [IN]

    TIME_DOKF   = .false.
    TIME_RES_KF = TIME_RES_KF + 1
    if ( TIME_RES_KF == TIME_DSTEP_KF ) then
       TIME_DOKF   = .true.
       TIME_RES_KF = 0
    endif

    if ( TIME_DOKF ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** KF Convection Check '

       ! convert variables
       call CP_kf_init( DENS  (:,:,:),   & ! [IN]
                        MOMZ  (:,:,:),   & ! [IN]
                        MOMX  (:,:,:),   & ! [IN]
                        MOMY  (:,:,:),   & ! [IN]
                        RHOT  (:,:,:),   & ! [IN]
                        QTRC  (:,:,:,:), & ! [IN]
                        U     (:,:,:),   & ! [OUT]
                        V     (:,:,:),   & ! [OUT]
                        TEMP  (:,:,:),   & ! [OUT]
                        PRES  (:,:,:),   & ! [OUT]
                        QV    (:,:,:),   & ! [OUT]
                        QHYD  (:,:,:,:), & ! [OUT]
                        QSAT  (:,:,:),   & ! [OUT]
                        rh    (:,:,:),   & ! [OUT]
                        deltap(:,:,:),   & ! [OUT]
                        deltaz(:,:,:),   & ! [OUT]
                        Z     (:,:,:),   & ! [OUT]
                        deltax           ) ! [OUT]

       call PROF_rapstart('CP_kf', 3)

       do j = JS, JE
       do i = IS, IE
          nca(i,j) = nca(i,j) - real(TIME_DSTEP_KF,RP) * TIME_DTSEC_ATMOS_PHY_CP

          if ( nca(i,j) < 0.5_DP * TIME_DTSEC_ATMOS_PHY_CP ) then ! check convection

             ! calc cumulus convection
             call CP_kf_main( DENS         (:,i,j),   & ! [IN]
                              RHOT         (:,i,j),   & ! [IN]
                              QTRC         (:,i,j,:), & ! [IN]
                              w0avg        (:,i,j),   & ! [IN]
                              U            (:,i,j),   & ! [IN]
                              V            (:,i,j),   & ! [IN]
                              TEMP         (:,i,j),   & ! [IN]
                              PRES         (:,i,j),   & ! [IN]
                              QV           (:,i,j),   & ! [IN]
                              QHYD         (:,i,j,:), & ! [INOUT]
                              QSAT         (:,i,j),   & ! [IN]
                              rh           (:,i,j),   & ! [IN]
                              deltap       (:,i,j),   & ! [IN]
                              deltaz       (:,i,j),   & ! [IN]
                              Z            (:,i,j),   & ! [IN]
                              deltax,                 & ! [IN]
                              nca          (i,j),     & ! [INOUT]
                              DENS_t_CP    (:,i,j),   & ! [OUT]
                              RHOT_t_CP    (:,i,j),   & ! [OUT]
                              RHOQ_t_CP    (:,i,j,:), & ! [OUT]
                              SFLX_convrain(i,j),     & ! [OUT]
                              cldfrac      (:,:),     & ! [OUT]
                              lifetime     (i,j),     & ! [OUT]
                              cloudtop     (i,j),     & ! [OUT]
                              cloudbase    (i,j),     & ! [OUT]
                              I_convflag   (i,j)      ) ! [OUT]

             cldfrac_sh(:,i,j) = cldfrac(:,1)
             cldfrac_dp(:,i,j) = cldfrac(:,2)

          endif
       enddo
       enddo

       call PROF_rapend('CP_kf', 3)
    endif

    call HIST_in( lifetime(:,:),            'KF_LIFETIME', 'lifetime of KF scheme', 's' )
    call HIST_in( real(I_convflag(:,:),RP), 'KF_CONVFLAG', 'CONVECTION FLAG',       ''  )

    return
  end subroutine ATMOS_PHY_CP_kf

  !------------------------------------------------------------------------------
  !> convert variables
  subroutine CP_kf_init( &
       DENS,   &
       MOMZ,   &
       MOMX,   &
       MOMY,   &
       RHOT,   &
       QTRC,   &
       Z,      &
       U,      &
       V,      &
       TEMP,   &
       PRES,   &
       QV,     &
       QHYD,   &
       QSAT,   &
       rh,     &
       deltap, &
       deltaz, &
       deltax  )
    use scale_tracer, only: &
       MP_QA,    &
       I_MP2ALL, &
       I_QV
    use scale_const, only: &
       GRAV => CONST_GRAV, &
       R    => CONST_Rdry
    use scale_grid,only: &
       DX => DX, &
       DY => DY
    use scale_grid_real, only: &
       CZ => REAL_CZ, &
       FZ => REAL_FZ
    use scale_atmos_thermodyn, only: &
       THERMODYN_temp_pres   => ATMOS_THERMODYN_temp_pres,   &
       THERMODYN_rhoe        => ATMOS_THERMODYN_rhoe,        &
       THERMODYN_temp_pres_E => ATMOS_THERMODYN_temp_pres_E, &
       THERMODYN_qd          => ATMOS_THERMODYN_qd
    use scale_atmos_saturation ,only :&
       SATURATION_psat_liq => ATMOS_SATURATION_psat_liq
    implicit none

    real(RP), intent(in)  :: DENS  (KA,IA,JA)
    real(RP), intent(in)  :: MOMX  (KA,IA,JA)
    real(RP), intent(in)  :: MOMY  (KA,IA,JA)
    real(RP), intent(in)  :: MOMZ  (KA,IA,JA)
    real(RP), intent(in)  :: RHOT  (KA,IA,JA)
    real(RP), intent(in)  :: QTRC  (KA,IA,JA,QA)
    real(RP), intent(out) :: Z     (KA,IA,JA)       ! z level (center level)
    real(RP), intent(out) :: U     (KA,IA,JA)       ! x-direction velocity
    real(RP), intent(out) :: V     (KA,IA,JA)       ! y-direction velocity
    real(RP), intent(out) :: TEMP  (KA,IA,JA)       ! temperature [K]
    real(RP), intent(out) :: PRES  (KA,IA,JA)       ! pressure [Pa]
    real(RP), intent(out) :: QV    (KA,IA,JA)       ! water vaper mixing ratio [kg/kg]
    real(RP), intent(out) :: QHYD  (KA,IA,JA,MP_QA) ! water mixing ratio [kg/kg] not warervaper
    real(RP), intent(out) :: QSAT  (KA,IA,JA)       ! saturation vaper mixing ratio
    real(RP), intent(out) :: rh    (KA,IA,JA)       ! relative humidity
    real(RP), intent(out) :: deltap(KA,IA,JA)       ! pressure difference  calculated by hydrostatic balance)
    real(RP), intent(out) :: deltaz(KA,IA,JA)       ! delta z (center level)
    real(RP), intent(out) :: deltax                 ! delta x

    real(RP) :: PSAT(KA,IA,JA) ! saturation vaper pressure
    real(RP) :: QDRY(KA,IA,JA) ! dry mixing ratio

    integer  :: k, i, j, iq
    !----------------------------------------------------------------------------

    ! calculate u(x-directin velocity ), v(y-direction velocity)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       u(k,i,j) = 0.5_RP * ( MOMX(k,i,j) + MOMX(k,i-1,j) ) / DENS(k,i,j)
       v(k,i,j) = 0.5_RP * ( MOMY(k,i,j) + MOMY(k,i,j-1) ) / DENS(k,i,j)
    enddo
    enddo
    enddo

    call THERMODYN_temp_pres( TEMP(:,:,:),  & ! [OUT]
                              PRES(:,:,:),  & ! [OUT]
                              DENS(:,:,:),  & ! [IN]
                              RHOT(:,:,:),  & ! [IN]
                              QTRC(:,:,:,:) ) ! [IN]

    ! calculate water vaper and relative humidity
    call THERMODYN_qd( QDRY(:,:,:), QTRC(:,:,:,:) )

    call SATURATION_psat_liq( PSAT(:,:,:), TEMP(:,:,:) )

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       QSAT(k,i,j) = 0.622_RP * PSAT(k,i,j) / ( PRES(k,i,j) - ( 1.0_RP-0.622_RP ) * PSAT(k,i,j) )
       QV  (k,i,j) = QTRC(k,i,j,I_QV) / QDRY(k,i,j)
       QV  (k,i,j) = max( 0.000001_RP, min( QSAT(k,i,j), QV(k,i,j) ) ) ! conpare QSAT and QV, guess lower limit
       rh  (k,i,j) = QV(k,i,j) / QSAT(k,i,j)
    enddo
    enddo
    enddo

    do iq = 1, MP_QA
    do j  = JS, JE
    do i  = IS, IE
    do k  = KS, KE
       QHYD(k,i,j,iq) = QTRC(k,i,j,I_MP2ALL(iq)) / QDRY(k,i,j)
    enddo
    enddo
    enddo
    enddo

    Z(:,:,:) = CZ(:,:,:) ! becouse scale_atmos_phy_cp interface ,not use scale_grid

    ! calculate delta P by hydrostatic balance
    ! deltap is the pressure interval between half levels(face levels) @ SCALE
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       deltap(k,i,j) = DENS(k,i,j) * GRAV * ( FZ(k+1,i,j) - FZ(k,i,j) ) ! rho*g*dz
    enddo
    enddo
    enddo

    ! deltaz is the interval of between model full levels(scalar point )
    deltaz(:,:,:) = 0.0_RP
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       deltaz(k,i,j) = CZ(k+1,i,j) - CZ(k,i,j)
    enddo
    enddo
    enddo
    deltaz(KE,:,:) = 0.0_RP

    deltax = sqrt( DX*DY )

    return
  end subroutine CP_kf_init

  !-----------------------------------------------------------------------------
  !> running mean vertical wind speed
  ! WRF comment out for W0
  !...TST IS THE NUMBER OF TIME STEPS IN 10 MINUTES...W0AVG IS CLOSE TO A
  !...RUNNING MEAN VERTICAL VELOCITY...NOTE THAT IF YOU CHANGE TST, IT WIL
  !...CHANGE THE FREQUENCY OF THE CONVECTIVE INTITIATION CHECK (SEE BELOW)
  !...NOTE THAT THE ORDERING OF VERTICAL LAYERS MUST BE REVERSED FOR W0AVG
  !...BECAUSE THE ORDERING IS REVERSED IN KFPARA...
  subroutine KF_wmean( &
       W0_avg, &
       DENS,   &
       MOMZ    )
    use scale_time , only :&
       TIME_DTSEC
    use scale_atmos_phy_cp_kf_sub, only: &
       KF_DT
    implicit none

    real(RP), intent(inout) :: W0_avg(KA,IA,JA)
    real(RP), intent(in)    :: DENS  (KA,IA,JA)
    real(RP), intent(in)    :: MOMZ  (KA,IA,JA)

    real(RP) :: W0
    real(RP) :: fact_W0_avg, fact_W0

    integer :: k, i, j
    !---------------------------------------------------------------------------

    if ( PARAM_ATMOS_PHY_CP_kf_wadapt ) then
       fact_W0_avg = 2.0_RP * max(KF_DT,TIME_DTSEC) - TIME_DTSEC
       fact_W0     = TIME_DTSEC
    else ! w_time is tuning parameter
       fact_W0_avg = real(PARAM_ATMOS_PHY_CP_kf_w_time,RP)
       fact_W0     = 1.0_RP
    endif

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       W0 = 0.5_RP * ( MOMZ(k,i,j) + MOMZ(k-1,i,j) ) / DENS(k,i,j)

       W0_avg(k,i,j) = ( W0_avg(k,i,j) * fact_W0_avg &
                       + W0            * fact_W0     ) / ( fact_W0_avg + fact_W0 )
    enddo
    enddo
    enddo

    return
  end subroutine KF_wmean

end module scale_atmos_phy_cp_kf
