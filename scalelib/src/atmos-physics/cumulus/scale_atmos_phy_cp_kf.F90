!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cumulus Parameterization
!!
!! @par Description
!!          Cumulus convection by Kain-Fritsch parameterization
!!          Reference: Kain and Fritsch(1990)
!!                     Kain (2004)
!!                     Narita and Ohmori (2007)
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2016-06-27 (S.Matsugishi) [new]
!!
!!
!! This file was originally copied from WRF.
!! The original file was published with the following notice.
!!
!! WRF was developed at the National Center for Atmospheric Research (NCAR) which is operated by the University Corporation for Atmospheric Research (UCAR). NCAR and UCAR make no proprietary claims, either statutory or otherwise, to this version and release of WRF and consider WRF to be in the public domain for use by any person or entity for any purpose without any fee or charge. UCAR requests that any WRF user include this notice on any partial or full copies of WRF. WRF is provided on an "AS IS" basis and any warranties, either express or implied, including but not limited to implied warranties of non-infringement, originality, merchantability and fitness for a particular purpose, are disclaimed. In no event shall UCAR be liable for any damages, whatsoever, whether direct, indirect, consequential or special, that arise out of or in connection with the access, use or performance of WRF, including infringement actions.
!!
!! WRF® is a registered trademark of the University Corporation for Atmospheric Research (UCAR).
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
  use scale_const, only: &
       TEM00 => CONST_TEM00

  use scale_tracer, only: QA
  !------------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_CP_kf_setup
  public :: ATMOS_PHY_CP_kf

  !-----------------------------------------------------------------------------
  !++ Public pparameters & variabeles
  real(RP), public,  save       :: KF_DT           = 5._RP    ! kf time scale [min!!!!!]

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: CP_kf_trigger, CP_kf_updraft, CP_kf_downdraft, CP_kf_compensational
  private :: calcexn
  private :: precipitation_OC1973
  private :: precipitation_Kessler
  private :: TPMIX2, DTFRZNEW, PROF5, TPMIX2DD, ENVIRTHT

  abstract interface
     subroutine precipitation( &
          QLIQ,QICE,WTW,DZ,BOTERM,ENTERM,QNEWLQ, &
          QNEWIC,QLQOUT,QICOUT,G)
       use scale_precision
       real(RP), INTENT(IN   )   :: G
       real(RP), INTENT(IN   )   :: DZ,BOTERM,ENTERM!,RATE to be local variablebles
       real(RP), INTENT(INOUT)   :: QLQOUT,QICOUT,WTW,QLIQ,QICE,QNEWLQ,QNEWIC
     end subroutine precipitation
  end interface

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  ! KF subroutine look up tables VV
  integer,  private,PARAMETER    :: KFNT=250,KFNP=220
  real(RP), private,SAVE         :: TTAB(KFNT,KFNP),QSTAB(KFNT,KFNP)
  real(RP), private,SAVE         :: THE0K(KFNP)
  real(RP), private,SAVE         :: ALU(200)
  real(RP), private,SAVE         :: RDPR,RDTHK,PLUTOP
  !^^^^
  real(RP), private,SAVE         :: GdCP ! GRAV/CP_dry
  !
  ! ALIQ saturrate watervape
  ! BLIQ is WRF SVP2
  ! DLIQ is WRF SVP3
  real(RP),private,parameter     :: ALIQ=6.112e2_RP ! Saturate pressure of water vapor [Pa]
  real(RP),private,parameter     :: BLIQ=17.68_RP   ! emanuel 1994 (4.6.2) 17.67
  real(RP),private,parameter     :: CLIQ=BLIQ*TEM00 ! convert degree to kelvin @ dew point temperature
  real(RP),private,parameter     :: DLIQ=29.65_RP   ! 273.15 - 243.5
  real(RP),private,parameter     :: XLV1=2370._RP,XLV0=3.15e6_RP
  ! Naming list tuning parameters
  ! RATE is used subroutine precipitation_OC1973
  real(RP),private, save       :: RATE            = 0.03_RP  ! ratio of cloud water and precipitation (Ogura and Cho 1973)
  integer ,private, save       :: TRIGGER         = 3        ! triger select will be modifid
  logical ,private, save       :: FLAG_QS         = .true.   ! FLAG_OS:  qs is allowed or not
  logical ,private, save       :: FLAG_QI         = .true.   ! FLAG_QI:  qi is allowe or not
  real(RP),private, save       :: DELCAPE         = 0.1_RP   ! cape decleace rate
  real(RP),private, save       :: DEEPLIFETIME    = 1800._RP ! minimum lifetimescale of deep convection
  real(RP),private, save       :: SHALLOWLIFETIME = 2400._RP ! shallow convection liftime 
  real(RP),private, save       :: DEPTH_USL       = 300._RP  ! depth of updraft source layer  [hPa]!!
  logical ,private, save       :: WARMRAIN        = .false.  ! QA<=3 then ??
  logical ,private, save       :: KF_LOG          = .false.  ! KF infomation output to log file(not ERROR messeage)
  real(RP),private, save       :: kf_threshold    = 1.e-3_RP ! kessler type autoconversion rate
  integer ,private, save       :: stepkf                     !! triger select will be modifid
  integer ,private, save       :: KF_prec         = 1        ! precipitation select 1. Ogura and Cho (1973), 2. Kessler
  procedure(precipitation), pointer,private :: p_precipitation => NULL()

  real(RP), private, allocatable :: lifetime  (:,:) ! convectime lifetime [s]
  integer , private, allocatable :: I_convflag(:,:) ! convection type 0:deep convection 1:shallow convection 2: no convection

  real(RP), private, allocatable :: deltaz (:,:,:) ! height interval (center level) [m]
  real(RP), private, allocatable :: Z      (:,:,:) ! centerlevel real height [m]
  real(RP)                       :: deltax         ! delta x [m]


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
    use scale_grid_real, only: &
       CZ => REAL_CZ
    use scale_grid,only: &
       DX => DX, &
       DY => DY
    implicit none

    character(len=*), intent(in) :: CP_TYPE

    ! tunning parameters, original parameter set is from KF2004 and NO2007
    real(RP) :: PARAM_ATMOS_PHY_CP_kf_rate      =   0.03_RP ! ratio of cloud water and precipitation (Ogura and Cho 1973)
    integer  :: PARAM_ATMOS_PHY_CP_kf_trigger   = 1         ! trigger function type 1:KF2004 3:NO2007
    logical  :: PARAM_ATMOS_PHY_CP_kf_qs        = .true.    ! qs is allowed?
    logical  :: PARAM_ATMOS_PHY_CP_kf_qi        = .true.    ! qi is allowed?
    real(DP) :: PARAM_ATMOS_PHY_CP_kf_dt        =    5.0_DP ! KF convection check time interval [min]
    real(RP) :: PARAM_ATMOS_PHY_CP_kf_dlcape    =    0.1_RP ! cape decleace rate
    real(RP) :: PARAM_ATMOS_PHY_CP_kf_dlifetime = 1800.0_RP ! minimum lifetime scale of deep convection [sec]
    real(RP) :: PARAM_ATMOS_PHY_CP_kf_slifetime = 2400.0_RP ! lifetime of shallow convection [sec]
    real(RP) :: PARAM_ATMOS_PHY_CP_kf_DEPTH_USL =  300.0_RP ! depth of updraft source layer [hPa]
    integer  :: PARAM_ATMOS_PHY_CP_kf_prec      = 1         ! precipitation type 1:Ogura-Cho(1973) 2:Kessler
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
       PARAM_ATMOS_PHY_CP_kf_prec,      &
       PARAM_ATMOS_PHY_CP_kf_thres,     &
       PARAM_ATMOS_PHY_CP_kf_warmrain,  &
       PARAM_ATMOS_PHY_CP_kf_LOG,       &
       PARAM_ATMOS_PHY_CP_kf_wadapt, &
       PARAM_ATMOS_PHY_CP_kf_w_time

    integer :: k, i, j
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
                      PARAM_ATMOS_PHY_CP_kf_prec,      &
                      PARAM_ATMOS_PHY_CP_kf_thres,     &
                      PARAM_ATMOS_PHY_CP_kf_warmrain,  &
                      PARAM_ATMOS_PHY_CP_kf_LOG ,      &
                      TIME_DSTEP_KF                    )

    ! output parameter lists
    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) "*** Interval for check [step]                       : ", TIME_DSTEP_KF
    if( IO_L ) write(IO_FID_LOG,*) "*** Ogura-Cho condense material convert rate        : ", PARAM_ATMOS_PHY_CP_kf_rate
    if( IO_L ) write(IO_FID_LOG,*) "*** Trigger function type, 1:KF2004 3:NO2007        : ", PARAM_ATMOS_PHY_CP_kf_trigger
    if( IO_L ) write(IO_FID_LOG,*) "*** Exist qi?                                       : ", PARAM_ATMOS_PHY_CP_kf_qi
    if( IO_L ) write(IO_FID_LOG,*) "*** Exist qs?                                       : ", PARAM_ATMOS_PHY_CP_kf_qs
    if( IO_L ) write(IO_FID_LOG,*) "*** CAPE decrease rate                              : ", PARAM_ATMOS_PHY_CP_kf_dlcape
    if( IO_L ) write(IO_FID_LOG,*) "*** Minimum lifetime scale of deep convection [sec] : ", PARAM_ATMOS_PHY_CP_kf_dlifetime
    if( IO_L ) write(IO_FID_LOG,*) "*** Lifetime of shallow convection            [sec] : ", PARAM_ATMOS_PHY_CP_kf_slifetime
    if( IO_L ) write(IO_FID_LOG,*) "*** Updraft souce layer depth                 [hPa] : ", PARAM_ATMOS_PHY_CP_kf_DEPTH_USL
    if( IO_L ) write(IO_FID_LOG,*) "*** Precipitation type 1:Ogura-Cho(1973) 2:Kessler  : ", PARAM_ATMOS_PHY_CP_kf_prec
    if( IO_L ) write(IO_FID_LOG,*) "*** Kessler type precipitation's threshold          : ", PARAM_ATMOS_PHY_CP_kf_thres
    if( IO_L ) write(IO_FID_LOG,*) "*** Warm rain?                                      : ", PARAM_ATMOS_PHY_CP_kf_warmrain
    if( IO_L ) write(IO_FID_LOG,*) "*** Use running mean of w in adaptive timestep?     : ", PARAM_ATMOS_PHY_CP_kf_wadapt
    if( IO_L ) write(IO_FID_LOG,*) "*** Fixed time scale for running mean of w          : ", PARAM_ATMOS_PHY_CP_kf_w_time
    if( IO_L ) write(IO_FID_LOG,*) "*** Output log?                                     : ", PARAM_ATMOS_PHY_CP_kf_LOG

    ! output variables
    allocate( lifetime  (IA,JA) )
    allocate( I_convflag(IA,JA) )
    lifetime  (:,:) = 0.0_RP
    I_convflag(:,:) = 2

    allocate( Z(KA,IA,JA) )
    Z(:,:,:) = CZ(:,:,:) ! becouse scale_atmos_phy_cp interface ,not use scale_grid

    allocate( deltaz(KA,IA,JA) )
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

    real(RP) :: cldfrac(KA,2)           ! 1 shallow , 2 deep

    !---------------------------------------------------------------------------

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

       call PROF_rapstart('CP_kf', 3)

       ! calc cumulus convection
       call CP_kf_main( DENS         (:,:,:),   & ! [IN]
                        MOMZ         (:,:,:),   & ! [IN]
                        MOMX         (:,:,:),   & ! [IN]
                        MOMY         (:,:,:),   & ! [IN]
                        RHOT         (:,:,:),   & ! [IN]
                        QTRC         (:,:,:,:), & ! [IN]
                        w0avg        (:,:,:),   & ! [IN]
                        nca          (:,:),     & ! [INOUT]
                        DENS_t_CP    (:,:,:),   & ! [OUT]
                        RHOT_t_CP    (:,:,:),   & ! [OUT]
                        RHOQ_t_CP    (:,:,:,:), & ! [OUT]
                        SFLX_convrain(:,:),     & ! [OUT]
                        cldfrac_sh   (:,:,:),   & ! [OUT]
                        cldfrac_dp   (:,:,:),   & ! [OUT]
                        lifetime     (:,:),     & ! [OUT]
                        cloudtop     (:,:),     & ! [OUT]
                        cloudbase    (:,:),     & ! [OUT]
                        I_convflag   (:,:)      ) ! [OUT]

       call PROF_rapend('CP_kf', 3)
    endif

    call HIST_in( lifetime(:,:),            'KF_LIFETIME', 'lifetime of KF scheme', 's' )
    call HIST_in( real(I_convflag(:,:),RP), 'KF_CONVFLAG', 'CONVECTION FLAG',       ''  )

    return
  end subroutine ATMOS_PHY_CP_kf

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

  !-----------------------------------------------------------------------------
  subroutine CP_kf_param( & ! set kf tuning parametres
       ![IN]
       RATE_in,            &
       TRIGGER_in,         & ! INOUT
       FLAG_QS_in,         &
       FLAG_QI_in,         &
       KF_DT_in,           &
       DELCAPE_in ,        &
       DEEPLIFETIME_in,    &
       SHALLOWLIFETIME_in, &
       DEPTH_USL_in,       &
       KF_prec_in,         &
       KF_threshold_in,    &
       WARMRAIN_in,        &
       KF_LOG_in,          &
       stepkf_in     )
    use scale_process, only: &
         PRC_MPIstop
    implicit none
    real(RP),intent(in) :: RATE_in
    integer, intent(inout) :: TRIGGER_in
    logical, intent(in) :: FLAG_QS_in,FLAG_QI_in
    real(RP),intent(in) :: KF_DT_in
    real(RP),intent(in) :: DELCAPE_in
    real(RP),intent(in) :: DEEPLIFETIME_in
    real(RP),intent(in) :: SHALLOWLIFETIME_in
    real(RP),intent(in) :: DEPTH_USL_in
    integer, intent(in) :: KF_prec_in
    real(RP),intent(in) :: KF_threshold_in
    logical, intent(in) :: WARMRAIN_in
    logical, intent(in) :: KF_LOG_in
    integer, intent(in) :: stepkf_in
    !
    RATE            = RATE_in
    ! TRIGGER must be 1 or 3
    if (TRIGGER_in /= 1 .and. TRIGGER_in /=3) then
       if (IO_L) write(IO_FID_LOG,*) "TRIGGER must be 1 or 3 but now :",TRIGGER_in
       if (IO_L) write(IO_FID_LOG,*) "CHAGNGE ",TRIGGER_in," to 3"
       TRIGGER_in = 3
    end if
    TRIGGER         = TRIGGER_in
    FLAG_QS         = FLAG_QS_in
    FLAG_QI         = FLAG_QI_in
    KF_DT           = KF_DT_in
    DELCAPE         = DELCAPE_in
    DEEPLIFETIME    = DEEPLIFETIME_in
    SHALLOWLIFETIME = SHALLOWLIFETIME_in
    DEPTH_USL       = DEPTH_USL_in
    WARMRAIN        = WARMRAIN_in
    KF_prec         = KF_prec_in
    KF_threshold    = KF_threshold_in
    KF_LOG          = KF_LOG_in
    stepkf          = stepkf_in
    if (KF_prec == 1) then
       p_precipitation => precipitation_OC1973 ! Ogura and Cho (1973)
    elseif( KF_prec == 2) then
       p_precipitation => precipitation_Kessler ! Kessler type
    else
       write(*,*) 'xxx ERROR at KF namelist'
       write(*,*) 'KF_prec must be 1 or 2'
       call PRC_MPIstop
    end if
    return
  end subroutine CP_kf_param

  subroutine CP_kf_main ( & ! main loop 
       ![IN]
       dens,        &
       MOMZ,        &
       MOMX,        &
       MOMY,        &
       RHOT,        &
       QTRC,        &
       w0avg,       &
       ! [INOUT]
       nca,         &
       ! [OUT]
       DENS_t_CP,   &
       DT_RHOT,     &
       DT_RHOQ,     &
       rainrate_cp, &
       cldfrac_sh,  &
       cldfrac_dp,  &
       timecp,      &
       cloudtop,    &
       zlcl,        &
       I_convflag   )
    use scale_precision
    use scale_grid_index
    use scale_tracer
    use scale_const, only: &
       GRAV => CONST_GRAV, &
       R    => CONST_Rdry
    use scale_atmos_phy_mp, only: &
       QA_MP, &
       QS_MP
    use scale_atmos_hydrometer, only: &
       I_QV, &
       I_QC, &
       I_QR, &
       I_QI, &
       I_QS, &
       I_QG
    use scale_time , only :&
       dt => TIME_DTSEC_ATMOS_PHY_CP
    use scale_grid_real, only: &
       FZ => REAL_FZ
    use scale_atmos_thermodyn, only: &
       THERMODYN_temp_pres   => ATMOS_THERMODYN_temp_pres,   &
       THERMODYN_rhoe        => ATMOS_THERMODYN_rhoe,        &
       THERMODYN_temp_pres_E => ATMOS_THERMODYN_temp_pres_E, &
       THERMODYN_qd          => ATMOS_THERMODYN_qd,          &
       THERMODYN_pott        => ATMOS_THERMODYN_pott
    use scale_atmos_saturation ,only :&
       SATURATION_psat_liq => ATMOS_SATURATION_psat_liq
    implicit none
    ![IN]
    real(RP),intent(in) :: DENS(KA,IA,JA)          ! density [kg/m**3]
    real(RP),intent(in) :: MOMZ(KA,IA,JA)          ! momentum
    real(RP),intent(in) :: MOMX(KA,IA,JA)          ! momentum
    real(RP),intent(in) :: MOMY(KA,IA,JA)          ! momentum
    real(RP),intent(in) :: RHOT(KA,IA,JA)          ! density*PT
    real(RP),intent(in) :: QTRC(KA,IA,JA,QA)       ! raito of water elements
    real(RP),intent(in) :: w0avg(KA,IA,JA)         ! running mean vertical wind velocity [m/s]
    ! [INOUT]
    real(RP),intent(inout) :: nca(IA,JA)           ! num of step convection active [step]
    real(RP),intent(out)   :: DENS_t_CP(KA,IA,JA)  ! dens/dt
    real(RP),intent(out)   :: DT_RHOT(KA,IA,JA)    ! dens*PT/dt
    real(RP),intent(out)   :: DT_RHOQ(KA,IA,JA,QA) ! dens*q/dt
    real(RP),intent(out)   :: rainrate_cp(IA,JA)   ! rain  PPTFLX(prcp_flux)/deltax**2*dt ! convective rain
    real(RP),intent(out)   :: cldfrac_sh(KA,IA,JA) ! cloud fraction
    real(RP),intent(out)   :: cldfrac_dp(KA,IA,JA) ! cloud fraction
    real(RP),intent(out)   :: timecp(IA,JA)        ! timescale of cumulus parameterization
    real(RP),intent(out)   :: cloudtop(IA,JA)      ! cloud top height
    real(RP),intent(out)   :: zlcl(IA,JA)          ! hight of lcl cloud bottom height[m]
    integer, intent(out)   :: I_convflag(IA,JA)    ! convection type
    !> I_convflag = 0  ==> deep convection
    !>            = 1  ==> shallow convection
    !>            = 2  ==> NONE !!

    ! [Internal Work]
    real(RP) :: u     (KA)       ! x-direction wind velocity [m/s]
    real(RP) :: v     (KA)       ! y-direction wind velocity [m/s]
    real(RP) :: temp  (KA)       ! temperature [K]
    real(RP) :: pres  (KA)       ! pressure [Pa]
    real(RP) :: qv    (KA)       ! water vapor mixing ratio[kg/kg]
    real(RP) :: QDRY  (KA)       ! ratio of dry air
    real(RP) :: qes   (KA)       ! saturate vapor [kg/kg]
    real(RP) :: PSAT  (KA)       ! saturation vaper pressure
    real(RP) :: QSAT  (KA)       ! saturate water vaper mixing ratio [kg/kg]
    real(RP) :: rh    (KA)       ! saturate vapor [%]
    real(RP) :: deltap(KA)       ! delta Pressure [Pa]

    real(RP) :: q_hyd(KA,QA_MP)  ! water mixing ratio [kg/kg]
    real(RP) :: dens_nw(KA)      ! density [kg/m**3]
    integer  :: nic
    real(RP) :: time_advec       ! advection timescale
    real(RP) :: umf(KA)          ! Updraft Mass Flux
    real(RP) :: umflcl           ! Updraft Mass Flux
    real(RP) :: upent(KA)        ! Updraft Mass flux ENTrainment
    real(RP) :: updet(KA)        ! Updraft Mass flux DETrainment
    ! updraft value
    real(RP) :: temp_u(KA)       ! updraft temperature [K]
    real(RP) :: tempv(KA)        ! vertual temperature [K]
    real(RP) :: qv_u(KA)         ! updraft qv
    real(RP) :: cape             ! cape
    ! water
    real(RP) :: qc(KA)           ! cloud water
    real(RP) :: qi(KA)           ! cloud ice
    real(RP) :: qvdet(KA)        ! vapor detrainment
    real(RP) :: qcdet(KA)        ! liquit detrainment
    real(RP) :: qidet(KA)        ! ice detrainment
    real(RP) :: totalprcp        ! total precipitation
    real(RP) :: flux_qr(KA)      ! rain flux
    real(RP) :: flux_qs(KA)      ! snow flux
    ! upward theta
    real(RP) :: theta_eu(KA)     ! updraft equivalent  PT
    real(RP) :: theta_ee(KA)     ! environment equivalent PT
    ! convection type
    ! index valuables
    integer  :: k_lcl            ! index of LCL layer
    integer  :: k_top            ! index of cloud top hight
    integer  :: k_ml             ! index of melt layer (temp < tem00)
    integer  :: k_lc,k_let,k_pbl ! indexs
    real(RP) :: zmix             ! usl(up source layer) layer depth [m]
    real(RP) :: presmix          ! usl layer depth [Pa]
    real(RP) :: umfnewdold(KA)   ! umfnew/umfold
    real(RP) :: dpthmx           ! max depth of pressure
    ! move internalwark
    real(RP) :: ems(KA)          ! box weight[kg]
    real(RP) :: emsd(KA)         ! 1/ems
    ! [OUT]@downdraft
    real(RP) :: wspd(3)          ! horizontal wind speed 1 k_lcl, 2 k_z5,3 k_top
    real(RP) :: dmf(KA)          ! downdraft massflux
    real(RP) :: downent(KA)      ! downdraft entrainment
    real(RP) :: downdet(KA)      ! downdraft detrainment
    real(RP) :: theta_d(KA)      ! downdraft theta
    real(RP) :: qv_d(KA)         ! water vapor of downdraft
    real(RP) :: prcp_flux        ! precpitation flux
    real(RP) :: tder             ! temperature change from evap
    real(RP) :: CPR              ! all precipitation  before consider cloud bottom evaporation
    integer  :: k_lfs            ! LFS layer index
    ![OUT] @ compensasional subsidence
    real(RP) :: temp_g(KA)       ! temporarly temperature -> after iteration then new variable
    real(RP) :: qv_g(KA)         ! temporarly vapor -> after iteration then new variable
    real(RP) :: qc_nw(KA)        ! new cloud water mixing ratio [kg/kg]
    real(RP) :: qi_nw(KA)        ! new cloud ice  mixing ratio  [kg/kg]
    real(RP) :: qr_nw(KA)        ! new rain water mixing ratio  [kg/kg]
    real(RP) :: qs_nw(KA)        ! new snow water mixing ratio  [kg/kg]
    ! new variable
    real(RP) :: rhot_nw(KA)      ! rho*PT
    real(RP) :: qtrc_nw(KA,QA)   ! qv,qc,qr,qi,qs (qg not change)
    real(RP) :: pott_nw(KA)      ! new PT
    !
    real(RP) :: cldfrac_KF(KA,2) ! cloud fraction
    !
    real(RP) :: qd(KA) ! dry mixing ratio
    ! do loop variable
    integer  :: k, i, j, iq, iqa, ii


    do j = JS, JE
    do i = IS, IE

       nca(i,j) = nca(i,j) - real(TIME_DSTEP_KF,RP) * dt

       ! check convection
       if ( nca(i,j) .ge. 0.5_DP * dt ) cycle


       ! convert variables

       ! calculate u(x-directin velocity ), v(y-direction velocity)
       do k = KS, KE
          u(k) = 0.5_RP * ( MOMX(k,i,j) + MOMX(k,i-1,j) ) / DENS(k,i,j)
          v(k) = 0.5_RP * ( MOMY(k,i,j) + MOMY(k,i,j-1) ) / DENS(k,i,j)
       enddo

       do k = KS, KE
          call THERMODYN_temp_pres( TEMP(k),       & ! [OUT]
                                    PRES(k),       & ! [OUT]
                                    DENS(k,i,j),   & ! [IN]
                                    RHOT(k,i,j),   & ! [IN]
                                    QTRC(k,i,j,:), & ! [IN]
                                    TRACER_CV(:),  & ! [IN]
                                    TRACER_R(:),   & ! [IN]
                                    TRACER_MASS(:) ) ! [IN]

          ! calculate water vaper and relative humidity
          call THERMODYN_qd( QDRY(k), QTRC(k,i,j,:), TRACER_MASS(:) )

          call SATURATION_psat_liq( PSAT(k), TEMP(k) )

          QSAT(k) = 0.622_RP * PSAT(k) / ( PRES(k) - ( 1.0_RP-0.622_RP ) * PSAT(k) )
          QV  (k) = QTRC(k,i,j,I_QV) / QDRY(k)
          QV  (k) = max( 0.000001_RP, min( QSAT(k), QV(k) ) ) ! conpare QSAT and QV, guess lower limit
          rh  (k) = QV(k) / QSAT(k)
       enddo

       ! calculate delta P by hydrostatic balance
       ! deltap is the pressure interval between half levels(face levels) @ SCALE
       do k = KS, KE
          deltap(k) = DENS(k,i,j) * GRAV * ( FZ(k+1,i,j) - FZ(k,i,j) ) ! rho*g*dz
       enddo


       DENS_t_CP(KS:KE,i,j) = 0.0_RP
       DT_RHOT  (KS:KE,i,j) = 0.0_RP
       do iq = 1, QA
          DT_RHOQ(KS:KE,i,j,iq) = 0.0_RP
       end do
       cldfrac_KF (KS:KE,:) = 0.0_RP
       rainrate_cp(i,j) = 0.0_RP
       timecp     (i,j) = 0.0_RP
       cloudtop   (i,j) = 0.0_RP
       zlcl       (i,j) = 0.0_RP

       do iq = 1, QA_MP-1
          iqa = iq + QS_MP
       do k  = KS, KE
          q_hyd(k,iq) = QTRC(k,i,j,iqa) / QDRY(k)
       enddo
       enddo


       call CP_kf_trigger ( &
            ! [IN]
            deltaz(:,i,j), Z(:,i,j), & ! deltaz and height[m]
            qv    (:),       & ! water vapor mixing ratio [kg/kg]
            qes   (:),       & ! saturation water vapor mixing ratio
            pres  (:),       & ! pressure [Pa]
            deltap(:),       & ! interval of pressure [Pa]
            temp  (:),       & ! temperature [K]
            w0avg (:,i,j),   & ! average w
            ! [OUT]
            I_convflag(i,j), & ! convection flag
            cloudtop  (i,j), & ! cloud top height [m]
            temp_u    (:),   & ! updraft temperature
            tempv     (:),   & ! virtual temperature
            ! water values
            qv_u      (:),   & ! updraft water vapor mixing ratio
            qc        (:),   & ! cloud water mixing ratio
            qi        (:),   & ! ice mixing ratio
            qvdet     (:),   & ! updraft detrain vapor
            qcdet     (:),   & ! updraft detrain water
            qidet     (:),   & ! updraft detrain ice
            flux_qr   (:),   & ! rain flux
            flux_qs   (:),   & ! snow flux
            totalprcp,       & ! total precipitation(rain and snow ) fall
            ! thermodynamics values
            theta_eu  (:),   & ! updraft equivalent PT
            theta_ee  (:),   & ! environment equivalent PT
            cape,            & ! cape
            ! massflux values
            umf       (:),   & ! updraft mass flux
            umflcl,          & ! updraft mass flux@lcl
            upent     (:),   & ! updraft entrainment
            updet     (:),   & ! updraft detrainment
            ! layer nambers
            k_lcl,           & ! index of LCL layer
            k_lc,            & ! index of USL layer botom
            k_pbl,           & ! index of USL layer top
            k_top,           & ! index of cloud top
            k_let,           & ! index of below detrain layer
            k_ml,            & ! index of temprerure < tem00
            ! pbl layer values
            presmix,         & ! USL layer pressure [Pa]
            dpthmx,          & ! USL layer depth ( [Pa])
            zlcl      (i,j), & ! lcl hight
            zmix,            & ! USL layer depth
            umfnewdold(:)    )
       !------------------------------------------------------------------------
       ! calc ems(box weight[kg])
       !------------------------------------------------------------------------
       if(I_convflag(i,j) /= 2) then
          ems(k_top+1:KE) = 0._RP
          emsd(k_top+1:KE) = 0._RP
          do k = KS, k_top
             ems(k) = deltap(k)*deltax*deltax/GRAV
             emsd(k) = 1._RP/ems(k)
             umfnewdold(k) = 1._RP/umfnewdold(k)
          end do
       end if
       !------------------------------------------------------------------------
       call  CP_kf_downdraft ( &
            ! [IN]
            deltaz(:,i,j), Z(:,i,j) ,& ! deltaz and height [m]
            u(:)             ,& ! u-wind m/s
            v(:)             ,& ! v-wind m/s
            zlcl(i,j)        ,& ! lcl height [m]
            rh(:)            ,& ! relative humidity
            deltap(:)        ,& ! interval of pressure [Pa]
            pres(:)          ,& ! pressure [Pa]
            qv(:)            ,& ! water vapor mixing ratio [kg/kg]
            ems(:)           ,& ! ems(box weight[kg])
            emsd(:)          ,& ! 1/ems
            theta_ee(:)      ,& ! environment theta_E
            umf(:)           ,& ! updraft mass flux
            totalprcp        ,& ! total precipitation
            flux_qs(:)       ,& ! snow flux
            tempv(:)         ,& ! vertual temperature [K]
            I_convflag(i,j)       ,& ! convection type
            ! layer index 
            k_lcl            ,& ! index of LCL layer
            k_ml             ,& ! index of melt layer 
            k_top            ,& ! index of cloud top
            k_pbl            ,& ! index of USL layer top
            k_let            ,& ! index of below detrain layer
            k_lc             ,& ! index of USL layer botom
            ! [out]
            wspd(:)          ,& ! wind speed 1 k_lcl, 2 k_z5,3 k_top
            dmf(:)           ,& ! downward mass flux
            downent(:)       ,& ! downdraft entrainment
            downdet(:)       ,& ! downdraft detrainment
            theta_d(:)       ,& ! downdraft theta
            qv_d(:)          ,& ! downdraft water vaper mixng ratio
            prcp_flux        ,& ! precipitateion
            k_lfs            ,& ! index of LFS layer
            CPR              ,& ! all precipitation  before consider cloud bottom evaporation
            tder)
       !------------------------------------------------------------------------
       call CP_kf_compensational ( &
            ![IN]
            deltaz(:,i,j),Z(:,i,j) ,& ! deltaz and height [m]
            pres(:),          & ! pressure [Pa]
            deltap(:),        & ! deltap [Pa] 
            temp(:),          & ! temperature [K]
            qv(:),            & ! water vapor mixing ratio
            ! form kf_trigger
            ems(:),           & ! box weight[kg]
            emsd(:),          & ! 1/ems
            presmix,          & ! usl layer pressure depth
            zmix,             & ! usl layer pressure depth
            dpthmx,           & ! usl layer depth
            cape,             & ! CAPE
            temp_u(:),        & ! updraft temperature
            qvdet(:),         & ! detrainment water vaper mixing ratio
            umflcl,           & ! umf @LCL
            umf(:),           & ! upward mass flux (updraft)
            upent(:),         & ! updraft entrainment
            updet(:),         & ! updraft detrainment
            qcdet(:),         & ! updraft detrainment qc
            qidet(:),         & ! updraft detrainment qi
            qc(:),            & ! cloud water mixing ratio
            qi(:),            & ! cloud ice   mixing ratio
            flux_qr(:),       & ! rain flux
            flux_qs(:),       & ! snow flux
            umfnewdold(:),    & ! umfnew/umfold ratio
            ! from kf_downdraft
            wspd(:),          & ! wind speed 1. 2 is used
            dmf(:),           & ! downward  mass flux (downdraft)
            downent(:),       & ! downdraft entrainment 
            downdet(:),       & ! downdraft detrainment
            qv_d(:),          & ! downdraft detrainment qv
            theta_d(:),       & ! potential temperature @downdraft
            prcp_flux,        & ! precipitation
            tder,             & ! evapolation
            cpr,              & ! all precipitation  before consider cloud bottom evaporation
            I_convflag(i,j),  & ! intent inout convective flag
            ! layer index 
            ! from kf_triger
            k_top,            & ! index of cloud top hight
            k_lcl,            & ! index of LCL layer 
            k_lc ,            & ! index of LCL layer botoom
            k_pbl,            & ! index of USL layer top
            k_ml,             & ! index of melt layer (temp < tem00)
            ! from kf_downdraft
            k_lfs,            & ! index of LFS layer
            ![OUT] new valuavles after timestep
            temp_g(:),        & ! update temperature [K]
            qv_g(:),          & ! update water vapor mixing ratio
            qc_nw(:),         & ! update cloud water mixing ratio
            qi_nw(:),         & ! update cloud ice   mixing ratio
            qr_nw(:),         & ! update rain  water mixing ratio
            qs_nw(:),         & ! update snow  water mixing ratio
            rainrate_cp(i,j), & ! rain  PPTFLX(prcp_flux)/deltax**2*dt ! convective rain
            nic,              & ! (timescale of cumulus parameterization)/dt integer
            cldfrac_KF,       & ! cloud fraction
            timecp(i,j),      & ! timescale of cumulus parameterization
            time_advec) ! advection timescale
       
       !------------------------------------------------------------------------
       ! compute tendencys
       !------------------------------------------------------------------------

       if(I_convflag(i,j) == 2) then ! no convection
          rainrate_cp(i,j)    = 0.0_RP
          cldfrac_KF(KS:KE,:) = 0.0_RP
          timecp(i,j)         = 0.0_RP
          cloudtop(i,j)       = 0.0_RP
          zlcl(i,j)           = 0.0_RP
       else ! convection allowed I_convflag=0 or 1
          ! chek
          !
          !...FEEDBACK TO RESOLVABLE SCALE TENDENCIES.
          !
          !...IF THE ADVECTIVE TIME PERIOD (time_advec) IS LESS THAN SPECIFIED MINIMUM
          !...TIMECP, ALLOW FEEDBACK TO OCCUR ONLY DURING TADVEC...
          !
          if (I_convflag(i,j) == 0) then ! deep
             if (time_advec < timecp(i,j)) nic=nint(time_advec/dt)
             nca(i,j) = real(nic,RP)*dt ! convection feed back act this time span
          elseif (I_convflag(i,j) == 1) then ! shallow
             timecp(i,j) = SHALLOWLIFETIME
             nca   (i,j) = KF_DT*60._RP ! convection feed back act this time span
          end if
          ! update qd
          if ( I_QC>0 ) q_hyd(KS:KE,I_QC-QS_MP) = qc_nw(KS:KE) + q_hyd(KS:KE,I_QC-QS_MP)
          if ( I_QR>0 ) q_hyd(KS:KE,I_QR-QS_MP) = qr_nw(KS:KE) + q_hyd(KS:KE,I_QR-QS_MP)
          if ( I_QI>0 ) q_hyd(KS:KE,I_QI-QS_MP) = qi_nw(KS:KE) + q_hyd(KS:KE,I_QI-QS_MP)
          if ( I_QS>0 ) q_hyd(KS:KE,I_QS-QS_MP) = qs_nw(KS:KE) + q_hyd(KS:KE,I_QS-QS_MP)
          if ( I_QG>0 ) q_hyd(KS:KE,I_QG-QS_MP) = qs_nw(KS:KE) + q_hyd(KS:KE,I_QG-QS_MP)
          qd(KS:KE) = 1.0_RP / ( 1.0_RP + qv_g(KS:KE) + sum(q_hyd(KS:KE,:),2)) ! new qdry
          ! new qtrc
          qtrc_nw(KS:KE,I_QV) = qv_g(KS:KE)*qd(KS:kE)
          do iq = 1, QA_MP-1
             qtrc_nw(KS:KE,QS_MP+iq) = q_hyd(KS:KE,iq) * qd(KS:KE)
          end do
          ! new density
          dens_nw(KS:KE) = dens(KS:KE,i,j)
          do ii = 1, QA_MP
             iq = QS_MP + ii - 1
             dens_nw(KS:KE) = dens_nw(KS:KE) + ( qtrc_nw(KS:KE,iq) - QTRC(KS:KE,i,j,iq) )*dens(KS:KE,i,j)
          end do
          ! dens_(n+1) = dens_(n) + tend_q*dens_(n)

          ! calc new potential temperature
          call THERMODYN_pott(&
               pott_nw(:), &
               temp_g(:), pres(:), qtrc_nw(:,:), &
               TRACER_CV(:), TRACER_R(:), TRACER_MASS(:) )
          ! update rhot
          rhot_nw(KS:KE) = dens_nw(KS:KE)*pott_nw(KS:KE)
          ! calc tendency
          DENS_t_CP(KS:KE,i,j) = (dens_nw(KS:KE) - dens(KS:KE,i,j))/timecp(i,j)
          DT_RHOT(KS:KE,i,j)   = (rhot_nw(KS:KE) - RHOT(KS:KE,i,j))/timecp(i,j)
          do ii = 1, QA_MP
             iq = QS_MP + ii - 1
             DT_RHOQ(KS:KE,i,j,iq) = ( dens_nw(KS:KE) * qtrc_nw(KS:KE,iq) - DENS(KS:KE,i,j) * QTRC(KS:KE,i,j,iq) ) &
                  /timecp(i,j)
          end do
          ! if noconvection then nca is same value before call. nca only modifyed convectioned
       end if

       cldfrac_sh(KS:KE,i,j) = cldfrac_KF(KS:KE,1)
       cldfrac_dp(KS:KE,i,j) = cldfrac_KF(KS:KE,2)

    end do
    end do

    return
  end subroutine CP_kf_main
  !------------------------------------------------------------------------------
  ! calc trigger function and upward mass flux
  !------------------------------------------------------------------------------  
  subroutine CP_kf_trigger ( & !> 1d model ! triger function
       ! [IN]
       dz_kf,z_kf ,& ! deltaz and height [m]
       qv         ,& ! water vapor mixing ratio [kg/kg]
       qes        ,& ! saturated vapor [kg/kg]
       pres       ,& ! pressure [Pa]
       deltap     ,& ! interval of pressure [Pa]
       temp       ,& ! temperature [K]
       w0avg      ,& ! running mean w
       ! [OUT]
       I_convflag ,& ! convection flag
       cloudtop   ,& ! cloud top height [m]
       temp_u     ,& ! updraft temperature
       tempv      ,& ! virtual temperature
       ! water variable
       qv_u       ,& ! updraft water vapor mixing ratio
       qc         ,& ! cloud water mixing ratio
       qi         ,& ! ice mixingratio
       qvdet      ,& ! updraft detrain vapor
       qcdet      ,& ! updraft detrain water
       qidet      ,& ! updraft detrain ice
       flux_qr    ,& ! rain flux
       flux_qs    ,& ! snow flux
       totalprcp  ,& ! total precipitation (rain and snow fall)
       ! thermodynamics variables
       theta_eu   ,& ! updraft equivalent theta
       theta_ee   ,& ! environment equivalent theta
       cape       ,& ! cape
       ! massflux variables
       umf        ,& ! updraft mass flux
       umflcl     ,& ! updraft mass flux@LCL
       upent      ,& ! updraft entrainment
       updet      ,& ! updraft detrainment
       ! layer indexs
       k_lcl      ,& ! index of LCL layer
       k_lc       ,& ! index of USL layer botom
       k_pbl      ,& ! index of USL layer top
       k_top      ,& ! index of cloud top
       k_let      ,& ! index of below detrain layer
       k_ml       ,& ! index of temprerure < tem00 (melting layer)
       ! pbl layer variables
       presmix    ,& ! USL layer pressure
       dpthmx     ,& ! USL layer depth (pressure [Pa])
       zlcl       ,& ! LCL layer height[m]
       zmix       ,& ! usl layer depth [m]
       umfnewdold ) ! ratio of umf/umfold(see updraft)
    use scale_precision
    use scale_grid_index
    use scale_const,only :&
         CP => CONST_CPdry   ,  &
         PRE00 => CONST_PRE00,  &
         TEM00 => CONST_TEM00, &
         GRAV  => CONST_GRAV
    use scale_process, only: &
         PRC_MPIstop
    implicit none
    ! [IN]
    real(RP), intent(in)      :: dz_kf(KA),z_kf(KA) ! delta z and height [m]
    real(RP), intent(in)      :: qv(KA)             ! water vapor
    real(RP), intent(in)      :: qes(KA)            ! saturation water vapor
    real(RP), intent(in)      :: pres(KA)           ! pressure [Pa]
    real(RP), intent(in)      :: deltap(KA)         ! delta pressure [Pa]
    real(RP), intent(in)      :: temp(KA)           ! temperature
    real(RP), intent(in)      :: w0avg(KA)          ! running mean w
    ! [OUT]
    ! mass flux
    real(RP), intent(out)     :: umf(KA)            ! upward mass flux
    real(RP), intent(out)     :: umflcl             ! upward mass flux @lcl
    real(RP), intent(out)     :: upent(KA)          ! upward mass flux entrainment
    real(RP), intent(out)     :: updet(KA)          ! upward mass flux detrainment
    ! updraft variable
    real(RP), intent(out)     :: temp_u(KA)         ! updraft temperature
    real(RP), intent(out)     :: tempv(KA)          ! vertual temperature 
    !    real(RP)             :: tempvq_u(KA)
    real(RP), intent(out)     :: qv_u(KA)           ! updraft qv
    real(RP), intent(out)     :: totalprcp          ! total prcpitation (rain+snow)
    real(RP), intent(out)     :: cape               ! CAPE
    real(RP), intent(out)     :: cloudtop           ! cloud top height
    real(RP)                  :: cloudhight         ! cloud depth (cloud top - cloud base)
    ! water
    real(RP), intent(out)     :: qc(KA)             ! cloud water mixing ratio
    real(RP), intent(out)     :: qi(KA)             ! cloud ice   mixing ratio
    real(RP)                  :: qrout(KA)          ! rain
    real(RP)                  :: qsout(KA)          ! snow
    real(RP), intent(out)     :: qvdet(KA)          ! detrainment water vapor
    real(RP), intent(out)     :: qcdet(KA)          ! detrainment cloud water
    real(RP), intent(out)     :: qidet(KA)          ! detrainment cloud ice
    real(RP), intent(out)     :: flux_qr(KA)        !
    real(RP), intent(out)     :: flux_qs(KA)        !
    ! upward theta
    real(RP), intent(out)     :: theta_eu(KA)
    real(RP), intent(out)     :: theta_ee(KA)
    ! convection type
    integer,  intent(out)     :: I_convflag
    !> I_convflag = 0  ==> deep convection
    !>              = 1 ==> shallow convection
    !>              = 2  ==> NONE !!
    ! index valuables
    integer,  intent(out)     :: k_lcl             ! index of LCL layer
    integer,  intent(out)     :: k_top             ! index of cloud top hight
    integer,  intent(out)     :: k_ml               ! index of melt layer (temp < tem00)
    real(RP), intent(out)     :: zlcl               ! hight of lcl
    integer,  intent(out)     :: k_lc,k_let,k_pbl  ! indexs
    real(RP), intent(out)      :: zmix              ! usl layer depth [m]
    real(RP), intent(out)      :: presmix           ! usl layer depth [Pa]
    real(RP), intent(out)      :: umfnewdold(KA)    !! umfnew/umfold
    real(RP)                   :: fbfrc             !< precipitation to be feed back ( 0.0-> deep, or 1.0 ->shallow )
    real(RP), intent(out)      :: dpthmx            !< max depth of pressure
    !---------------------------------------------------------------------------

    integer, parameter :: itr_max = 10000

    !> [Internal work]
    real(RP) :: pres300       ! pressure sfc-300hpa
    !> hight variable
    integer  :: kk,nk         ! kk is "do loop index" and nk is counter of  "do loop  "
    integer  :: k_llfc        ! k of LFC height haight of lowest
    !> usl(upward source layer) check variables
    integer  :: n_uslcheck    ! usl chek layer number
    real(RP) :: pres15        ! temporaly valuables pressure 15 hpa interval
    ! calculate mix tempreature (usl has 50mb or so)
    ! usl layer variable VV
    real(RP) :: theta_mix     ! theta in usl layer
    real(RP) :: temp_mix      ! temperature in usl layer
    real(RP) :: qv_mix        ! water vaper mixing ratio in usl layer
    real(RP) :: emix          ! saturate water vaper mixing ratio in usl layer
    !
    real(RP) :: temp_dew      ! dew point temperature
    real(RP) :: templog       ! log emix/aliq
    real(RP) :: deltaz        ! dz used for interpolation of lcl
    real(RP) :: temp_env      ! environment temperature
    real(RP) :: tempv_env     ! environment temperature vertual
    real(RP) :: qvenv         ! environment q vapor
    real(RP) :: w_cz          ! velocity Kain(2004) c(z)
    real(RP) :: w_g           ! velocity Kain(2004) wg
    real(RP) :: w_lcl         ! vertical velocity @lcl
    real(RP) :: temp_lcl      ! temperature @lcl
    real(RP) :: tempv_lcl     ! vertal temperature @lcl
    real(RP) :: pres_lcl      ! pressure @ lcl
    real(RP) :: dtvv          ! delta Tvv Kain(2004)
    real(RP) :: dtrh          ! trigger rh 
    real(RP) :: radius        ! cloud radius 
    real(RP) :: dptt          ! pressure thickness 
    real(RP) :: dumfdp        ! temp var
    real(RP) :: d_min         !< minimum cloud  hight
    ! calc in subroutin kf_updraft
    real(RP) :: umfnew,umfold ! from updraft
    integer  :: k_lclm1       ! k_lcl -1
    integer  :: k_start       ! tempraly val
    integer  :: k_check(KA)   ! check layer index (because of 15mb interbal)
    integer  :: n_check       ! num of check
    integer  :: n_layers      ! num of USL layer
    integer  :: nchm          ! used shallow convection layer index
    integer  :: nuchm         ! (used shallow convection check) 
    real(RP) :: CHMAX         ! max cloud height used in shallow convection
    integer  :: NNN           ! internal work
    real(RP) :: CLDHGT(KA)    ! used for shallow convection
    real(RP) :: dpthmin       !< minimum depth of calc usl layer ??? check below
    ! trigger 3 variables
    real(RP) :: rh_lcl
    real(RP) :: U00
    real(RP) :: DQSSDT
    real(RP) :: qs_lcl

    integer :: itr
    !!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !! start cood
    !!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    pres300 = pres(KS) - DEPTH_USL*100._RP ! pressure @ surface - 300 mb. maybe 700mb or so default depth_usl is 300hPa
    tempv(:)   = temp(:)*(1._RP + 0.608_RP*qv(:)) ! vertual temperature
    ! search above 300 hPa index  to "k_llfc"
    do kk = KS, KE
       if (pres(kk) >= pres300)  k_llfc = kk
    end do
    ! usl(updraft sourcer layer) has interval for 15mb
    n_check     = KS ! first layer
    k_check(KS) = KS ! first layer
    pres15      = pres(KS) - 15.e2_RP !< pressure above 15mb
!    k_check     = KS
    ! calc 15 hpa interval Num of layer(n_check) and index of layer(k_check)
    do kk = KS+1, k_llfc
       if(pres(kk) < pres15 ) then
          n_check = n_check + 1
          k_check(n_check) = kk
          pres15 = pres15 - 15.e2_RP
       end if
    end do
    ! main routine
    ! set initial
    fbfrc      = 0._RP    ! set inital 0 
    dpthmin    = 50.e2_RP ! set initialy USL layer depth (50 hPa)
    I_convflag = 2        ! no convection set initialy
    NUCHM      = KS-1     ! initialize (used shallow convection check) !like "0"
    n_uslcheck = KS-1     ! initialize index of usl check like "0"
    nchm       = KS-1     ! initializelike "0"
    CLDHGT(:)  = 0._RP    ! cloud hight is all 0
    do itr = 1, itr_max ! usl
       n_uslcheck = n_uslcheck + 1 ! nu in WRF
       ! calc shallow convection after all layer checked
       ! NOTE:
       ! This 'if block' is only used
       if (n_uslcheck > n_check) then ! if over potentially usl layer
          if (I_convflag == 1) then   ! if schallow convection then
             chmax = 0._RP ! initialize Cloud Height Max
             NCHM  = KS-1     ! initialize index of Cloud Height Max "like 0"
             do kk = KS, n_check
                NNN = k_check(kk) ! index of checked layer (15 hpa interval layer)
                if (CLDHGT(NNN) > chmax) then ! get max cloud height in shallow convection
                   NCHM  = NNN
                   NUCHM = kk
                   CHMAX = CLDHGT(NNN) 
                end if
             end do
             n_uslcheck = nuchm - 1
             fbfrc      = 1._RP ! shallow convection is no precipitation
             cycle ! usl ! recalc updraft for shallow convection
          else ! not shallow convection then
             I_convflag = 2 ! no convection
             return
          end if ! end if schallow convection then
       end if ! end if over potentially usl layer
       k_lc     = k_check(n_uslcheck) !< condithin of Updraft Source Layer 
       n_layers = 0                   ! number of USL layer contain discretization layers
       dpthmx   = 0._RP               ! depth of USL (pressure [Pa])
       !       nk       = k_lc - KS  !< check k_lc -KS is not surface or bottom of ground
       nk       = k_lc - 1            !< check k_lc -KS is not surface or bottom of ground
       ! calculate above k_lc layers depth of pressure
       ! usl layer depth is nessesary 50mb or more
       if ( nk + 1 < KS ) then !< check k_lc -1 is not surface or bottom of ground
          if(KF_LOG) then
             if(IO_L) write(IO_FID_LOG,*) 'would go off bottom: cu_kf',pres(KS),k_lc,nk+1,n_uslcheck !,nk i,j
          end if
       else
          do ! serach USL layer index. USL layer is nessesally more 50hPa
             nk = nk +1
             if ( nk > KE ) then !< check k_lc is not index of top layer
                if(KF_LOG) then
                   if(IO_L) write(IO_FID_LOG,*) 'would go off top: cu_kf'!,nk i,j
                end if
                exit
             end if
             dpthmx   = dpthmx + deltap(nk) ! depth of pressure
             n_layers = n_layers + 1
             if (dpthmx  > dpthmin) then 
                exit
             end if
          end do
       end if
       if (dpthmx  < dpthmin) then ! if dpthmx(USL layer depth) < dptmin(minimum USL layerdepth)
          I_convflag = 2
          return ! chenge at 3d version
       end if
       k_pbl = k_lc + n_layers -1 !< pbl is index  of top of "possible" USL (until determined USL )
       ! initialize 
       theta_mix = 0._RP; temp_mix = 0._RP; qv_mix = 0._RP; zmix = 0._RP;presmix = 0._RP
       ! clc mix variable (USL layer average valuable)
       do kk = k_lc, k_pbl
          temp_mix = temp_mix + deltap(kk)*temp(kk)
          qv_mix    = qv_mix + deltap(kk)*qv(kk)
          zmix     = zmix + deltap(kk)*z_kf(kk)
          presmix  = presmix + deltap(kk)*pres(kk)
       end do
       ! mix tmp calculate
       temp_mix = temp_mix/dpthmx
       qv_mix   = qv_mix/dpthmx
       presmix  = presmix/dpthmx
       zmix     = zmix/dpthmx 
       emix     = qv_mix*presmix/(0.622_RP + qv_mix) ! water vapor pressure
       ! calculate dewpoint and LCL temperature not use look up table
       ! LCL is Lifted condensation level
       ! this dewpoint formulation is Bolton's formuration (see Emanuel 1994 4.6.2)
       templog  = log(emix/ALIQ)
       temp_dew  = (CLIQ - DLIQ*templog)/(BLIQ - templog) ! dew point temperature
       ! temperature @ lcl setting need dewpoit
       ! this algolizm proposed by Davies-Jones (1983)       
       temp_lcl = temp_dew - (0.212_RP + 1.571E-3_RP*(temp_dew - TEM00) &
            - 4.36E-4_RP*(temp_mix - TEM00))*(temp_mix - temp_dew) ! LCL temperature
       temp_lcl  = min(temp_lcl,temp_mix)
       tempv_lcl = temp_lcl*(1._RP + 0.608_RP*qv_mix) ! LCL vertual temperature
       zlcl      = zmix + (temp_lcl - temp_mix)/GdCP ! height of LCL
       ! index of z@lcl
       ! z(k_lclm1) < zlcl < z(k_lcl)
       do kk = k_lc, KE 
          k_lcl = kk
          if( zlcl <= z_kf(kk) )  exit 
       end do
       k_lclm1   = k_lcl - 1
       if (zlcl > z_kf(KE)) then   !< if not found lcl layer
          I_convflag = 2 ! no convection
          return 
       end if
       !! interpolate environment temperature and vapor at LCL
       deltaz    = ( zlcl - z_kf(k_lclm1) )/( z_kf(k_lcl)- z_kf(k_lclm1 ) )
       temp_env  = temp(k_lclm1) + ( temp(k_lcl) - temp(k_lclm1) )*deltaz
       qvenv     = qv(k_lclm1) + ( qv(k_lcl) - qv(k_lclm1) )*deltaz
       tempv_env = temp_env*( 1._RP + 0.608_RP*qvenv )
       ! w_cz set review Kain(2004) eq.2
       ! w_g set 
       ! dtlcl setting
       ! wg is 
       if (zlcl < 2.e3_RP) then ! Kain 2004 eq.2
!          w_cz = 0.02_RP*zlcl/2.e3_RP!
          w_cz = 1.e-5_RP*zlcl !
       else
          w_cz = 0.02_RP !< units m/s
       end if
       !! wg is iapproximate running mean grid resolved vertical velocity at the lcl(m/s)
       !! need w0avg
       !!wg = wg-c(z)
       w_g = (w0avg(k_lclm1) + (w0avg(k_lcl) - w0avg(k_lclm1))*deltaz )*deltax/25.e3_RP - w_cz ! need w0avg
       if ( w_g < 1.e-4_RP) then ! too small wg
          dtvv = 0._RP
       else
          dtvv = 4.64_RP*w_g**0.33_RP ! Kain (2004) eq.1 **1/3
       end if
       !! triggers will be make
       dtrh = 0._RP
       !! in WRF has 3 type of trigger function
       if ( TRIGGER == 2) then ! new function none
       else if ( TRIGGER == 3) then ! NO2007 trigger function
          !...for ETA model, give parcel an extra temperature perturbation based
          !...the threshold RH for condensation (U00)...
          ! as described in Narita and Ohmori (2007), 12th Mesoscale Conf.
          !
          !...for now, just assume U00=0.75...
          !...!!!!!! for MM5, SET DTRH = 0. !!!!!!!!
          U00 = 0.75_RP
          if(U00 < 1._RP) then
            qs_lcl=QES(k_lclm1)+(QES(k_lcl)-QES(k_lclm1))*deltaz
            rh_lcl = qvenv/qs_lcl
            DQSSDT = qv_mix*(CLIQ-BLIQ*DLIQ)/((temp_lcl-DLIQ)*(temp_lcl-DLIQ))
            if(rh_lcl >= 0.75_RP .and. rh_lcl <=0.95_RP)then
              dtrh = 0.25*(rh_lcl-0.75)*qv_mix/DQSSDT
            elseif(rh_lcl > 0.95_RP) then
               dtrh = (1._RP/rh_lcl-1._RP)*qv_mix/DQSSDT
            else
               dtrh = 0._RP
            end if
         end if
       end if
!@================================================================================          ! check
       if (temp_lcl + dtvv + dtrh < temp_env) then ! kf triggerfucn dtrh is used @ NHM trigger func
          ! parcel is not bouyant
          ! cycle and check one more up layer(15 hPa )
          cycle
       else ! parcel is bouyant ,determin updraft
          ! theta_lclm1 is need

          ! calc updraft theta_E
          call envirtht(presmix,temp_mix,qv_mix,theta_eu(k_lclm1))
          !
          if (dtvv + dtrh > 1.e-4_RP ) then
             w_lcl = 1._RP + 0.5_RP*sqrt(2._RP*GRAV*(dtvv + dtrh)*500._RP/tempv_env)! Kain(2004) eq. 3??
             w_lcl = min(w_lcl,3._RP)
          else
             w_lcl = 1._RP
          end if
          k_let = k_lcl ! add k_let
          ! interpolate pressure
          pres_lcl = pres(k_lclm1) + (pres(k_lcl) - pres(k_lclm1))*deltaz
          !          denslcl = pres_lcl/(R*tempv_lcl)
          ! temp_lcl is already caliculated
          if (w_g < 0 ) then !< Kain(2004) eq.6
             radius = 1000._RP
          elseif (w_g > 0.1_RP) then
             radius = 2000._RP
          else
             radius = 1000._RP  + 1000._RP*w_g*10._RP ! wg/0.1 -> w_g*10
          end if
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@ Compute updraft properties
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
          call CP_kf_updraft ( &
               ! [OUT]
               umf(:)           ,& ! upward mass flux
               umflcl           ,& ! updraft mass flux@LCL
               upent(:)         ,& ! updraft entrainment
               updet(:)         ,& ! updraft detrainment
               umfnewdold(:)    ,& ! ratio of umfnew/umfold
               umfnew,umfold    ,& !
               temp_u(:)        ,& ! updraft temperature internal work??
               theta_eu(:)      ,& ! updraft theta_E
               theta_ee(:)      ,& ! environment theta_E inout
               cloudhight       ,& ! cloud depth
               totalprcp        ,& ! total amount of precipitation and snow
               cloudtop         ,& ! cloud top hight
               qv_u(:)          ,& ! updraft qv
               qc(:)            ,& ! cloud water 
               qi(:)            ,& ! cloud ice
               qrout(:)         ,& ! rain
               qsout(:)         ,& ! snow
               qvdet(:)         ,& ! detrainment of qv
               qcdet(:)         ,& ! detrainment of qc
               qidet(:)         ,& ! detrainment of qc
               cape             ,& ! cape !
               flux_qr(:)       ,& ! flux of qr
               flux_qs(:)       ,& ! flux of qi
               k_top            ,& ! index of cloud top 
               k_let            ,& ! index of below detrain layer
               ! [IN]
               dz_kf(:),z_kf(:) ,&
               w_lcl            ,& ! vertical velocity of LCL layer
               temp_lcl         ,& ! temperature of LCL layer
               tempv_lcl        ,& ! vertural temperature of LCL layer
               pres_lcl         ,& ! pressure of LCL layer
               qv_mix           ,& ! vapor of USL layer
               qv(:)            ,& ! vapor
               temp(:)          ,& ! temperature
               tempv_env        ,& ! lcl temprature
               zlcl             ,& ! z[m]@LCL will be intent in
               pres(:)          ,& ! preassure
               deltap(:)        ,& ! dp
               radius           ,& ! cloud radius
               dpthmx           ,& ! pressure of depth of USL layer
               k_lcl            ,& ! index of LCL layer
               k_pbl       & ! index of USL layer
               )
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@ Compute updraft properties
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

          !! temporaly cloud hight
          !! for shallow convection
          CLDHGT(k_lc) = cloudhight
          !! minimum cloud hight for deep convection
          !! Kain (2004) eq.7
          if(temp_lcl > 293._RP) then
             D_min = 4.e3_RP
          elseif (temp_lcl < 293._RP .and. temp_lcl > 273._RP) then
             d_min = 2.e3_RP + 100._RP*(temp_lcl - tem00)
          else
             d_min = 2.e3_RP
          end if

          !! convection type check
          if ( k_top   <= k_lcl  .or. &  
               k_top   <= k_pbl  .or. &
               k_let+1 <= k_pbl ) then ! no convection allowd
             cloudhight = 0._RP
             CLDHGT(k_lc) = cloudhight ! 0
             do kk = k_lclm1,k_top
                umf(kk)     = 0._RP
                upent(kk)   = 0._RP
                updet(kk)   = 0._RP
                qcdet(kk)   = 0._RP
                qidet(kk)   = 0._RP
                flux_qr(kk) = 0._RP
                flux_qs(kk) = 0._RP
             end do
          elseif (cloudhight > d_min .and. cape > 1._RP) then ! shallow convection
             I_convflag = 0 ! deep convection
             exit ! exit usl loop
          else ! shallow convection
             !! shallow convection after determin
             !! remember this layer (by virture of non-zero CLDHGT) as potential shallo-cloudlayer
             !! detern shallow convection is after checked all potentialy USL layer
             !! when no deep convection is found (this is deterned first if block in do loop)
             I_convflag = 1
             if(n_uslcheck == NUCHM) then ! 
                exit ! exit usl loop
             else
                do kk = k_lclm1,k_top
                   umf(kk)     = 0._RP
                   upent(kk)   = 0._RP
                   updet(kk)   = 0._RP
                   qcdet(kk)   = 0._RP
                   qidet(kk)   = 0._RP
                   flux_qr(kk) = 0._RP
                   flux_qs(kk) = 0._RP
                end do
             end if ! if(n_uslcheck == NUCHM) then 
          end if ! convection type 
       end if ! triggeer
    end do ! usl
    if ( itr .ge. itr_max ) then
       write(*,*) 'xxx iteration max count was reached in the USL loop in the KF scheme'
       call PRC_MPIstop
    end if


    if (I_convflag == 1) then ! shallow convection
       k_start = max(k_pbl,k_lcl)
       k_let = k_start
    end if
    !
    !...IF THE LET AND LTOP ARE THE SAME, DETRAIN ALL OF THE UPDRAFT MASS FL
    !   THIS LEVEL...
    !
    if (k_let == k_top) then
       updet(k_top) = umf(k_top) + updet(k_top) - upent(k_top)
       qcdet(k_top) = qc(k_top)*updet(k_top)*umfnew/umfold
       qidet(k_top) = qi(k_top)*updet(k_top)*umfnew/umfold
       upent(k_top) = 0._RP
       umf(k_top)   = 0._RP
    else
       !! begin total detrainment at the level above the k_let
       dptt = 0._RP
       do kk = k_let+1,k_top
          dptt = dptt + deltap(kk)
       end do
       dumfdp = umf(k_let)/dptt
       !!
       !!... adjust mass flux profiles, detrainment rates, and precipitation fall
       !!... rates to reflect the linear decrease in mass flux between the k_let and
       !!
       !!...entrainment is allowed at every level except for K_TOP, so disallow
       !!...entrainment at k_TOP and adjust entrainment rates between k_LET and k_TOP
       !!...so the the dilution factor due to entrainment is not changed but 
       !!...the actual entrainment rate will change due due forced total 
       !!...detrainment in this layer...
       do kk = k_let+1,k_top
          if(kk == k_top) then
             updet(kk) = umf(kk-1)
             upent(kk) = 0._RP
             qcdet(kk) = updet(kk)*qc(kk)*umfnewdold(kk)
             qidet(kk) = updet(kk)*qi(kk)*umfnewdold(kk)
          else
             umf(kk)   = umf(kk-1) - deltap(kk)*dumfdp
             upent(kk) = umf(kk)*(1._RP - 1._RP/umfnewdold(kk))
             updet(kk) = umf(kk-1) - umf(kk) + upent(kk)
             qcdet(kk) = updet(kk)*qc(kk)*umfnewdold(kk)
             qidet(kk) = updet(kk)*qi(kk)*umfnewdold(kk)
          end if
          if (kk >= k_let+2) then
             totalprcp   = totalprcp - flux_qr(kk) - flux_qs(kk)
             flux_qr(kk) = umf(kk-1)*qrout(kk)
             flux_qs(kk) = umf(kk-1)*qsout(kk)
             totalprcp   = totalprcp + flux_qr(kk) + flux_qs(kk)
          end if
          !
       end do

    end if ! let== ktop

    !! initialize some arrays below cloud base and above cloud top
    !! below cloud base
    !! melt layer setting
    do kk = KS,k_top
       if(temp(kk) > TEM00) k_ml = kk !! melt layer
    end do
    !
    do kk = KS,k_lclm1
       !!
       if(kk >= k_lc) then
          !
          if (kk == k_lc) then ! if bototm of USL(pbl) 
             upent(kk) = umflcl*deltap(kk)/dpthmx 
             umf(kk)   = umflcl*deltap(kk)/dpthmx
          elseif (kk <= k_pbl) then ! if in USL(pbl)
             upent(kk) = umflcl*deltap(kk)/dpthmx
             umf(kk)   = umf(kk-1) + upent(kk) !! assume no detrain
          else
             upent(kk) = 0._RP
             umf(kk)   = umflcl
          end if
          !
          temp_u(kk) = temp_mix + (z_kf(kk) - zmix)*GdCP ! layner assumption of temp_u
          qv_u(kk)   = qv_mix 
          !             wu(kk) = w_lcl no use @wrf
       else ! below USL layer no updraft
          temp_u(kk) = 0._RP
          qv_u(kk)   = 0._RP
          umf(kk)    = 0._RP
          upent(kk)  = 0._RP
       end if
       !
       updet(kk)   = 0._RP
       qvdet(kk)   = 0._RP
       qc(kk)      = 0._RP
       qi(kk)      = 0._RP
       qrout(kk)   = 0._RP
       qsout(kk)   = 0._RP
       flux_qr(kk) = 0._RP
       flux_qs(kk) = 0._RP
       qcdet(kk)   = 0._RP
       qidet(kk)   = 0._RP
       !!calc theta_E environment
       call envirtht(pres(KK),temp(kk),qv(kk),theta_ee(kk))
       !!
    end do

    ! define variables above cloud top
    do kk = k_top+1,KE
       umf(kk)     = 0._RP
       upent(kk)   = 0._RP
       updet(kk)   = 0._RP
       qvdet(kk)   = 0._RP
       qc(kk)      = 0._RP
       qi(kk)      = 0._RP
       qrout(kk)   = 0._RP
       qsout(kk)   = 0._RP
       flux_qr(kk) = 0._RP
       flux_qs(kk) = 0._RP
       qcdet(kk)   = 0._RP
       qidet(kk)   = 0._RP
    end do
    do kk = k_top+2,KE
       temp_u(kk) = 0._RP
       qv_u(kk)   = 0._RP
    end do

    return
  end subroutine CP_kf_trigger
  !!-------------------------------------------------------------------------------
  !! updraft
  !!-------------------------------------------------------------------------------
  subroutine CP_kf_updraft (& ! 1d model updraft
       ![OUT]
       umf           ,& ! upward mass flux(updraft)
       umflcl        ,& ! upward mass flux@LCL
       upent         ,& ! updraft entrainment
       updet         ,& ! updraft detrainment
       umfnewdold    ,& ! UMF side entrainment/detrainment calculation before(old) or after(new)
       umfnew,umfold ,& ! updraft
       temp_u        ,& ! updraft temperature internal work??
       theta_eu      ,& ! updraft theta_E
       theta_ee      ,& ! thera_e @env
       cloudhight    ,& ! cloud depth [m]
       totalprcp     ,& ! total amount of precipitation and snow
       cloudtop      ,& ! cloud top hight [m]
       qv_u          ,& ! updraft qv
       qc            ,& ! cloud water 
       qi            ,& ! cloud ice
       qrout         ,& ! rain
       qsout         ,& ! snow
       qvdet         ,& ! detrainment of qv
       qcdet         ,& ! detrainment of qc
       qidet         ,& ! detrainment of qc
       cape          ,& ! cape
       flux_qr       ,& ! flux of qr
       flux_qs       ,& ! flux of qi
       k_top         ,& ! index of cloud top layer
       k_let         ,& ! index of below detrain layer
       ! [IN]
       dz_kf,z_kf    ,& ! delataz and z height
       w_lcl         ,& ! vertical velocity of LCL layer
       temp_lcl      ,& ! temperature of LCL layer
       tempv_lcl     ,& ! vertural temperature of LCL layer
       pres_lcl      ,& ! pressure of LCL layer
       qv_mix        ,& ! water vapor of USL layer
       qv            ,& ! water vapor
       temp          ,& ! temperature
       tempv_env     ,& ! temperature @ environmetnt
       zlcl          ,& ! z[m]@LCL will be intent in
       pres          ,& ! preassure [hPa]
       deltap        ,& ! dp [pa]
       radius        ,& ! cloud radius [m]
       dpthmx        ,& ! pressure of depth of USL layer
       k_lcl         ,& ! index of LCL layer
       k_pbl )          ! index of USL layer
    use scale_precision
    use scale_grid_index
    use scale_const,only :&
         CP => CONST_CPdry    , &
         PRE00 => CONST_PRE00 , &
         GRAV  => CONST_GRAV  , &
         R     => CONST_Rdry
    implicit none
    ![OUT]
    real(RP),intent(out)    :: umf(KA)        ! upward mass flux 
    real(RP),intent(out)    :: umflcl         ! upward mass flux @lcl
    real(RP),intent(out)    :: upent(KA)      ! updraft entrainment
    real(RP),intent(out)    :: updet(KA)      ! updraft detrainment
    real(RP),intent(out)    :: umfnewdold(KA) ! ratio of below umfnew/umfold
    real(RP),intent(out)    :: umfnew,umfold  ! UMF side entrainment/detrainment calculation before(old) or after(new)
    real(RP),intent(out)    :: temp_u(KA)     ! updraft temperature internal work
    real(RP),intent(out)    :: theta_ee(KA)   !< environment theta_E
    real(RP),intent(inout)  :: theta_eu(KA)   !< updraft theta_E
    real(RP),intent(out)    :: cloudhight     ! cloud depth 
    real(RP),intent(out)    :: totalprcp      ! precipitation
    real(RP),intent(out)    :: cloudtop       ! cloud top height
    real(RP),intent(out)    :: qv_u(KA)       ! updraft water vape mixing ratio
    real(RP),intent(out)    :: qc(KA)         ! liquit  QLIQ in WRF 
    real(RP),intent(out)    :: qi(KA)         ! QICE in WRF
    real(RP),intent(out)    :: qrout(KA)      ! QLQOUT in WRF
    real(RP),intent(out)    :: qsout(KA)      ! QICOUT in WRF
    real(RP),intent(out)    :: qvdet(KA)      ! will be intent out
    real(RP),intent(out)    :: qcdet(KA)      ! detrainment of qc
    real(RP),intent(out)    :: qidet(KA)      ! detrainment of qc
    real(RP),intent(out)    :: cape           ! CAPE
    real(RP),intent(out)    :: flux_qr(KA)    ! ratio of generation of liquit fallout(RAIN)
    real(RP),intent(out)    :: flux_qs(KA)    ! ratio of generation of ice fallout(SNOW)
    integer, intent(out)    :: k_top          ! top of convection layer index
    integer, intent(inout)  :: k_let          ! top of convection layer not detrain only layer
    ! [IN]
    real(RP),intent(in)     :: dz_kf(KA),z_kf(KA)
    ! lcl layer valiable
    real(RP),intent(in)     :: w_lcl 
    real(RP),intent(in)     :: temp_lcl
    real(RP),intent(in)     :: tempv_lcl
    real(RP),intent(in)     :: pres_lcl
    real(RP),intent(in)     :: qv_mix 
    real(RP),intent(in)     :: qv(KA)
    real(RP),intent(in)     :: temp(KA)
    real(RP),intent(in)     :: tempv_env ! environment vertual temperature
    real(RP),intent(in)     :: zlcl ! z[m]@LCL will be intent in
    real(RP),intent(in)     :: pres(KA)
    real(RP),intent(in)     :: deltap(KA)
    real(RP),intent(in)     :: radius
    real(RP),intent(in)     :: dpthmx ! pressure of depth of USL layer
    integer ,intent(in)     :: k_lcl
    integer ,intent(in)     :: k_pbl
    ![INTERNAL WORK]
    integer  :: kk,kkp1                         ! kk : do loop ,kkp1: kk+1
    integer  :: k_lclm1                         ! k_lcl -1
    real(RP) :: tempv_u(KA)                     ! updraft vertual temperature internalwork
    real(RP) :: tempvq_u(KA)                    ! temperature vertial for updraft ,qv, qc, qi
    real(RP) :: denslcl                         ! density @LCL
    real(RP) :: ee1,ud1, ee2,ud2                ! entrainment and detrainment calc valiables
    real(RP) :: f_eq(KA)                        ! prof5 variable
    real(RP) :: f_mix1,f_mix2                   ! factor of mixed fraction 
    real(RP) :: REI,DILBE                       ! REI KF(1990) Eq.1 , DILBE is tempvar for calc cape
    real(RP) :: qcnew,qinew                     ! qcnew is qc new var , qinew is qi newver
    real(RP) :: qfrz                            !< all frozn water
    real(RP) :: f_frozen1                       !< factor of frozen 0 to 1
    real(RP) :: temptmp                         ! temporaly temperature
    real(RP) :: temptmp_ice                     ! temporaly temperature for ice or liquid face calculateion
    real(RP) :: tempv(KA)                       ! virtual temperature  [K]
    real(RP) :: wtw                             ! w**2
    real(RP) :: boeff                           ! bouyancy effect
    real(RP) :: boterm                          ! bouyancy term for calc vertical vilocity
    real(RP) :: dztmp                           !< temporary dz
    real(RP) :: entterm                         ! entrainment term for calc vertical vilocity
    real(RP) :: theta_tmp                       ! tmporaly temperature
    real(RP) :: wu(KA)                          ! vertical velocity of updraft
    real(RP) :: qvtmp, qctmp, qitmp             !< temporaly qv
    real(RP) :: temp_u95, temp_u10              ! temporaly Temperature value use determin Mixed Fraction
    real(RP),parameter :: temp_frzT = 268.16_RP !< frozen temperature start frozen -5degC
    real(RP),parameter :: temp_frzB = 248.16_RP !< frozen temperature all frozen  -25degC
    !!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !! start cood
    !!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    tempv = temp*(1._RP + 0.608_RP*qv)
    ! initial updraft mass flux
    umfnewdold(:)    = 1._RP
    k_lclm1          = k_lcl - 1
    wtw              = w_lcl*w_lcl
    denslcl          = pres_lcl/(R*tempv_lcl) 
    umf(k_lclm1)     = denslcl*1.e-2_RP*deltax*deltax ! (A0)
    umflcl           = umf(k_lclm1)
    umfold           = umflcl
    umfnew           = umfold
    ! initial vas setup
    temp_u(k_lclm1)  = temp_lcl
    tempv_u(k_lclm1) = tempv_lcl
    qv_u(k_lclm1)    = qv_mix
    ! initialize
    ! only qv exist other water variables is none
    f_eq(k_lclm1)    = 0._RP
    upent            = 0._RP
    qc(k_lclm1)      = 0._RP
    qi(k_lclm1)      = 0._RP
    qrout(k_lclm1)   = 0._RP
    qsout(k_lclm1)   = 0._RP
    qcdet(k_lclm1)   = 0._RP
    qidet(k_lclm1)   = 0._RP
    flux_qr          = 0._RP
    flux_qs          = 0._RP
    totalprcp        = 0._RP !< total rain and snow
    !
    ee1              = 1._RP 
    ud1              = 0._RP
    rei              = 0._RP
    dilbe            = 0._RP
    cape             = 0._RP
    !! temptmp is used during clculation of the linearglaciation
    !! process. it is initially set to the temperature at which freezing
    !! is specified to begin. within the glaciation
    !! interval, it is set equal to the updraft temp at the
    !! previous model level...
    temptmp_ice = temp_frzT 
    !< updraft main loop
    !< lcl-1 to z coodinate top-1
    !< if w<0 is cloud top -> exit
    do kk = k_lclm1,KE-1 ! up_main original(wrf cood K is k_lclm1)
       kkp1 = kk + 1
       ! temporaly use below layer valuables
       f_frozen1      = 0._RP ! frozen rate initialize this variable 0 (water)to 1(ice)
       temp_u(kkp1)   = temp(kkp1)   ! up parcel temperature assumed env temperature
       theta_eu(kkp1) = theta_eu(kk) ! up parcel theta_E 
       qv_u(kkp1)     = qv_u(kk)     ! up parcel vapor
       qc(kkp1)       = qc(kk)       ! up parcel water
       qi(kkp1)       = qi(kk)       ! up parcel ice
       !!# calqulate satulation and cleate q and qic
       !!# and calc updraft pacel temperature
       !!# it is use determinnent of frozn or not
       !!#
       call tpmix2(pres(kkp1),theta_eu(kkp1),temp_u(kkp1),qv_u(kkp1),qc(kkp1),qi(kkp1),qcnew,qinew)
       !!# check to see if updraft temperature is avove the temperature at which
       !!# glaciation is assumed to initiate. if it is, calculate the
       !!# fraction of remainning liquid water to freeze... temp_frzT is the
       !!# temperature at which freezing begins, temp_frzB isthe temperature below which all
       !!# liquid water is frozen at each level...
       if (temp_u(kkp1) <= temp_frzT ) then ! if frozen temperature then ice is made
          if (temp_u(kkp1) > temp_frzB) then
             if (temptmp_ice > temp_frzT)  temptmp_ice = temp_frzT
             !!# liner assumption determin frozen ratio
             f_frozen1 = (temptmp_ice - temp_u(kkp1))/(temptmp_ice - temp_frzB)
          else
             f_frozen1 = 1._RP ! all water is frozen
          endif
          temptmp_ice = temp_u(kkp1)
          !!# calc how much ice is a layer ?
          qfrz     = (qc(kkp1) + qcnew)*f_frozen1  ! all ice
          qinew    = qinew + qcnew*f_frozen1       ! new create ice
          qcnew    = qcnew - qcnew*f_frozen1       ! new create liquit
          qi(kkp1) = qi(kkp1) + qc(kkp1)*f_frozen1 ! ice old + convert liquit to ice
          qc(kkp1) = qc(kkp1) - qc(kkp1)*f_frozen1 ! liquit old - convet liquit to ice
          !!# calculate effect of freezing
          !!# and determin new create frozen
          call DTFRZNEW(temp_u(kkp1), pres(kkp1),theta_eu(kkp1),qv_u(kkp1),qfrz,qi(kkp1) )
       end if
       tempv_u(kkp1) = temp_u(kkp1)*(1._RP + 0.608_RP*qv_u(kkp1)) ! updraft vertual temperature
       ! calc bouyancy term  for verticl velocity
       if (kk == k_lclm1) then !! lcl layer exist  between kk and kk+1 layer then use interporate value
          boeff  = (tempv_lcl + tempv_u(kkp1))/(tempv_env + tempv(kkp1)) - 1._RP
          boterm = 2._RP*(z_kf(kkp1) - zlcl)*GRAV*boeff/1.5_RP
          dztmp  = z_kf(kkp1) - zlcl
       else
          boeff  = (tempv_u(kk) + tempv_u(kkp1))/(tempv(kk) + tempv(kkp1)) - 1._RP
          boterm = 2._RP*(dz_kf(kk)        )*GRAV*boeff/1.5_RP
          dztmp  = dz_kf(kk)
       end if
       !!# entrainment term
       entterm = 2._RP*rei*wtw/umfold
       !! calc precp(qr) , snow(qr) and vertical virocity
       call p_precipitation(&
            qc(kkp1), qi(kkp1), WTW, dztmp, boterm, entterm, qcnew, qinew, qrout(kkp1), &
            qsout(kkp1),GRAV )
       !! precipitation_OC1973 or precipitation_Kessler
       !! in Ogura and Cho (1973), 'RATE' is nessesary but this rate is parameter in subroutine "precipitation_OC1973"
       !! RATE is tuning parameter
       !! read Kain 1990
       if (wtw < 1.e-3_RP ) then ! if no vertical velocity
          exit ! end calc updraft
       else
          wu(kkp1) = sqrt(wtw)
       end if
       !!# calc tehta_e in environment to entrain into updraft
       call envirtht(pres(kkp1),temp(kkp1),qv(kkp1),theta_ee(kkp1))
       ! rei is the rate of environment inflow
       rei = umflcl*deltap(kkp1)*0.03_RP/radius !!# Kain 1990 eq.1 ;Kain 2004 eq.5 
       !!--------------------------------------------------
       !! calc cape
       !!--------------------------------------------------
       tempvq_u(kkp1) = temp_u(kkp1)*(1._RP + 0.608_RP*qv_u(kkp1) - qc(kkp1) - qi(kkp1))
       if (kk == k_lclm1) then!! lcl layer exist  between kk and kk+1 then use interporate value
          dilbe = ((tempv_lcl + tempvq_u(kkp1))/(tempv_env + tempv(kkp1)) - 1._RP)*dztmp
       else
          dilbe = ((tempvq_u(kk) + tempvq_u(kkp1))/(tempv(kk) + tempv(kkp1)) - 1._RP)*dztmp
       end if
       if(dilbe > 0._RP) cape = cape + dilbe*GRAV
       !!--------------------------------------------------
       !! if cloud parcels are virtually colder then the entrainment , minimal
       !! entrainment 0.5*rei is imposed...
       !! read Kain 2004
       !!--------------------------------------------------
       !! calc entrainment/detrainment
       if(tempvq_u(kkp1) <= tempv(kkp1)) then ! if entrain and detrain
          ! original KF90 no entrainment allow
          ee2        = 0.5_RP ! Kain (2004) eq.4
          ud2        = 1._RP
          f_eq(kkp1) = 0._RP
       else
          k_let = kkp1
          temptmp   = tempvq_u(kkp1)
          !! determine teh critical mixed fraction of updraft and environmental air ...
          !! if few mix moisture air 
          f_mix1    = 0.95_RP
          f_mix2    = 1._RP - f_mix1
          theta_tmp = f_mix1*theta_ee(kkp1) + f_mix2*theta_eu(kkp1)
          qvtmp = f_mix1*qv(kkp1) + f_mix2*qv_u(kkp1)
          qctmp = f_mix2*qc(kkp1)
          qitmp = f_mix2*qi(kkp1)
          !! need only temptmp because calc bouyancy
          call tpmix2(pres(kkp1),theta_tmp,temptmp,qvtmp,qctmp,qitmp,qcnew,qinew)
          !! qinew and qcnew is damy valuavle(not use )
          temp_u95 = temptmp*(1._RP + 0.608_RP*qvtmp - qctmp - qitmp)
          ! TU95 in old coad
          if ( temp_u95 > tempv(kkp1)) then ! few mix but bouyant then ! if95
             ee2        = 1._RP ! rate of entrain is 1 -> all entrain
             ud2        = 0._RP
             f_eq(kkp1) = 1._RP
          else
             f_mix1    = 0.10_RP
             f_mix2    = 1._RP - f_mix1
             theta_tmp = f_mix1*theta_ee(kkp1) + f_mix2*theta_eu(kkp1)
             qvtmp     = f_mix1*qv(kkp1) + f_mix2*qv_u(kkp1)
             qctmp     = f_mix2*qc(kkp1)
             qitmp     = f_mix2*qi(kkp1)
             !! need only temptmp because calc bouyancy
             call tpmix2(pres(kkp1),theta_tmp,temptmp,qvtmp,qctmp,qitmp,qcnew,qinew)
             !! qinew and qcnew is damy valuavle(not use )
             temp_u10 = temptmp*(1._RP + 0.608_RP*qvtmp - qctmp - qitmp)
             if (abs(temp_u10 - tempvq_u(kkp1)) < 1.e-3_RP ) then !if10%
                ee2        = 1._RP ! all entrain
                ud2        = 0._RP
                f_eq(kkp1) = 1._RP
             else 
                f_eq(kkp1) = (tempv(kkp1) - tempvq_u(kkp1))*f_mix1 &
                     &    /(temp_u10 - tempvq_u(kkp1))
                f_eq(kkp1) = max(0._RP,f_eq(kkp1) )
                f_eq(kkp1) = min(1._RP,f_eq(kkp1) )
                if (f_eq(kkp1) == 1._RP) then ! f_eq
                   ee2 = 1._RP
                   ud2 = 0._RP
                elseif (f_eq(kkp1) == 0._RP) then
                   ee2 = 0._RP
                   ud2 = 1._RP
                else
                   !!# subroutine prof5 integrates over the gaussian dist to determine
                   !!# the fractional entrainment and detrainment rates...
                   call prof5(f_eq(kkp1),ee2,ud2)
                end if ! end of f_iq
             end if ! end of if10%
          end if ! end of if95%
       end if! end of entrain/detrain
       ee2 = max(ee2,0.5_RP)
       ud2 = 1.5_RP*ud2
       !
       upent(kkp1) = 0.5*rei*(ee1 + ee2)
       updet(kkp1) = 0.5*rei*(ud1 + ud2)
       !! if the calculated updraft detrainment rate is greater then the total
       !! updraft mass flux(umf), all cloud mass detrains
       if (umf(kk) - updet(kkp1) < 10._RP) then ! why 10.???
          !! if the calculated detrained mass flux is greater than totla upd mass
          !! flux, impose total detrainment of updraft mass at the previous model level..
          !! first, correct cape calculation if neede
          if (dilbe > 0._RP) then !< colect cape 
             cape = cape - dilbe*GRAV
          end if
          k_let = kk ! then k_let = k_top
          exit
       else
          ee1              = ee2
          ud1              = ud2
          umfold           = umf(kk) - updet(kkp1)
          umfnew           = umfold  + upent(kkp1)
          umf(kkp1)        = umfnew
          umfnewdold(kkp1) = umfnew/umfold
          qcdet(kkp1)      = qc(kkp1)*updet(kkp1)
          qidet(kkp1)      = qi(kkp1)*updet(kkp1)
          qvdet(kkp1)      = qv_u(kkp1)
          !! below layer updraft qv and entrain q /new updraft massflux
          qv_u(kkp1)       = (umfold*qv_u(kkp1) + upent(kkp1)*qv(kkp1))/umfnew
          theta_eu(kkp1)   = (umfold*theta_eu(kkp1) + upent(kkp1)*theta_ee(kkp1))/umfnew
          ! in WRF
          ! PPTLIQ(KKP11)  = qlqout(KKP11)*UMF(KKP1)
          qc(kkp1)         = qc(kkp1)*umfold/umfnew
          qi(kkp1)         = qi(kkp1)*umfold/umfnew
          ! flux_qr is ratio of generation of liquit fallout(RAIN)
          ! flux_qi is ratio of generation of ice fallout(SNOW)
          flux_qr(kkp1)    = qrout(kkp1)*umf(kk)
          flux_qs(kkp1)    = qsout(kkp1)*umf(kk)
          totalprcp        = totalprcp + flux_qr(kkp1) + flux_qs(kkp1) ! total prcp vars
          !! if below cloud base then
          !! updraft entrainment is umf@LCL*dp/depth
          if(kkp1 <= k_pbl ) upent(kkp1) = upent(kkp1) + umflcl*deltap(kkp1)/dpthmx
       end if
    end do        ! this loop is exit at  w=0.
    k_top = kk ! vertical coodinate index @ w=0
    ! cloud height
    cloudhight = z_kf(k_top) - zlcl
    cloudtop   = z_kf(k_top)
    return
  end subroutine CP_kf_updraft
  !!------------------------------------------------------------------------------
  !! Down draft
  !!------------------------------------------------------------------------------
  subroutine CP_kf_downdraft ( &
       ![IN]
       dz_kf,z_kf ,& ! 
       u          ,& ! x-direction wind
       v          ,& ! meridional wind
       zlcl       ,& ! lcl_hight [m]
       rh         ,& ! relative humidity initial make kf_init 
       deltap     ,& ! delta P
       pres       ,& ! pressure [Pa]
       qv         ,& ! watervapor
       ems        ,& ! 
       emsd       ,& !
       theta_ee   ,& ! environment equivalent theta
       umf        ,& ! upward mass flux
       totalprcp  ,& ! total precipitation
       flux_qs    ,& !
       tempv      ,&
       I_convflag ,& ! convection flag
       ! index
       k_lcl      ,&
       k_ml       ,&
       k_top      ,&
       k_pbl      ,&
       k_let      ,&
       k_lc       ,&
       ![OUT]
       wspd       ,& ! wind speed 1 k_lcl, 2 k_z5,3 k_top
       dmf        ,& ! down
       downent    ,&
       downdet    ,&
       theta_d    ,&
       qv_d       ,&
       prcp_flux  ,&
       k_lfs      ,&
       CPR        ,&
       tder)
    use scale_precision
    use scale_grid_index
    use scale_const,only :&
         CP => CONST_CPdry    , &
         PRE00 => CONST_PRE00 , &
         EMELT => CONST_EMELT , &
         GRAV  => CONST_GRAV  , &
         R     => CONST_Rdry
    use scale_atmos_saturation ,only :&
         ATMOS_SATURATION_psat_liq
    use scale_process, only: &
         PRC_MPIstop
    implicit none
    ! [IN]
    ! from init
    real(RP),intent(in) :: dz_kf(KA),z_kf(KA)
    real(RP),intent(in) :: u(KA)      ! x velocity 
    real(RP),intent(in) :: v(KA)      ! y velocity
    real(RP),intent(in) :: rh(KA)     ! ! initial make kf_init 
    real(RP),intent(in) :: deltap(KA) ! delta pressure
    real(RP),intent(in) :: pres(KA)   ! pressure
    real(RP),intent(in) :: qv(KA)     ! watervapor
    real(RP),intent(in) :: ems(KA)    ! dp*(deltax)**2/g
    real(RP),intent(in) :: emsd(KA)   ! 1._RP/ems
    ! from kf_triger
    real(RP),intent(in) :: zlcl         ! lcl_hight
    real(RP),intent(in) :: umf(KA)      ! upward mass flux
    real(RP),intent(in) :: totalprcp    ! in !total precipitation
    real(RP),intent(in) :: flux_qs(KA)  ! snow fall
    real(RP),intent(in) :: tempv(KA)    ! virtual env temperature  internal work
    real(RP),intent(in) :: theta_ee(KA) ! in ! environment equivalent theta
    integer  :: k_z5                    ! internalwork 500hPa layer index
    ! from kf_triger
    integer,intent(in)  :: I_convflag ! index of convection
    integer,intent(in)  :: k_lcl      ! index of lcl layer
    integer,intent(in)  :: k_top      ! index of cloud top layer
    integer,intent(in)  :: k_pbl      ! index of USL layer top
    integer,intent(in)  :: k_let      ! index of updraft only detrainment 
    integer,intent(in)  :: k_lc       !  index of USL layer bottom
    integer,intent(in)  :: k_ml       ! index of melt layer 
    ! [OUT]
    real(RP),intent(out) :: wspd(3)     ! wind speed 1 k_lcl, 2 k_z5,3 k_top
    real(RP),intent(out) :: dmf(KA)     ! downdraft massflux
    real(RP),intent(out) :: downent(KA) ! downdraft entrainment
    real(RP),intent(out) :: downdet(KA) ! downdraft detrainment
    real(RP),intent(out) :: theta_d(KA) ! downdraft theta
    real(RP),intent(out) :: qv_d(KA)    ! vapor of downdraft
    real(RP),intent(out) :: prcp_flux   ! out! precpitation flux
    real(RP),intent(out) :: tder        ! temperature change from evap
    real(RP),intent(out) :: CPR         ! all precipitation  before consider cloud bottom evaporation
    integer,intent(out) :: k_lfs        ! LFS layer index
    ! [INTERNAL WORK]
    integer :: kk, kkp1   ! do loop variable
    real(RP) :: shsign    ! frip variable tmp use
    real(RP) :: vws       ! vertical wind shear
    real(RP) :: pef       ! precipitation efficiency
    real(RP) :: pef_wind  ! precipitation efficiency from wind shear
    real(RP) :: pef_cloud ! precipitation efficiency from cloud hight
    real(RP) :: cbh       ! cloud base hight
    real(RP) :: rcbh      ! E @ Zhang and Fritsch 1986 eq.13
    !
    real(RP) :: dens_d    ! density of downdraft @ LFS
    real(RP) :: rhbar     ! rerative humidity RH, bar is mean of
    real(RP) :: rhh       ! decreace rhh
    !! humidity adjust ment tempolaly variables VV
    real(RP) :: dssdt
    real(RP) :: dtmp
    real(RP) :: qsrh
    real(RP) :: RL
    real(RP) :: T1rh
    !! ^^^^^^^^^^^^^
    real(RP) :: dpthtmp      ! temporaly depth of downdraft source layer
    real(RP) :: dpthdet      ! downdraft detrainment  depth (Pa)
    real(RP) :: temp_d(KA)   ! downdraft temperature 
    real(RP) :: tempv_d(KA)  ! downdraft virtual temperature
    real(RP) :: theta_ed(KA) ! downdraft equivalent theta
    real(RP) :: qvsd(KA)     ! saturate watervapor in downdraft
    real(RP) :: qvs_tmp      ! saturate qv temp
    real(RP) :: exn(KA)      ! exner function
    real(RP) :: f_dmf        ! factor of dmf Kain(2004) eq.11
    real(RP) :: dtempmlt     ! melting effect of temperatture
    real(RP) :: es           ! saturate vapor pressure
    real(RP) :: ddinc        ! downdraft ??atode
    real(RP) :: prcpmlt      ! melt precp internal
    integer :: k_dstart      ! 150mb below of downdraft start (near USL top)
    integer :: k_lfstmp      ! downdraft start layer
    integer :: k_ldt         ! index top of  downdrraft detraiment layer
    integer :: k_ldb         ! index botom of downdraft detrainmengt layer
    !!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !! start cood
    !!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    ! if no convection 
    if (I_convflag == 2) return ! if 3d, may be change
    ! detremin 
    do kk = KS, KE
       if (pres(kk) >= pres(KS)*0.5_RP)   k_z5 = kk
    end do
    ! calc  wind speed at LCL and cloud top and
    wspd(1) = sqrt(u(k_lcl)*u(k_lcl) + v(k_lcl)*v(k_lcl))
    wspd(2) = sqrt(u(k_z5)*u(k_z5)   + v(k_z5)*v(k_z5))
    wspd(3) = sqrt(u(k_top)*u(k_top) + v(k_top)*v(k_top))
    ! calc wind shear and precipitation efficiency
    ! precipitation efficiency is calc two way
    ! 1. wind shear 2.cloud base hight
    ! we use mean of them
    if(wspd(3) > wspd(1)) then
       shsign = 1._RP
    else
       shsign = -1._RP
    end if
    vws = (u(k_top) - u(k_lcl))*(u(k_top) - u(k_lcl)) &
         + (v(k_top) - v(k_lcl))*(v(k_top) - v(k_lcl))

    vws = 1.e3_RP*shsign*sqrt(vws)/(z_kf(k_top) - z_kf(k_lcl))
    !! PEF(precipitation efficiency ) from Fritsch and Chappel 1980
    ! wind shear type
    pef_wind = 1.591_RP + vws*(-0.639_RP + vws*(9.53e-2_RP - vws*4.96e-3_RP) )
    ! 0.2 < pef_wind< 0.9
    pef_wind = max(pef_wind,0.2_RP)
    pef_wind = min(pef_wind,0.9_RP)

    !! 2. cloud base height
    ! refarence Zhang and Fritsch 1986
    cbh = (zlcl - z_kf(KS))*3.281e-3_RP ! convert m-> feet that's Amaerican
    if( cbh < 3._RP) then
       rcbh = 0.02_RP
    else
       rcbh = 0.96729352_RP + CBH*(-0.70034167_RP + CBH*(0.162179896_RP + CBH*(-            &
            1.2569798E-2_RP + CBH*(4.2772E-4_RP - CBH*5.44E-6_RP))))
    end if
    if(cbh > 0.25) rcbh = 2.4_RP
    pef_cloud = 1._RP/(1._RP + rcbh)
    pef_cloud = min(pef_cloud,0.9_RP)
    ! pef is mean of wind shear type and  cloud base heigt type
    pef = 0.5_RP*(pef_wind + pef_cloud)
    !!--------------------------------------------------
    !! compute downdraft propeties
    !!--------------------------------------------------
    tder = 0._RP ! initialize evapolation valuavle
    if(I_convflag == 1) then ! down devap
       k_lfs = KS ! shallow convection
    else
       !! start downdraft about 150 hPa above cloud base
       k_dstart = k_pbl + 1 ! index of usl layer top  +1
       k_lfstmp = k_let - 1
       do kk = k_dstart+1, KE
          dpthtmp = pres(k_dstart) - pres(kk)
          if(dpthtmp > 150.e2_RP) then ! if dpth > 150hpa then
             k_lfstmp = kk
             exit
          end if
       end do
       k_lfstmp = min(k_lfstmp, k_let - 1) ! k_let  is equivalent temperature layer index then
       k_lfs    = k_lfstmp

       !... if LFS in not at least 50 mb above cloud base (implying that the
       !... level of equil temp, let ,is just above cloud base) do not allow downdraft
       ! dondraft layer
       ! downdraft is saturated 
       if(pres(k_dstart) - pres(k_lfs) > 50.e2_RP) then   ! LFS > 50mb(minimum of downdraft source layer)
          theta_ed(k_lfs) = theta_ee(k_lfs)
          qv_d(k_lfs)     = qv(k_lfs)
          !! call tpmix2dd to find wet-bulb temperature and qv
          call tpmix2dd(pres(k_lfs),theta_ed(k_lfs),temp_d(k_lfs),qvs_tmp)
          call calcexn(pres(k_lfs),qvs_tmp,exn(k_lfs))
          !!          exn(kk) = (PRE00/pres(k_lfs))**(0.2854*(1._RP - 0.28_RP*qv_d(k_lfs)))
          !!                theta_d(k_lfs) = temp_d(k_lfs)*(PRE00/pres(k_lfs))**(0.2854*(1._RP -
          !!                0.28*qv_d(k_lfs)))
          theta_d(k_lfs) = temp_d(k_lfs)*exn(k_lfs)
          !! take a first guess at hte initial downdraft mass flux
          tempv_d(k_lfs) = temp_d(k_lfs)*(1._RP + 0.608_RP*qvs_tmp)
          dens_d         = pres(k_lfs)/(R*tempv_d(k_lfs))
          dmf(k_lfs)     = -(1._RP - pef)*1.e-2_RP*deltax*deltax*dens_d ! AU0 = 1.e-2*dx**2
          downent(k_lfs) = dmf(k_lfs)
          downdet(k_lfs) = 0._RP
          rhbar          = rh(k_lfs)*deltap(k_lfs)
          dpthtmp        = deltap(k_lfs)
          !! calc downdraft entrainment and (detrainment=0.)
          !! and dmf
          !! downdraft theta and q
          do kk = k_lfs-1,k_dstart,-1
             kkp1         = kk + 1
             downent(kk)  = downent(k_lfs)*ems(kk)/ems(k_lfs)
             downdet(kk)  = 0._RP
             dmf(kk)      = dmf(kkp1) + downent(kk)
             theta_ed(kk) = ( theta_ed(kkp1)*dmf(kkp1) + theta_ee(kk)*downent(kk) )/dmf(kk)
             qv_d(kk)     = ( qv_d(kkp1)*dmf(kkp1)  + qv(kk)*downent(kk) )/dmf(kk)
             dpthtmp      = dpthtmp + deltap(kk)
             rhbar        = rhbar + rh(kk)*deltap(kk) ! rh average@ downdraft layer
          end do
          rhbar   = rhbar/dpthtmp
          f_dmf   = 2._RP*(1._RP - rhbar) ! Kain 2004 eq.11
          !
          !! calc melting effect
          !! first, cocalc total frozen precipitation generated..
          prcpmlt = 0._RP
          do kk = k_lcl,k_top
             prcpmlt = prcpmlt + flux_qs(kk)
          end do
          if (k_lc < k_ml ) then ! if below melt level layer then
             !!             RLF is EMELT
             !!             dtempmlt = RLF*prcpmlt/(CP*umf(k_lcl))
             dtempmlt = EMELT*prcpmlt/(CP*umf(k_lcl))
          else
             dtempmlt = 0._RP
          end if
          call tpmix2dd(pres(k_dstart),theta_ed(k_dstart),temp_d(k_dstart),qvs_tmp)
          temp_d(k_dstart) = temp_d(k_dstart) - dtempmlt
          !! use check theis subroutine is this
          call ATMOS_SATURATION_psat_liq(es,temp_d(k_dstart)) !saturation vapar pressure
          qvs_tmp = 0.622_RP*es/(pres(k_dstart) - es )
          !! Bolton 1980 pseudoequivalent potential temperature
          theta_ed(k_dstart) = temp_d(k_dstart)*(PRE00/pres(k_dstart))**(0.2854_RP*(1._RP - 0.28_RP*qvs_tmp))*   &
               &    exp((3374.6525_RP/temp_d(k_dstart)-2.5403_RP)*qvs_tmp*(1._RP + 0.81_RP*qvs_tmp))
          k_ldt   = min(k_lfs-1, k_dstart-1)
          dpthdet = 0._RP
          do kk = k_ldt,KS,-1 ! start calc downdraft detrain layer index
             !!
             dpthdet      = dpthdet + deltap(kk)
             theta_ed(kk) = theta_ed(k_dstart)
             qv_d(kk)     = qv_d(k_dstart)
             call tpmix2dd(pres(kk),theta_ed(kk),temp_d(kk),qvs_tmp)
             qvsd(kk)     = qvs_tmp
             !... specify RH decrease of 20%/km indowndraft
             rhh = 1._RP - 2.E-4_RP*(z_kf(k_dstart) -z_kf(kk) ) ! 0.2/1000.
             !!
             !!... adjust downdraft temp,qv to secified RH :
             !!
             !! Kain(2004)
             if(rhh < 1._RP) then
                !!
                dssdt = (cliq - bliq*dliq)/((temp_d(kk) - dliq)*(temp_d(kk) - dliq))
                RL    = XLV0 - XLV1*temp_d(kk)
                dtmp  = RL*qvs_tmp*(1._RP - rhh )/(CP + RL*rhh*qvs_tmp*dssdt )
                T1rh  = temp_d(kk) + dtmp
                call ATMOS_SATURATION_psat_liq(es,T1rh) !saturation vapar pressure
                es = RHH*es
                qsrh = 0.622_RP*es/(pres(kk) - es)
                if(qsrh < qv_d(kk) ) then
                   qsrh = qv_d(kk)
                   t1rh = temp_d(kk) + (qvs_tmp - qsrh)*RL/CP
                end if

                temp_d(kk) = t1rh
                qvs_tmp    = qsrh
                qvsd(kk)   = qvs_tmp
                !!
             end if
             !!
             tempv_d(kk) = temp_d(kk)*( 1._RP + 0.608_RP*qvsd(kk) )
             if(tempv_d(kk) > tempv(kk) .or. kk == KS) then
                k_ldb = kk
                exit
             end if
             !!
          end do ! end calc downdraft detrain layer index
          if((pres(k_ldb)-pres(k_lfs)) > 50.e2_RP ) then ! minimum downdraft depth !
             do kk = k_ldt,k_ldb,-1 
                kkp1        = kk + 1
                downdet(kk) = -dmf(k_dstart)*deltap(kk)/dpthdet
                downent(kk) = 0._RP
                dmf(kk)     = dmf(kkp1) + downdet(kk)
                tder        = tder + (qvsd(kk) - qv_d(kk))*downdet(kk)
                qv_d(kk)    = qvsd(kk)
                call calcexn(pres(kk),qv_d(kk),exn(kk))
                theta_d(kk) = temp_d(kk)*exn(kk)
             end do
          end if
       end if ! LFS>50mb
    end if ! down devap
    !!
    !!... if downdraft does not evaporate any water for specified relative
    !!...humidity, no downdraft is allowed ..
    !!and shallow convection does not have downdraft
    if (tder < 1._RP) then ! dmf modify
       prcp_flux = totalprcp  
       cpr       = totalprcp
       tder      = 0._RP
       k_ldb     = k_lfs
       do kk = KS, k_top 
          dmf(kk)     = 0._Rp
          downdet(kk) = 0._Rp
          downent(kk) = 0._Rp
          theta_d(kk) = 0._Rp
          temp_d(kk)  = 0._Rp
          qv_d(kk)    = 0._Rp
       end do
    else
       ddinc = -f_dmf*umf(k_lcl)/dmf(k_dstart) ! downdraft keisuu Kain 2004 eq.11
       if(tder*ddinc > totalprcp) then
          ddinc = totalprcp/tder ! rate of prcp/evap
       end if
       tder = tder*ddinc
       do kk = k_ldb,k_lfs
          dmf(kk)     = dmf(kk)*ddinc
          downent(kk) = downent(kk)*ddinc
          downdet(kk) = downdet(kk)*ddinc
       end do
       cpr       = totalprcp
       prcp_flux = totalprcp - tder ! precipitation - evapolate water
       pef       = prcp_flux/totalprcp
       !
       !...ADJUST UPDRAFT MASS FLUX, MASS DETRAINMENT RATE, AND LIQUID WATER AN
       !   DETRAINMENT RATES TO BE CONSISTENT WITH THE TRANSFER OF THE ESTIMATE
       !   FROM THE UPDRAFT TO THE DOWNDRAFT AT THE LFS...
       !     
       !         DO NK=LC,LFS
       !           UMF(NK)=UMF(NK)*UPDINC
       !           UDR(NK)=UDR(NK)*UPDINC
       !           UER(NK)=UER(NK)*UPDINC
       !           PPTLIQ(NK)=PPTLIQ(NK)*UPDINC
       !           PPTICE(NK)=PPTICE(NK)*UPDINC
       !           DETLQ(NK)=DETLQ(NK)*UPDINC
       !           DETIC(NK)=DETIC(NK)*UPDINC
       !         ENDDO
       !     
       !...ZERO OUT THE ARRAYS FOR DOWNDRAFT DATA AT LEVELS ABOVE AND BELOW THE
       !...DOWNDRAFT...
       !
       if (k_ldb > KS) then
          do kk = KS,k_ldb-1
             dmf(kk)     = 0._RP
             downdet(kk) = 0._RP
             downent(kk) = 0._RP
             theta_d(kk) = 0._RP
             temp_d(kk)  = 0._RP
             qv_d(kk)    = 0._RP
          end do
       end if
       ! no dmf is above LFS layer
       do kk = k_lfs+1,KE
          dmf(kk)     = 0._RP
          downdet(kk) = 0._RP
          downent(kk) = 0._RP
          theta_d(kk) = 0._RP
          temp_d(kk)  = 0._RP
          qv_d(kk)    = 0._RP
       end do
       !! no temperature and qv avave downdraft detrainment startlayer (k_ldt)
       do kk = k_ldt+1,k_lfs-1
          temp_d(kk)  = 0._RP
          qv_d(kk)    = 0._RP
          theta_d(kk) = 0._RP
       end do
    end if ! dmf modify

    return
  end subroutine CP_kf_downdraft
  !!------------------------------------------------------------------------------
  !! compute properties for compensational subsidence
  !! return new valuavles
  !!------------------------------------------------------------------------------
  !!compute properties for compensational subsidence
  !!and determin dt variables
  subroutine CP_kf_compensational (&
       !![in]
       dz_kf,z_kf  ,&
       pres        ,& ! pressure
       deltap      ,& ! deltap
       temp_bf     ,& ! temperature
       qv          ,& ! water vapor
       ems         ,& 
       emsd        ,&
       !! trigger
       presmix     ,& ! usl layer pressure depth
       zmix        ,& ! usl layer pressure depth
       dpthmx      ,& ! usl layer depth
       cape        ,& !cape
       temp_u      ,& ! updraft temperature
       qvdet       ,& ! updraft water vapor
       umflcl      ,& ! umf @LCL
       umf         ,& ! UMF
       upent       ,& ! updraft entrainment
       updet       ,& ! updraft detrainment
       qcdet       ,& ! updraft detrainment qc
       qidet       ,& ! updraft detrainment qi
       qc          ,& ! cloud water
       qi          ,& ! cloud ice
       flux_qr     ,& ! rain flux
       flux_qs     ,& ! snow flux
       umfnewdold  ,& ! 1/(umfnew/umfold ratio) !
       !!downdraft
       wspd        ,& ! wind speed 1. 2 is used
       dmf         ,& ! DMF
       downent     ,& ! downdraft entrainment 
       downdet     ,& ! downdraft detrainment
       qv_d        ,& ! downdraft detrainment qv
       theta_d     ,& ! potential temperature @downdraft
       prcp_flux   ,& ! precipitation
       tder        ,& ! evapolation
       cpr         ,& ! all precipitation  before consider cloud bottom evaporation
       ! 
       I_convflag  ,& ! intent inout convective flag
       !triger
       k_top       ,&
       k_lcl_bf    ,&
       k_lc        ,&
       k_pbl       ,&
       k_ml        ,&
       !downdraft
       k_lfs       ,&
       !![OUT] new valuavles after timestep
       temp_g      ,&
       qv_g        ,&
       qc_nw       ,&
       qi_nw       ,&
       qr_nw       ,&
       qs_nw       ,&
       rainrate_cp ,&
       nic         ,&
       cldfrac_KF  ,& 
       timecp      ,&
       time_advec)
    use scale_precision
    use scale_grid_index
    use scale_const,only :&
         CP => CONST_CPdry    , &
         PRE00 => CONST_PRE00 , &
         GRAV  => CONST_GRAV  , &
         EMELT => CONST_EMELT
    use scale_atmos_saturation ,only :&
         ATMOS_SATURATION_psat_liq
    use scale_time , only :&
         dt => TIME_DTSEC_ATMOS_PHY_CP
    use scale_process, only: &
         PRC_MPIstop
    implicit none
    ! in
    ! form init
    real(RP),intent(in)    :: dz_kf(KA),z_kf(KA) ! delta Z, and haight [m] full point
    real(RP),intent(in)    :: pres(KA)           ! pressure [Pa]
    real(RP),intent(in)    :: deltap(KA)         ! delta pressure
    real(RP),intent(in)    :: temp_bf(KA)        ! temperature berore
    real(RP),intent(in)    :: qv(KA)             ! cloud vaper mixing ratio
    real(RP),intent(in)    :: ems(KA)
    real(RP),intent(in)    :: emsd(KA)
    ! from trigger
    real(RP),intent(in)    :: presmix
    real(RP),intent(in)    :: zmix
    real(RP),intent(in)    :: dpthmx             ! usl layer depth
    real(RP),intent(in)    :: cape
    real(RP),intent(in)    :: temp_u(KA)
    real(RP),intent(in)    :: qvdet(KA)
    real(RP),intent(in)    :: umflcl             ! umf@LCL internal work(umf(k_lcl - 1)
    real(RP),intent(inout) :: umf(KA)            ! UMF
    real(RP),intent(inout) :: upent(KA)          ! updraft entrainment
    real(RP),intent(inout) :: updet(KA)          ! updraft detrainment
    real(RP),intent(inout) :: qcdet(KA)          ! updraft detrainment qc
    real(RP),intent(inout) :: qidet(KA)          ! updraft detrainment qi
    real(RP),intent(in)    :: qc(KA)             ! cloud water mixing ratio
    real(RP),intent(in)    :: qi(KA)             ! cloud ice mixing ratio
    real(RP),intent(in)    :: flux_qr(KA)
    real(RP),intent(in)    :: flux_qs(KA)
    real(RP),intent(in)    :: umfnewdold(KA)     ! umfnew/umfold ratio
    ! from dwndraft 
    real(RP),intent(in)    :: wspd(3)            ! wind speed 1. 2 is used
    real(RP),intent(inout) :: dmf(KA)            ! DMF
    real(RP),intent(inout) :: downent(KA)        ! downdraft entrainment 
    real(RP),intent(inout) :: downdet(KA)        ! downdraft detrainment
    real(RP),intent(in)    :: qv_d(KA)           ! downdraft detrainment qv
    real(RP),intent(in)    :: theta_d(KA)
    real(RP),intent(inout) :: prcp_flux
    real(RP),intent(inout) :: tder               ! evapolation
    real(RP),intent(in)    :: cpr
    integer,intent(inout)  :: I_convflag         ! intent inout
    ! from trigger
    integer,intent(in)     :: k_top
    integer,intent(inout)  :: k_lcl_bf
    integer,intent(in)     :: k_lc
    integer,intent(in)     :: k_pbl
    integer,intent(in)     :: k_ml
    integer,intent(in)     :: k_lfs
    ! out
    real(RP),intent(out)   :: temp_g(KA)         ! temporarly temperature -> after iteration then new variable 
    real(RP),intent(out)   :: qv_g(KA)
    real(RP),intent(out)   :: qc_nw(KA)
    real(RP),intent(out)   :: qi_nw(KA)
    real(RP),intent(out)   :: qr_nw(KA)
    real(RP),intent(out)   :: qs_nw(KA)
    integer,intent (out)   :: nic                ! num of step convection active [step]
    ! dt(cumulus parametrization deltat)/cumulus parameterization time scale
    real(RP),intent(out)   :: rainrate_cp        !PPTFLX(prcp_flux)/deltax**2/dt convective rain rate
    real(RP),intent(out)   :: cldfrac_KF(KA,2)   ! 1.shallow and 2.deep
    real(RP),intent(out)   :: timecp             ! timescale of cumulus parameterization
    real(RP),intent(out)   :: time_advec         ! advection timescale
    ![INTERNAL WORK]
    ! iteratevar
    real(RP) :: umf2(KA)     ! UMF
    real(RP) :: upent2(KA)   ! updraft entrainment
    real(RP) :: updet2(KA)   ! updraft detrainment
    real(RP) :: qcdet2(KA)   ! updraft detrainment qc
    real(RP) :: qidet2(KA)   ! updraft detrainment qi
    real(RP) :: dmf2(KA)     ! DMF
    real(RP) :: downent2(KA) ! downdraft entrainment 
    real(RP) :: downdet2(KA) ! downdraft detrainment
    real(RP) :: prcp_flux2   ! precpitation flux
    real(RP) :: tder2        ! evaporation
    !
    real(RP) :: tkemax       !tkemax tuning prameter
    real(RP) :: z_lcl        ! lcl layer hight
    real(RP) :: theta(KA)    ! theta is not same as SCALE theta. This theta is only assume qv
    real(RP) :: theta_u(KA)  ! theta in updraft
    real(RP) :: theta_eu(KA) ! equivalent PT updraft
    real(RP) :: theta_eg(KA) ! equivalent PT environment
    real(RP) :: exn(KA)      ! exner function
    real(RP) :: qv_env       ! environment qv
    real(RP) :: qv_mix       ! USL layer mean  qv
    real(RP) :: qv_gu(KA)    ! updraft qv (used calc CAPE)
    real(RP) :: temp_env     ! temporarly environment temperature lcl layer
    real(RP) :: tempv_env    ! temporarly environment virtual temperature lcl layer
    real(RP) :: temp_lcl     ! LCL temperature used calcurate CAPE
    real(RP) :: tempv_lcl    ! temporarly environment virtual temperature lcl layer
    real(RP) :: temp_mix     ! temporarly environment temperature USL layer mean
    real(RP) :: temp_gu(KA)  ! temporarly updraft temperature
    real(RP) :: tempv_g(KA)  ! temporaly virtual
    real(RP) :: tempvq_u(KA) ! temporaly virtual
    real(RP) :: es           ! saturate vapor pressure
    real(RP) :: qvss         ! saturate vapor pressure mixingratio
    ! calc dew point temperature etc.
    real(RP) :: DQ,TDPT,DSSDT,emix,RL,TLOG
    !
    real(RP) :: vconv
    real(RP) :: dzz           ! lcl layer depth  used calcurate CAPE for interpolate lcl layer
    real(RP) :: deltaz        ! lcl layer depth  used calcurate CAPE for interpolate lcl layer
    real(RP) :: dilbe         ! used calculate CAPE
    integer  :: k_lcl,k_lclm1 ! LCL and LCL-1 layer index used calucrat CAPE
    ! inter porlatevariable
    ! new variables
    real(RP) :: theta_nw(KA)      ! new PT (used iteration
    real(RP) :: theta_g(KA)       ! new PT
    real(RP) :: qv_nw(KA)         ! tempolaly new qv, becouse itaration
    real(RP) :: dpth              ! pressure depth of cloud used check 
    real(RP) :: cape_g            ! new cape clculate after timestep 
    real(RP) :: dcape             ! deltacape compared 10% of original
    real(RP) :: fxm(KA)           ! flux factor
    real(RP) :: f_cape ,f_capeold ! cape ratio new/old
    real(RP) :: stab              ! 0.95 stablevariable
    real(RP) :: dtt_tmp
    real(RP) :: dtt               ! deltat of cumulus prameterization determined by layer depth/omega
    real(RP) :: deltat            ! deltat of cumulus parameterization same as dtt but reduced
    integer :: istop              !  effor flag 1 mass blance check
    integer :: ncount             ! countor for iteration
    integer :: ntimecount         ! timecount do loop index
    integer :: nstep              ! max step of ntimecount 
    integer :: noiter             ! iteration flag
    ! temporaly
    integer :: kk,kkp1                       ! index of do loop
    integer :: k_lmax                        ! max of k_lcl,k_lfs
    real(RP) :: tma,tmb,tmm                  ! eff
    real(RP) :: acoeff,bcoeff                ! eff
    real(RP) :: evac                         ! shallow convection TKE factor
    real(RP) :: ainc,ainctmp, aincmx,aincold ! factors ainctmp is tmpvariable; aincmx is max of ainc
    real(RP) :: omg(KA)                      ! pressure velocity
    real(RP) :: topomg                       ! cloud top omg calc by updraft
    real(RP) :: fbfrc                        ! precpitation  to be fedback 0.0 -1.0 shallo-> 1.0(no rain) deep->0.0
    real(RP) :: frc2
    real(RP) :: dfda
    real(RP) :: domg_dp(KA)                  ! d omega/dp
    real(RP) :: absomgtc,absomg
    real(RP) :: f_dp
    ! fluxs
    real(RP) :: theta_fxin(KA)  ! cmpensational subsidence flux form
    real(RP) :: theta_fxout(KA) ! cmpensational subsidence flux form
    real(RP) :: qv_fxin(KA)     ! cmpensational subsidence flux form
    real(RP) :: qv_fxout(KA)    ! cmpensational subsidence flux form
    real(RP) :: qc_fxin(KA)     ! cmpensational subsidence flux form
    real(RP) :: qc_fxout(KA)    ! cmpensational subsidence flux form
    real(RP) :: qi_fxin(KA)     ! cmpensational subsidence flux form
    real(RP) :: qi_fxout(KA)    ! cmpensational subsidence flux form
    real(RP) :: qr_fxin(KA)     ! cmpensational subsidence flux form
    real(RP) :: qr_fxout(KA)    ! cmpensational subsidence flux form
    real(RP) :: qs_fxin(KA)     ! cmpensational subsidence flux form
    real(RP) :: qs_fxout(KA)    ! cmpensational subsidence flux form
    real(RP) :: rainfb(KA)      ! rain fall
    real(RP) :: snowfb(KA)      ! snow fall
    ! evariable moisture budeget
    real(RP) :: err      ! error (tmp var)
    real(RP) :: qinit    ! init water condensation (only qv)
    real(RP) :: qfinl    ! fineal water condensation (producted by cumulus parameterization)
    ! warm rain
    real(RP) :: cpm      ! tempolaly variable
    real(RP) :: UMF_tmp  ! tempolaly variable 
    real(RP) :: xcldfrac ! tempolaly variable
    !!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !! start cood
    !!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    !! initialize
    if(I_convflag == 2) return     ! noconvection
    ! start cood
    k_lcl = k_lcl_bf
    fbfrc = 0._RP
    if(I_convflag == 1) fbfrc = 1._RP
    !! set timecp(TIMEC)
    ! time scale of adjustment 
    vconv = 0.5_RP*(wspd(1) + wspd(2)) ! k_lcl + k_z5
    timecp = deltax/vconv
    time_advec = timecp ! advection time sclale
    ! 30 minutes < timecp < 60 minutes
    timecp = max(DEEPLIFETIME, timecp)
    timecp = min(3600._RP, timecp)
    if(I_convflag == 1) timecp = SHALLOWLIFETIME ! shallow convection timescale is 40 minutes
    nic = nint(timecp/dt)
    ! timecp is KF timescale
    ! must be bigger than "cumplus parameterization timestep " given by namelist
    timecp = real( nic,RP )*dt! determin timecp not change below
    !
    ! maximam of ainc calculate
    aincmx = 1000._RP
    k_lmax = max(k_lcl,k_lfs)
    do kk = k_lc,k_lmax
       if((upent(kk)-downent(kk)) > 1.e-3_RP) then
          ainctmp = ems(kk)/( (upent(kk) -downent(kk))*timecp )
          aincmx = min(aincmx,ainctmp)
       end if
    end do
    ainc = 1._RP
    if(aincmx < ainc ) ainc = aincmx
    !! for interpolation save original variable
    tder2 = tder
    prcp_flux2 = prcp_flux
    do kk = KS,k_top
       umf2(kk)     = umf(kk)
       upent2(kk)   = upent(kk)
       updet2(kk)   = updet(kk)
       qcdet2(kk)   = qcdet(kk)
       qidet2(kk)   = qidet(kk)
       dmf2(kk)     = dmf(kk)
       downent2(kk) = downent(kk)
       downdet2(kk) = downdet(kk)
    end do
    f_cape = 1._RP
    stab = 1.05_RP  -  DELCAPE ! default 0.95
    noiter = 0 ! noiter=1 then stop iteration
    istop = 0
    if (I_convflag == 1) then ! shallow convection
       ! refarence Kain 2004
       tkemax = 5._RP
       evac = 0.50_RP*tkemax*1.e-1_RP
       ainc = evac*dpthmx*deltax*deltax/( umflcl*GRAV*timecp)
       tder = tder2*ainc
       prcp_flux = prcp_flux2*ainc
       do kk = KS,k_top
          umf(kk)     = umf2(kk)*ainc
          upent(kk)   = upent2(kk)*ainc
          updet(kk)   = updet2(kk)*ainc
          qcdet(kk)   = qcdet2(kk)*ainc
          qidet(kk)   = qidet2(kk)*ainc
          dmf(kk)     = dmf2(kk)*ainc
          downent(kk) = downent2(kk)*ainc
          downdet(kk) = downdet2(kk)*ainc
       end do
    end if
    !! theta set up by Emanuel Atomospheric convection, 1994 111p
    !! original KF theta is calced  apploximatly.
    do kk = KS,k_top
       call calcexn(pres(kk),qv(kk),exn(kk))
       theta(kk) = temp_bf(kk)*exn(kk)
       call calcexn(pres(kk),qvdet(kk),exn(kk)) 
       theta_u(kk) = temp_u(kk)*exn(kk) 
    end do
    temp_g(k_top+1:KE) = temp_bf(k_top+1:KE)
    qv_g(k_top+1:KE)   = qv(k_top+1:KE)
    omg(:) = 0._RP ! initialize
    !!
    !!@@@  start main loop @@@
    do ncount = 1,10 ! itaration
       dtt = timecp
       do kk = KS,k_top
          domg_dp(kk) = -(upent(kk) - downent(kk) - updet(kk) - downdet(kk))*emsd(kk)
          if(kk > KS)then
             omg(kk) = omg(kk-1) - deltap(kk-1)*domg_dp(kk-1)
             absomg = abs(omg(kk))
             absomgtc  = absomg*timecp
             f_dp = 0.75*deltap(kk-1)
             if(absomgtc > f_dp)THEN
                dtt_tmp = f_dp/abs(omg(kk))
                dtt=min(dtt,dtt_tmp) ! it is use determin nstep
             end if
          end if
       end do
       !! theta_nw and qv_nw has valus only in 1 to k_top
       do kk = KS, k_top 
          theta_nw(kk) = theta(kk)
          qv_nw(kk)    = qv(kk)
          fxm(kk)      = omg(kk)*deltax*deltax/GRAV ! fluxmass
       end do
       nstep  = nint(timecp/dtt + 1) ! how many time step forwad 
       deltat = timecp/real(nstep,RP) ! deltat*nstep = timecp
       do ntimecount = 1, nstep
          ! tempolaly time forward start iteration 
          !!... assign theta and q variables at the top and bottom of each layer based
          !!... on sign of omega
          do kk = KS, k_top   ! initialize
             theta_fxin(kk)  = 0._RP
             theta_fxout(kk) = 0._RP
             qv_fxin(kk)     = 0._RP
             qv_fxout(kk)    = 0._RP
          end do
          do kk = KS+1,k_top ! calc flux variable
             if(omg(kk) <= 0._RP) then 
                theta_fxin(kk)      = -fxm(kk)*theta_nw(kk-1)
                qv_fxin(kk)         = -fxm(kk)*qv_nw(kk-1)
                theta_fxout(kk - 1) = theta_fxout(kk-1) + theta_fxin(kk)
                qv_fxout(kk - 1)    = qv_fxout(kk-1)    + qv_fxin(kk)
             else
                theta_fxout(kk)    = fxm(kk)*theta_nw(kk) 
                qv_fxout(kk)       = fxm(kk)*qv_nw(kk)
                theta_fxin(kk - 1) = theta_fxin(kk-1) + theta_fxout(kk)
                qv_fxin(kk - 1)    = qv_fxin(kk-1)    + qv_fxout(kk)
             end if
          end do
          !!   update the theta and qv variables at each level
          !!   only theta and qv calc cape below and if cape > 10%of old cape then iterate
          do kk = KS, k_top
             theta_nw(kk) = theta_nw(kk) &
                  + (theta_fxin(kk) - theta_fxout(kk)  &
                  + updet(kk)*theta_u(kk) + downdet(kk)*theta_d(kk)  &
                  - ( upent(kk) - downent(kk) )*theta(kk) ) *deltat*emsd(kk)
             qv_nw(kk) = qv_nw(kk) &
                  + (qv_fxin(kk) - qv_fxout(kk)  &
                  + updet(kk)*qvdet(kk) + downdet(kk)*qv_d(kk)  &
                  - ( upent(kk) - downent(kk) )*qv(kk) )*deltat*emsd(kk)
          end do
       end do ! ntimecount

       do kk = KS, k_top
          theta_g(kk) = theta_nw(kk)
          qv_g(kk)    = qv_nw(kk)
       end do
       ! Check to see if mixing ratio dips below zero anywere ; if so, borrow
       ! moisture forom adjacent layers to bring it back up abobe zero
       do kk = KS,k_top
          if(qv_g(kk) < 0._RP ) then ! negative moisture
             if(kk == KS) then
                write(*,*) "error qv<0 @ Kain-Fritsch cumulus parameterization"
                write(*,*) "@sub scale_atmos_phy_cp_kf",__FILE__,__LINE__
                call PRC_MPIstop
             end if
             kkp1 = kk + 1
             if(kk == k_top) then
                kkp1 = k_lcl
             end if
             tma        = qv_g(kkp1)*ems(kkp1)
             tmb        = qv_g(kk-1)*ems(kk-1)
             tmm        = (qv_g(kk) - 1.e-9_RP )*ems(kk)
             bcoeff     = -tmm/((tma*tma)/tmb + tmb)
             acoeff     = bcoeff*tma/tmb
             tmb        = tmb*(1._RP - bcoeff)
             tma        = tma*(1._RP - acoeff)
             qv_g(kk)   = 1.e-9_RP
             qv_g(kkp1) = tma*emsd(kkp1)
             qv_g(kk-1) = tmb*emsd(kk-1)
          end if
       end do
       ! clculate top layer omega and conpare determ omg
       topomg = (updet(k_top) - upent(k_top))*deltap(k_top)*emsd(k_top)
       if( abs(topomg - omg(k_top)) > 1.e-3_RP) then ! not same omega velocity error
          istop = 1
          write(*,*) "xxxERROR@KF omega is not consistent",ncount
          write(*,*) "omega error",abs(topomg - omg(k_top)),k_top,topomg,omg(k_top)
          call PRC_MPIstop
       end if
       ! convert theta to T
       do kk = KS,k_top
          call calcexn(pres(kk),qv_g(kk),exn(kk))
          temp_g(kk)  = theta_g(kk)/exn(kk)
          tempv_g(kk) = temp_g(kk)*(1._RP + 0.608_RP*qv_g(kk))
       end do
       !------------------------------------------------------------------------------
       ! compute new cloud and change in available bouyant energy(CAPE)
       !------------------------------------------------------------------------------
       if(I_convflag == 1) then
          exit ! if shallow convection , calc Cape  is not used
       end if
       temp_mix = 0._RP
       qv_mix   = 0._RP
       do kk = k_lc, k_pbl
          temp_mix = temp_mix + deltap(kk)*temp_g(kk)
          qv_mix   = qv_mix + deltap(kk)*qv_g(kk)
          !         presmix = presmix + deltap(kk)
       end do
       temp_mix = temp_mix/dpthmx
       qv_mix   = qv_mix/dpthmx
       ! calc saturate water vapor pressure
       call ATMOS_SATURATION_psat_liq(es,temp_mix)
       qvss = 0.622_RP*es/(presmix -es) ! saturate watervapor
       !!
       !!... Remove supersaturation for diagnostic purposes, if necessary..
       !!
       if (qv_mix > qvss) then ! saturate then
          RL       = XLV0 -XLV1*temp_mix
          CPM      = CP*(1._RP + 0.887_RP*qv_mix)
          DSSDT    = qvss*(CLIQ-BLIQ*DLIQ)/( (temp_mix-DLIQ)**2)
          DQ       = (qv_mix -qvss)/(1._RP + RL*DSSDT/CPM)
          temp_mix = temp_mix + RL/CP*DQ
          qv_mix   = qv_mix - DQ
          temp_lcl = temp_mix
       else !same as detern trigger
          qv_mix   = max(qv_mix,0._RP)
          emix     = qv_mix*presmix/(0.6222_RP + qv_mix)
          TLOG     = log(emix/ALIQ)
          ! dew point temperature Bolton(1980)
          TDPT     = (CLIQ - DLIQ*TLOG)/(BLIQ - TLOG)
          temp_lcl = TDPT - (0.212_RP + 1.57e-3_RP*(TDPT - TEM00) - 4.36e-4_RP*(temp_mix - TEM00))*(temp_mix - TDPT)
          temp_lcl = min(temp_lcl,temp_mix)
       end if
       tempv_lcl = temp_lcl*(1._RP + 0.608_RP*qv_mix)
       z_lcl     = zmix + (temp_lcl - temp_mix)/GdCP
       do kk = k_lc, KE
          k_lcl = kk
          if( z_lcl <= z_kf(kk) )  exit 
       end do
       ! estimate environmental tempeature and mixing ratio
       ! interpolate environment temperature and vapor at LCL
       k_lclm1   = k_lcl - 1
       deltaz    = ( z_lcl - z_kf(k_lclm1) )/( z_kf(k_lcl)- z_kf(k_lclm1 ) )
       temp_env  = temp_g(k_lclm1) + ( temp_g(k_lcl) - temp_g(k_lclm1) )*deltaz
       qv_env    = qv_g(k_lclm1) + ( qv_g(k_lcl) - qv_g(k_lclm1) )*deltaz
       tempv_env = temp_env*( 1._RP + 0.608_RP*qv_env )
       !!       pres_lcl=pres(k_lcl-1)+(pres(k_lcl-1)-pres(k_lcl-1))*deltaz
       theta_eu(k_lclm1)=temp_mix*(PRE00/presmix)**(0.2854_RP*(1._RP - 0.28_RP*qv_mix))*   &
            exp((3374.6525_RP/temp_lcl-2.5403_RP)*qv_mix*(1._RP + 0.81_RP*qv_mix))
       !!
       !!...COMPUTE ADJUSTED ABE(ABEG).(CAPE)
       !!
       cape_g = 0._RP !  cape "_g" add because 
       do kk=k_lclm1,k_top-1 ! LTOPM1
          kkp1=kk+1
          theta_eu(kkp1) = theta_eu(kk)
          call tpmix2dd(pres(kkp1),theta_eu(kkp1),temp_gu(kkp1),qv_gu(kkp1)) ! get temp_gu and qv_gu
          tempvq_u(kkp1) = temp_gu(kkp1)*(1._RP + 0.608_RP*qv_gu(kkp1) - qc(kkp1)- qi(kkp1))
          if(kk == k_lclm1) then !  interporate
             dzz = z_kf(k_lcl) - z_lcl
             dilbe = ((tempv_lcl + tempvq_u(kkp1))/(tempv_env + tempv_g(kkp1)) - 1._RP)*dzz
          else
             dzz = dz_kf(kk)
             dilbe = ((tempvq_u(kk) + tempvq_u(kkp1))/(tempv_g(kk) + tempv_g(kkp1)) - 1._RP)*dzz
          end if
          if(dilbe > 0._RP) cape_g = cape_g + dilbe*GRAV
          !!
          !!...DILUTE BY ENTRAINMENT BY THE RATE AS ORIGINAL UPDRAFT...
          !!
          call ENVIRTHT(pres(kkp1),temp_g(kkp1),qv_g(kkp1),theta_eg(kkp1)) ! calc get theta_eg
          !! theta_eg(environment theta_E )
          !          theta_eu(kkp1) = theta_eu(kkp1)*(1._RP/umfnewdold(kkp1)) + theta_eg(kkp1)*(1._RP - (1._RP/umfnewdold(kkp1)))
          theta_eu(kkp1) = theta_eu(kkp1)*(umfnewdold(kkp1)) + theta_eg(kkp1)*(1._RP - (umfnewdold(kkp1)))
       end do
       if (noiter == 1) exit ! noiteration
       dcape = max(cape - cape_g,cape*0.1_RP) ! delta cape
       f_cape = cape_g/cape ! ratio of cape new/old
       if(f_cape > 1._RP .and. I_convflag == 0) then ! if deep convection and cape is inclease this loop
          I_convflag = 2
          return
       end if
       if(ncount /= 1) then
          if(abs(ainc - aincold) < 1.e-4_RP) then !  IN cycle not change facter then exit iter next loop 
             noiter = 1 ! exit this loop in nex step
             ainc   = aincold
             cycle ! iter
          end if
          dfda = (f_cape - f_capeold)/(ainc - aincold)
          if (dfda > 0._RP ) then
             noiter = 1 ! exit this loop @next loop step
             ainc   = aincold
             cycle ! iter
          end if
       end if
       aincold   = ainc
       f_capeold = f_cape
       !! loop exit
       !! aincmx is upper limit of massflux factor
       !! if massflux factor 'ainc' is near "aincmax" then exit
       !! but need   CAPE is less than 90% of original
       if (ainc/aincmx > 0.999 .and. f_cape > 1.05-stab ) then
          exit
       end if
       !! loop exit 1. or 2.
       !! 1. NEW cape is less than 10% of oliginal cape
       !! 2. ncount = 10
       if( (f_cape <=  1.05-stab .and. f_cape >= 0.95-stab) .or. ncount == 10) then
          exit
       else ! no exit 
          if(ncount > 10) exit ! sayfty ??  ncount musn't over 10...
          if(f_cape == 0._RP) then ! f_cape is 0 -> new cape is 0 : too reducted
             ainc = ainc*0.5_RP
          else
             if(dcape < 1.e-4) then ! too small cape then exit at next loop(iter loop) step
                noiter = 1
                ainc   = aincold
                cycle
             else ! calculate new factor  ainc
                ainc = ainc*stab*cape/dcape
             end if
          end if
          ainc = min(aincmx,ainc) ! ainc must be less than aincmx
          !!  if ainc becomes very small, effects of convection will be minimal so just ignore it
          if (ainc < 0.05) then !! noconvection 
             I_convflag = 2 
             return
          end if
          !!update valuables use factar ainc
          tder      = tder2*ainc
          prcp_flux = prcp_flux2*ainc
          do kk = KS,k_top
             umf(kk)     = umf2(kk)*ainc
             upent(kk)   = upent2(kk)*ainc
             updet(kk)   = updet2(kk)*ainc
             qcdet(kk)   = qcdet2(kk)*ainc
             qidet(kk)   = qidet2(kk)*ainc
             dmf(kk)     = dmf2(kk)*ainc
             downent(kk) = downent2(kk)*ainc
             downdet(kk) = downdet2(kk)*ainc
          end do
          ! go back up for another iteration
       end if

    end do ! iter(ncount)
    ! get the cloud fraction
    cldfrac_KF(:,:) = 0._RP
    if (I_convflag == 1) then
       do kk = k_lcl-1, k_top
          umf_tmp = umf(kk)/(deltax*deltax)
          xcldfrac = 0.07*log(1._RP+(500._RP*UMF_tmp))
          xcldfrac = max(1.e-2_RP,xcldfrac)
          cldfrac_KF(kk,1) = min(2.e-2_RP,xcldfrac) ! shallow
       end do
    else
       do kk = k_lcl-1, k_top
          umf_tmp = umf(kk)/(deltax*deltax)
          xcldfrac = 0.14*log(1._RP+(500._RP*UMF_tmp))
          xcldfrac = max(1.e-2_RP,xcldfrac)
          cldfrac_KF(kk,2) = min(6.e-1_RP,xcldfrac) ! deep
       end do
    end if
    !!
    !!...compute hydrometeor tendencies as is done for T,qv...
    !!
    !! ...frc2 is the fraction of total condensate
    !!... generated that goes into precipitation
    !!
    !! Redistribute hydrometerors according to the final mass-flux variables
    !! dtheta/dt and dqv/dt is
    if (cpr > 0._RP) then
       frc2 = prcp_flux/(cpr*ainc)
    else
       frc2 = 0._RP
    end if
    !! no qc qi qr qs inputted in KF scheme
    qc_nw(KS:kE) = 0._RP
    qi_nw(KS:kE) = 0._RP
    qr_nw(KS:kE) = 0._RP
    qs_nw(KS:kE) = 0._RP
    do kk = KS,k_top
       rainfb(kk) = flux_qr(kk)*ainc*fbfrc*frc2
       snowfb(kk) = flux_qs(kk)*ainc*fbfrc*frc2
    end do

    do ntimecount = 1, nstep ! same as T, QV
       do kk = KS, k_top   ! initialize
          qc_fxin(kk)  = 0._RP
          qc_fxout(kk) = 0._RP
          qi_fxin(kk)  = 0._RP
          qi_fxout(kk) = 0._RP
          qr_fxin(kk)  = 0._RP
          qr_fxout(kk) = 0._RP
          qs_fxin(kk)  = 0._RP
          qs_fxout(kk) = 0._RP
       end do
       do kk = KS+1,k_top ! calc flux variable
          if(omg(kk) <= 0._RP) then
             qc_fxin(kk) = -fxm(kk)*qc_nw(kk-1)
             qi_fxin(kk) = -fxm(kk)*qi_nw(kk-1)
             qr_fxin(kk) = -fxm(kk)*qr_nw(kk-1)
             qs_fxin(kk) = -fxm(kk)*qs_nw(kk-1)
             !
             qc_fxout(kk-1) = qc_fxout(kk-1) + qc_fxin(kk)
             qi_fxout(kk-1) = qi_fxout(kk-1) + qi_fxin(kk)
             qr_fxout(kk-1) = qr_fxout(kk-1) + qr_fxin(kk)
             qs_fxout(kk-1) = qs_fxout(kk-1) + qs_fxin(kk)
          else
             qc_fxout(kk) = fxm(kk)*qc_nw(kk)
             qi_fxout(kk) = fxm(kk)*qi_nw(kk)
             qr_fxout(kk) = fxm(kk)*qr_nw(kk)
             qs_fxout(kk) = fxm(kk)*qs_nw(kk)
             !
             qc_fxin(kk-1) =  qc_fxin(kk-1) + qc_fxout(kk)
             qi_fxin(kk-1) =  qi_fxin(kk-1) + qi_fxout(kk)
             qr_fxin(kk-1) =  qr_fxin(kk-1) + qr_fxout(kk)
             qs_fxin(kk-1) =  qs_fxin(kk-1) + qs_fxout(kk)
          end if
       end do

       do kk = KS, k_top
          qc_nw(kk) = qc_nw(kk) + (qc_fxin(kk) - qc_fxout(kk) + qcdet(kk) )*deltat*emsd(kk)
          qi_nw(kk) = qi_nw(kk) + (qi_fxin(kk) - qi_fxout(kk) + qidet(kk) )*deltat*emsd(kk)
          qr_nw(kk) = qr_nw(kk) + (qr_fxin(kk) - qr_fxout(kk) + rainfb(kk) )*deltat*emsd(kk)
          qs_nw(kk) = qs_nw(kk) + (qs_fxin(kk) - qs_fxout(kk) + snowfb(kk) )*deltat*emsd(kk)
       end do
       ! in original qlg= qlpa but it is not nessesary
    end do

    !! cumulus parameterization rain (rainrate_cp)and rain rate (rainratecp) is detern
    rainrate_cp =  prcp_flux*(1._RP - fbfrc)/(deltax*deltax) ! if shallow convection then fbfrc = 1. -> noprcpitation
    !! evaluate moisuture budget
    qinit      = 0._RP ! initial qv 
    qfinl      = 0._RP ! final water subsidence qv,qi,qr,qs...
    dpth       = 0._RP  !
    do kk = KS,k_top
       dpth  = dpth + deltap(kk)
       qinit = qinit + qv(kk)*ems(kk)
       qfinl = qfinl + qv_g(kk)*ems(kk) ! qv
       qfinl = qfinl + (qc_nw(kk) + qi_nw(kk) + qr_nw(kk) + qs_nw(kk))*ems(kk)
    end do
    qfinl = qfinl + prcp_flux*timecp*(1._RP - fbfrc)
    err   = (qfinl - qinit )*100._RP/qinit
    if (abs(err) > 0.05_RP .and. istop == 0) then
       ! write error message
       ! moisture budget error
       istop = 1
       write(*,*) "XXXX ERROR@KF,MOISTURE"
       call PRC_MPIstop
    end if
    !! feed back to resolvable scale tendencies
    !! if the advective time period(time_advec) is less than specified minimum
    !! timec, allow feed back to occur only during time_advec
    !! calc tendency of Theta_,qv, qc, qi, qr, qs
    !! temp_g ,qv_g, qc_nw,qi_nw,qr_nw,qs_nw
    !! calc new theta (SCALE theta) 
    !! then return tendency
    !! warm rain
    if (WARMRAIN) then ! assume Kessuler scheem?
       !!
       !...IF ICE PHASE IS NOT ALLOWED, MELT ALL FROZEN HYDROMETEORS...         
       !!
       do kk = KS,KE
          cpm        = CP*(1._RP + 0.887_RP*qv_g(kk))
          temp_g(kk) = temp_g(kk) - (qi_nw(kk) + qs_nw(kk))*EMELT/CPM
          qc_nw(kk)  = qc_nw(kk) + qi_nw(kk)
          qr_nw(kk)  = qr_nw(kk) + qs_nw(kk)
          qi_nw(kk)  = 0._RP
          qs_nw(kk)  = 0._RP
       end do
       return
    elseif(.not. FLAG_QS ) then
       !!
       !!...IF ICE PHASE IS ALLOWED, BUT MIXED PHASE IS NOT, MELT FROZEN HYDROME   
       !!...BELOW THE MELTING LEVEL, FREEZE LIQUID WATER ABOVE THE MELTING LEVEL  
       !!
       do kk = KS,KE
          cpm = cp*(1._RP + 0.887*qv_g(kk))
          if(kk < k_ml) then
             temp_g(kk) = temp_g(kk) - (qi_nw(kk) + qs_nw(kk))*EMELT/CPM
          elseif(kk > k_ml) then! kk == k_ml no melt
             temp_g(kk) = temp_g(kk) + (qi_nw(kk) + qs_nw(kk))*EMELT/CPM
          end if
          qc_nw(kk) = qc_nw(kk) + qi_nw(kk)
          qr_nw(kk) = qr_nw(kk) + qs_nw(kk)
          qi_nw(kk) = 0._RP
          qs_nw(kk) = 0._RP
       end do
       return
       !!
       !!...IF MIXED PHASE HYDROMETEORS ARE ALLOWED, FEED BACK CONVECTIVE TENDENCIES
       !!...OF HYDROMETEORS DIRECTLY...    
       !!
    elseif( FLAG_QS ) then
       if(.not.FLAG_QI) then
          do kk = KS, KE
             qs_nw(kk) =  qs_nw(kk) + qi_nw(kk)
          end do
       end if
       return
    else ! not allow
       write(*,*) "xxx ERROR@KF,NOTallow namelist"
       write(*,*) "!!!!!Moiture setting error in KF check namelist"
       call PRC_MPIstop
    end if
    return
  end subroutine CP_kf_compensational
  !! calc exner  function for potential temperature(contain water vapor)
  subroutine calcexn( pres,qv,exn) ! emanuel 1994 111pp ! check potential temperature definition
    use scale_const,only :&
         CP_dry => CONST_CPdry , &
         PRE00 => CONST_PRE00  , &
         R_dry => CONST_Rdry   
    implicit none
    real(RP),intent(in) :: pres
    real(RP),intent(in) :: qv
    real(RP),intent(out) :: exn
    exn = (PRE00/pres)**(0.2854_RP*(1._RP - 0.28_RP*qv ))
    !! exect
    !    exn = (PRE00/pres)**(( R_dry + qv*R_vap)/(Cp_dry + qv*Cp_vap))
    return
  end subroutine calcexn
!@!------------------------------------------------------------------------------  
!@! kf subroutine from WRF cood
!@!------------------------------------------------------------------------------
  subroutine precipitation_OC1973( &
       QLIQ,QICE,WTW,DZ,BOTERM,ENTERM,QNEWLQ,           &
       QNEWIC,QLQOUT,QICOUT,G)
    !-----------------------------------------------------------------------
    implicit none
    !-----------------------------------------------------------------------
    !  9/18/88...THIS PRECIPITATION FALLOUT SCHEME IS BASED ON THE SCHEME US
    !  BY OGURA AND CHO (1973).  LIQUID WATER FALLOUT FROM A PARCEL IS CAL-
    !  CULATED USING THE EQUATION DQ=-RATE*Q*DT, BUT TO SIMULATE A QUASI-
    !  CONTINUOUS PROCESS, AND TO ELIMINATE A DEPENDENCY ON VERTICAL
    !  RESOLUTION THIS IS EXPRESSED AS Q=Q*EXP(-RATE*DZ).

    real(RP), intent(in)    :: G
    real(RP), intent(in)    :: DZ,BOTERM,ENTERM!,RATE to be local variablebles
    real(RP), intent(inout) :: QLQOUT,QICOUT,WTW,QLIQ,QICE,QNEWLQ,QNEWIC

    real(RP) :: QTOT,QNEW,QEST,G1,WAVG,CONV,RATIO3,OLDQ,RATIO4,DQ,PPTDRG

    QTOT=QLIQ+QICE
    QNEW=QNEWLQ+QNEWIC
    !
    !  ESTIMATE THE VERTICAL VELOCITY SO THAT AN AVERAGE VERTICAL VELOCITY 
    !  BE CALCULATED TO ESTIMATE THE TIME REQUIRED FOR ASCENT BETWEEN MODEL 
    !  LEVELS...
    !
    QEST=0.5_RP*(QTOT+QNEW)
    G1=WTW+BOTERM-ENTERM-2._RP*G*DZ*QEST/1.5_RP
    IF(G1.LT.0.0)G1=0._RP
    WAVG=0.5_RP*(SQRT(WTW)+SQRT(G1))
    CONV=RATE*DZ/WAVG               ! KF90  Eq. 9

    !
    !  RATIO3 IS THE FRACTION OF LIQUID WATER IN FRESH CONDENSATE, RATIO4 IS
    !  THE FRACTION OF LIQUID WATER IN THE TOTAL AMOUNT OF CONDENSATE INVOLV
    !  IN THE PRECIPITATION PROCESS - NOTE THAT ONLY 60% OF THE FRESH CONDEN
    !  SATE IS IS ALLOWED TO PARTICIPATE IN THE CONVERSION PROCESS...       
    !
    RATIO3=QNEWLQ/(QNEW+1.E-8_RP)
    !     OLDQ=QTOT
    QTOT=QTOT+0.6_RP*QNEW
    OLDQ=QTOT
    RATIO4=(0.6_RP*QNEWLQ+QLIQ)/(QTOT+1.E-8_RP)
    QTOT=QTOT*EXP(-CONV)            ! KF90  Eq. 9
    !
    !  DETERMINE THE AMOUNT OF PRECIPITATION THAT FALLS OUT OF THE UPDRAFT  
    !  PARCEL AT THIS LEVEL...
    !
    DQ=OLDQ-QTOT
    QLQOUT=RATIO4*DQ
    QICOUT=(1._RP-RATIO4)*DQ
    !
    !  ESTIMATE THE MEAN LOAD OF CONDENSATE ON THE UPDRAFT IN THE LAYER, CAL
    !  LATE VERTICAL VELOCITY
    !
    PPTDRG=0.5_RP*(OLDQ+QTOT-0.2_RP*QNEW)
    WTW=WTW+BOTERM-ENTERM-2._RP*G*DZ*PPTDRG/1.5_RP
    IF(ABS(WTW).LT.1.E-4_RP)WTW=1.E-4_RP
    !
    !  DETERMINE THE NEW LIQUID WATER AND ICE CONCENTRATIONS INCLUDING LOSSE
    !  DUE TO PRECIPITATION AND GAINS FROM CONDENSATION...
    !
    QLIQ=RATIO4*QTOT+RATIO3*0.4_RP*QNEW
    QICE=(1._RP-RATIO4)*QTOT+(1._RP-RATIO3)*0.4_RP*QNEW
    QNEWLQ=0._RP
    QNEWIC=0._RP
    return
  end subroutine precipitation_OC1973

  ! Kessler type auto conversion
  subroutine precipitation_Kessler( &
       QLIQ,QICE,WTW,DZ,BOTERM,ENTERM,QNEWLQ,           &
       QNEWIC,QLQOUT,QICOUT,G)
     implicit none
    ! kesseler type
    ! auto conversion
    real(RP), intent(in  )  :: G
    real(RP), intent(in  )  :: DZ,BOTERM,ENTERM!,RATE to be local variablebles
    real(RP), intent(inout) :: QLQOUT,QICOUT,WTW,QLIQ,QICE,QNEWLQ,QNEWIC
    real(RP) :: pptdrg
    real(RP) :: total_liq, total_ice
    ! parameter module value kf_threshold
    real(RP) :: auto_qc, auto_qi
    auto_qc = kf_threshold
    auto_qi = kf_threshold

    total_liq = QLIQ + QNEWLQ
    total_ice = QICE + QNEWIC

    ! condensate in convective updraft is converted into precipitation
    qlqout = max( total_liq - auto_qc, 0.0_RP )
    qicout = max( total_ice - auto_qi, 0.0_RP )

    pptdrg = max( 0.5_RP * ( total_liq + total_ice - qlqout - qicout ), 0.0_RP )
    WTW=WTW+BOTERM-ENTERM-2._RP*G*DZ*PPTDRG/1.5_RP
    IF(ABS(WTW).LT.1.E-4_RP)WTW=1.E-4_RP

    QLIQ = max( total_liq - qlqout, 0.0_RP )
    QLQOUT = total_liq - QLIQ

    QICE = max( total_ice - qicout, 0.0_RP )
    QICOUT = total_ice - QICE

    QNEWLQ=0.0_RP
    QNEWIC=0.0_RP

    return
  end subroutine precipitation_Kessler

  ! ----------------------------------------------------------------------
  subroutine TPMIX2(p,thes,tu,qu,qliq,qice,qnewlq,qnewic )!,XLV1,XLV0)
    !
    ! Lookup table variables:
    !     INTEGER, PARAMETER :: (KFNT=250,KFNP=220)
    !     REAL, SAVE, DIMENSION(1:KFNT,1:KFNP) :: TTAB,QSTAB
    !     REAL, SAVE, DIMENSION(1:KFNP) :: THE0K
    !     REAL, SAVE, DIMENSION(1:200) :: ALU
    !     REAL, SAVE :: RDPR,RDTHK,PLUTOP
    ! End of Lookup table variables:
    !-----------------------------------------------------------------------
    implicit none
    !-----------------------------------------------------------------------
    real(RP), intent(in)    :: P,THES!,XLV1,XLV0
    real(RP), intent(out)   :: QNEWLQ,QNEWIC
    real(RP), intent(inout) :: TU,QU,QLIQ,QICE

    real(RP) :: TP,QQ,BTH,TTH,PP,T00,T10,T01,T11,Q00,Q10,Q01,Q11
    real(RP) :: TEMP,QS,QNEW,DQ,QTOT,RLL,CPP
    integer  :: IPTB,ITHTB
    !-----------------------------------------------------------------------

    !c******** LOOKUP TABLE VARIABLES... ****************************
    !      parameter(kfnt=250,kfnp=220)
    !c
    !      COMMON/KFLUT/ ttab(kfnt,kfnp),qstab(kfnt,kfnp),the0k(kfnp),
    !     *              alu(200),rdpr,rdthk,plutop 
    !C*************************************************************** 
    !c
    !c***********************************************************************
    !c     scaling pressure and tt table index
    !c***********************************************************************
    !c
    tp=(p-plutop)*rdpr
    qq=tp-aint(tp)
    iptb=int(tp)+1

    !
    !***********************************************************************
    !              base and scaling factor for the
    !***********************************************************************
    !
    !  scaling the and tt table index
    bth=(the0k(iptb+1)-the0k(iptb))*qq+the0k(iptb)
    tth=(thes-bth)*rdthk
    pp   =tth-aint(tth)
    ithtb=int(tth)+1
    !      IF(IPTB.GE.220 .OR. IPTB.LE.1 .OR. ITHTB.GE.250 .OR. ITHTB.LE.1)THEN
    IF(IPTB.GE.kfnp .OR. IPTB.LE.1 .OR. ITHTB.GE.250 .OR. ITHTB.LE.1)THEN
       ! modify
       write(*,*)   '**** OUT OF BOUNDS *********',IPTB,ITHTB,P,THES
       !        flush(98)
    ENDIF
    !
    t00=ttab(ithtb  ,iptb  )
    t10=ttab(ithtb+1,iptb  )
    t01=ttab(ithtb  ,iptb+1)
    t11=ttab(ithtb+1,iptb+1)
    !
    q00=qstab(ithtb  ,iptb  )
    q10=qstab(ithtb+1,iptb  )
    q01=qstab(ithtb  ,iptb+1)
    q11=qstab(ithtb+1,iptb+1)
    !
    !***********************************************************************
    !              parcel temperature
    !***********************************************************************
    !
    temp=(t00+(t10-t00)*pp+(t01-t00)*qq+(t00-t10-t01+t11)*pp*qq)
    !
    qs=(q00+(q10-q00)*pp+(q01-q00)*qq+(q00-q10-q01+q11)*pp*qq)
    !
    DQ=QS-QU
    IF(DQ.LE.0._RP)THEN
       QNEW=QU-QS
       QU=QS
    ELSE 
       !
       !   IF THE PARCEL IS SUBSATURATED, TEMPERATURE AND MIXING RATIO MUST BE
       !   ADJUSTED...IF LIQUID WATER IS PRESENT, IT IS ALLOWED TO EVAPORATE
       ! 
       QNEW=0._RP
       QTOT=QLIQ+QICE
       !
       !   IF THERE IS ENOUGH LIQUID OR ICE TO SATURATE THE PARCEL, TEMP STAYS AT ITS
       !   WET BULB VALUE, VAPOR MIXING RATIO IS AT SATURATED LEVEL, AND THE MIXING
       !   RATIOS OF LIQUID AND ICE ARE ADJUSTED TO MAKE UP THE ORIGINAL SATURATION
       !   DEFICIT... OTHERWISE, ANY AVAILABLE LIQ OR ICE VAPORIZES AND APPROPRIATE
       !   ADJUSTMENTS TO PARCEL TEMP; VAPOR, LIQUID, AND ICE MIXING RATIOS ARE MADE.
       !
       !...subsaturated values only occur in calculations involving various mixtures of
       !...updraft and environmental air for estimation of entrainment and detrainment.
       !...For these purposes, assume that reasonable estimates can be given using 
       !...liquid water saturation calculations only - i.e., ignore the effect of the
       !...ice phase in this process only...will not affect conservative properties...
       !
       IF(QTOT.GE.DQ)THEN
          qliq=qliq-dq*qliq/(qtot+1.e-10_RP)
          qice=qice-dq*qice/(qtot+1.e-10_RP)
          QU=QS
       ELSE
          RLL=XLV0-XLV1*TEMP
          CPP=1004.5_RP*(1._RP+0.89_RP*QU)
          IF(QTOT.LT.1.E-10_RP)THEN
             !
             !...IF NO LIQUID WATER OR ICE IS AVAILABLE, TEMPERATURE IS GIVEN BY:
             TEMP=TEMP+RLL*(DQ/(1._RP+DQ))/CPP
          ELSE
             !
             !...IF SOME LIQ WATER/ICE IS AVAILABLE, BUT NOT ENOUGH TO ACHIEVE SATURATION,
             !   THE TEMPERATURE IS GIVEN BY:
             !
             TEMP=TEMP+RLL*((DQ-QTOT)/(1._RP+DQ-QTOT))/CPP
             QU=QU+QTOT
             QTOT=0._RP
             QLIQ=0._RP
             QICE=0._RP
          ENDIF
       ENDIF
    ENDIF
    TU=TEMP
    qnewlq=qnew
    qnewic=0._RP
    return
  end subroutine TPMIX2
  !******************************************************************************
  subroutine DTFRZNEW(TU,P,THTEU,QU,QFRZ,QICE)!,ALIQ,BLIQ,CLIQ,DLIQ)
    !-----------------------------------------------------------------------
    use scale_precision
    use scale_atmos_saturation ,only :&
         ATMOS_SATURATION_psat_liq
    implicit none
    !-----------------------------------------------------------------------
    real(RP), intent(in)    :: P,QFRZ!to module variable,ALIQ,BLIQ,CLIQ,DLIQ
    real(RP), intent(inout) :: TU,THTEU,QU,QICE

    real(RP) :: RLC,RLS,RLF,CPP,A,DTFRZ,ES,QS,DQEVAP,PII
    !-----------------------------------------------------------------------
    !
    !...ALLOW THE FREEZING OF LIQUID WATER IN THE UPDRAFT TO PROCEED AS AN 
    !...APPROXIMATELY LINEAR FUNCTION OF TEMPERATURE IN THE TEMPERATURE RANGE 
    !...TTFRZ TO TBFRZ...
    !...FOR COLDER TEMPERATURES, FREEZE ALL LIQUID WATER...
    !...THERMODYNAMIC PROPERTIES ARE STILL CALCULATED WITH RESPECT TO LIQUID WATER
    !...TO ALLOW THE USE OF LOOKUP TABLE TO EXTRACT TMP FROM THETAE...
    !

    RLC=2.5E6_RP-2369.276_RP*(TU-273.16_RP)
    !      RLC=2.5E6_RP-2369.276_RP*(TU-273.15_RP)   ! 273.16 -> 273.15 ??
    RLS=2833922._RP-259.532_RP*(TU-273.16_RP)
    !      RLS=2833922._RP-259.532_RP*(TU-273.15_RP) ! 273.16 -> 273.15 ??
    RLF=RLS-RLC
    CPP=1004.5_RP*(1._RP+0.89_RP*QU)
    !
    !  A = D(es)/DT IS THAT CALCULATED FROM BUCK (1981) EMPERICAL FORMULAS
    !  FOR SATURATION VAPOR PRESSURE...
    !
    A=(CLIQ-BLIQ*DLIQ)/((TU-DLIQ)*(TU-DLIQ))
    DTFRZ = RLF*QFRZ/(CPP+RLS*QU*A)
    TU = TU+DTFRZ
    !      call ATMOS_SATURATION_psat_liq(ES,TU) !saturation vapar pressure

    ES = ALIQ*EXP((BLIQ*TU-CLIQ)/(TU-DLIQ))
    QS = ES*0.622_RP/(P-ES)
    !
    !...FREEZING WARMS THE AIR AND IT BECOMES UNSATURATED...ASSUME THAT SOME OF THE 
    !...LIQUID WATER THAT IS AVAILABLE FOR FREEZING EVAPORATES TO MAINTAIN SATURA-
    !...TION...SINCE THIS WATER HAS ALREADY BEEN TRANSFERRED TO THE ICE CATEGORY,
    !...SUBTRACT IT FROM ICE CONCENTRATION, THEN SET UPDRAFT MIXING RATIO AT THE NEW
    !...TEMPERATURE TO THE SATURATION VARIABLE...
    !
    DQEVAP = QS-QU
    QICE = QICE-DQEVAP
    QU = QU+DQEVAP
    PII=(1.E5_RP/P)**(0.2854_RP*(1._RP-0.28_RP*QU))
    ! Bolton 1980
    ! Emanuel 1994 132p eq(4.7.9) pseudoequivalent PT
    THTEU = TU*PII*EXP((3374.6525_RP/TU - 2.5403_RP)*QU*(1._RP + 0.81_RP*QU))
    !
  end subroutine DTFRZNEW
  ! --------------------------------------------------------------------------------
  subroutine PROF5(EQ,EE,UD)
    !
    !***********************************************************************
    !*****    GAUSSIAN TYPE MIXING PROFILE....******************************
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !  THIS SUBROUTINE INTEGRATES THE AREA UNDER THE CURVE IN THE GAUSSIAN  
    !  DISTRIBUTION...THE NUMERICAL APPROXIMATION TO THE INTEGRAL IS TAKEN FROM
    !  "HANDBOOK OF MATHEMATICAL FUNCTIONS WITH FORMULAS, GRAPHS AND MATHEMATICS TABLES"
    !  ED. BY ABRAMOWITZ AND STEGUN, NATL BUREAU OF STANDARDS APPLIED
    !  MATHEMATICS SERIES.  JUNE, 1964., MAY, 1968.
    !                                     JACK KAIN
    !                                     7/6/89
    !  Solves for KF90 Eq. 2
    !
    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    !-----------------------------------------------------------------------
    implicit none
    !-----------------------------------------------------------------------
    real(RP), intent(in)    :: EQ
    real(RP), intent(inout) :: EE,UD

    real(RP) :: SQRT2P,A1,A2,A3,P,SIGMA,FE
    real(RP) :: X,Y,EY,E45,T1,T2,C1,C2

    DATA SQRT2P,A1,A2,A3,P,SIGMA,FE/2.506628_RP,0.4361836_RP,-0.1201676_RP,       &
         0.9372980_RP,0.33267_RP,0.166666667_RP,0.202765151_RP/
    X = (EQ - 0.5_RP)/SIGMA
    Y = 6._RP*EQ - 3._RP
    EY = EXP(Y*Y/(-2._RP))
    E45 = EXP(-4.5_RP)
    T2 = 1._RP/(1._RP + P*ABS(Y))
    T1 = 0.500498_RP
    C1 = A1*T1+A2*T1*T1+A3*T1*T1*T1
    C2 = A1*T2+A2*T2*T2+A3*T2*T2*T2
    IF(Y.GE.0._RP)THEN
       EE=SIGMA*(0.5_RP*(SQRT2P-E45*C1-EY*C2)+SIGMA*(E45-EY))-E45*EQ*EQ/2._RP
       UD=SIGMA*(0.5_RP*(EY*C2-E45*C1)+SIGMA*(E45-EY))-E45*(0.5_RP+EQ*EQ/2._RP-    &
            EQ)
    ELSE
       EE=SIGMA*(0.5_RP*(EY*C2-E45*C1)+SIGMA*(E45-EY))-E45*EQ*EQ/2._RP
       UD=SIGMA*(0.5_RP*(SQRT2P-E45*C1-EY*C2)+SIGMA*(E45-EY))-E45*(0.5_RP+EQ*   &
            EQ/2._RP-EQ)
    ENDIF
    EE=EE/FE
    UD=UD/FE
  end subroutine PROF5
  ! ------------------------------------------------------------------------
  subroutine TPMIX2DD(p,thes,ts,qs)!,i,j)
    !
    ! Lookup table variables:
    !     integer, PARAMETER :: (KFNT=250,KFNP=220)
    !     REAL, SAVE, DIMENSION(1:KFNT,1:KFNP) :: TTAB,QSTAB
    !     REAL, SAVE, DIMENSION(1:KFNP) :: THE0K
    !     REAL, SAVE, DIMENSION(1:200) :: ALU
    !     REAL, SAVE :: RDPR,RDTHK,PLUTOP
    ! End of Lookup table variables:
    !-----------------------------------------------------------------------
    implicit none
    !-----------------------------------------------------------------------
    real(RP), intent(in)    :: P,THES
    real(RP), intent(inout) :: TS,QS

    real(RP) :: TP,QQ,BTH,TTH,PP,T00,T10,T01,T11,Q00,Q10,Q01,Q11
    integer  :: IPTB,ITHTB
    !-----------------------------------------------------------------------

    !
    !******** LOOKUP TABLE VARIABLES (F77 format)... ****************************
    !     parameter(kfnt=250,kfnp=220)
    !
    !     COMMON/KFLUT/ ttab(kfnt,kfnp),qstab(kfnt,kfnp),the0k(kfnp),        &
    !                   alu(200),rdpr,rdthk,plutop
    !***************************************************************
    !
    !***********************************************************************
    !     scaling pressure and tt table index
    !***********************************************************************
    !
    tp=(p-plutop)*rdpr
    qq=tp-aint(tp)
    iptb=int(tp)+1
    !
    !***********************************************************************
    !              base and scaling factor for the
    !***********************************************************************
    !
    !  scaling the and tt table index
    bth=(the0k(iptb+1)-the0k(iptb))*qq+the0k(iptb)
    tth=(thes-bth)*rdthk
    pp   =tth-aint(tth)
    ithtb=int(tth)+1
    !
    t00=ttab(ithtb  ,iptb  )
    t10=ttab(ithtb+1,iptb  )
    t01=ttab(ithtb  ,iptb+1)
    t11=ttab(ithtb+1,iptb+1)
    !
    q00=qstab(ithtb  ,iptb  )
    q10=qstab(ithtb+1,iptb  )
    q01=qstab(ithtb  ,iptb+1)
    q11=qstab(ithtb+1,iptb+1)
    !
    !***********************************************************************
    !              parcel temperature and saturation mixing ratio
    !***********************************************************************
    !
    ts=(t00+(t10-t00)*pp+(t01-t00)*qq+(t00-t10-t01+t11)*pp*qq)
    !
    qs=(q00+(q10-q00)*pp+(q01-q00)*qq+(q00-q10-q01+q11)*pp*qq)
    !
    return
  end subroutine TPMIX2DD
  !!-----------------------------------------------------------------------
  subroutine ENVIRTHT(P1,T1,Q1,THT1)!,ALIQ,BLIQ,CLIQ,DLIQ)
    !
    !-----------------------------------------------------------------------
    use scale_precision
    use scale_const,only : &
         P00 => CONST_PRE00
    !!         C1  ->
    implicit none
    !-----------------------------------------------------------------------
    real(RP), intent(in)  :: P1,T1,Q1!,ALIQ,BLIQ,CLIQ,DLIQ module variables 
    real(RP), intent(out) :: THT1

    real(RP) :: EE,TLOG,ASTRT,AINC,A1,TP,VALUE,AINTRP,TDPT,TSAT,THT
!    real(RP) :: T00,P00,C1,C2,C3,C4,C5
    real(RP),parameter :: C1=3374.6525_RP
    real(RP),parameter :: C2=2.5403_RP
    !-----------------------------------------------------------------------
    !   DATA T00,P00,C1,C2,C3,C4,C5/273.16_RP,1.E5_RP,3374.6525_RP,2.5403_RP,3114.834_RP,   &
    !        0.278296_RP,1.0723E-3_RP/
    !
    !  CALCULATE ENVIRONMENTAL EQUIVALENT POTENTIAL TEMPERATURE...
    !
    ! NOTE: Calculations for mixed/ice phase no longer used...jsk 8/00
    !        For example, KF90 Eq. 10 no longer used
    !
    EE=Q1*P1/(0.622_RP+Q1)
    !     TLOG=ALOG(EE/ALIQ)
    ! ...calculate LOG term using lookup table...
    !
!      astrt=1.e-3_RP
!      ainc=0.075_RP
!      a1=ee/aliq
!      tp=(a1-astrt)/ainc
!      indlu=int(tp)+1
!      value=(indlu-1)*ainc+astrt
!      aintrp=(a1-value)/ainc
    !      tlog=aintrp*alu(indlu+1)+(1-aintrp)*alu(indlu)
    ! change nouse lookuptable
    tlog = log(EE/ALIQ)
    !! Bolton(1980) Dew point temperature[K]
    TDPT=(CLIQ-DLIQ*TLOG)/(BLIQ-TLOG)
    !! Bolton(1980) 
    TSAT=TDPT - (0.212_RP+1.571E-3_RP*(TDPT-TEM00)-4.36E-4_RP*(T1-TEM00))*(T1-TDPT)
    !      TSAT = 2840._RP/(3.5_RP - log(ee) -4.805_RP) +55
    !! Bolton(1980) emanuel 132p (4.7.9)
    THT=T1*(P00/P1)**(0.2854_RP*(1._RP-0.28_RP*Q1))
    THT1=THT*EXP((C1/TSAT-C2)*Q1*(1._RP+0.81_RP*Q1))
    !
    return
  end subroutine ENVIRTHT
  !***********************************************************************
  subroutine kf_lutab!(SVP1,SVP2,SVP3,SVPT0)
    !
    !  This subroutine is a lookup table.
    !  Given a series of series of saturation equivalent potential 
    !  temperatures, the temperature is calculated.
    !
    !--------------------------------------------------------------------
    use scale_const, only :&
         CP => CONST_CPdry   , &
         PRE00 => CONST_PRE00, &
         GRAV  => CONST_GRAV
    IMPLICIT NONE
    !--------------------------------------------------------------------
    ! Lookup table variables
    !     INTEGER, SAVE, PARAMETER :: KFNT=250,KFNP=220
    !     REAL, SAVE, DIMENSION(1:KFNT,1:KFNP) :: TTAB,QSTAB
    !     REAL, SAVE, DIMENSION(1:KFNP) :: THE0K
    !     REAL, SAVE, DIMENSION(1:200) :: ALU
    !     REAL, SAVE :: RDPR,RDTHK,PLUTOP
    ! End of Lookup table variables
    !!    use scale_const, only : &
    !!         SVPT0 -> CONST_TEM00,&
    integer :: KP,IT,ITCNT,I
    real(RP) :: DTH=1._RP,TMIN=150._RP,TOLER=0.001_RP
    real(RP) :: PBOT,DPR,                               &
         TEMP,P,ES,QS,PI,THES,TGUES,THGUES,F0,T1,T0,THGS,F1,DT, &
         ASTRT,AINC,A1,THTGS
    !    REAL    :: ALIQ,BLIQ,CLIQ,DLIQ,SVP1,SVP2,SVP3,SVPT0
    !! to module variables      
    !    real(RP)    :: ALIQ,BLIQ,CLIQ,DLIQ
    !!    real(RP), INTENT(IN)    :: SVP1,SVP2,SVP3,SVPT0
    !
    ! equivalent potential temperature increment
    !    data dth/1._RP/
    ! minimum starting temp 
    !    data tmin/150._RP/
    ! tolerance for accuracy of temperature 
    !    data toler/0.001_RP/
    ! top pressure (pascals)
    plutop=5000.0_RP
    ! bottom pressure (pascals)
    pbot=110000.0_RP
    !! to module variable
    !!    ALIQ = SVP1*1000.
    !!    BLIQ = SVP2
    !!    CLIQ = SVP2*SVPT0
    !!    DLIQ = SVP3

    !
    ! compute parameters
    !
    ! 1._over_(sat. equiv. theta increment)
    rdthk=1._RP/dth
    ! pressure increment
    !
    DPR=(PBOT-PLUTOP)/REAL(KFNP-1)
    !      dpr=(pbot-plutop)/REAL(kfnp-1)
    ! 1._over_(pressure increment)
    rdpr=1._RP/dpr
    ! compute the spread of thes
    !     thespd=dth*(kfnt-1)
    !
    ! calculate the starting sat. equiv. theta
    !
    temp=tmin 
    p=plutop-dpr
    do kp=1,kfnp
       p=p+dpr
       es=aliq*exp((bliq*temp-cliq)/(temp-dliq))
       qs=0.622_RP*es/(p-es)
       pi=(1.e5_RP/p)**(0.2854_RP*(1.-0.28_RP*qs))
       the0k(kp)=temp*pi*exp((3374.6525_RP/temp-2.5403_RP)*qs*        &
            (1._RP+0.81_RP*qs))
    enddo
    !
    ! compute temperatures for each sat. equiv. potential temp.
    !
    p=plutop-dpr
    do kp=1,kfnp
       thes=the0k(kp)-dth
       p=p+dpr
       do it=1,kfnt
          ! define sat. equiv. pot. temp.
          thes=thes+dth
          ! iterate to find temperature
          ! find initial guess
          if(it.eq.1) then
             tgues=tmin
          else
             tgues=ttab(it-1,kp)
          endif
          es=aliq*exp((bliq*tgues-cliq)/(tgues-dliq))
          qs=0.622_RP*es/(p-es)
          pi=(1.e5_RP/p)**(0.2854_RP*(1._RP-0.28_RP*qs))
          thgues=tgues*pi*exp((3374.6525_RP/tgues-2.5403_RP)*qs*      &
               (1._RP + 0.81_RP*qs))
          f0=thgues-thes
          t1=tgues-0.5_RP*f0
          t0=tgues
          itcnt=0
          ! iteration loop
          do itcnt=1,11
             es=aliq*exp((bliq*t1-cliq)/(t1-dliq))
             qs=0.622_RP*es/(p-es)
             pi=(1.e5_RP/p)**(0.2854_RP*(1._RP-0.28_RP*qs))
             thtgs=t1*pi*exp((3374.6525_RP/t1-2.5403_RP)*qs*(1._RP + 0.81_RP*qs))
             f1=thtgs-thes
             if(abs(f1).lt.toler)then
                exit
             endif
             !           itcnt=itcnt+1
             dt=f1*(t1-t0)/(f1-f0)
             t0=t1
             f0=f1
             t1=t1-dt
          enddo
          ttab(it,kp)=t1
          qstab(it,kp)=qs
       enddo
    enddo
    !
    ! lookup table for tlog(emix/aliq)
    !
    ! set up intial variable for lookup tables
    !
    astrt=1.e-3_RP
    ainc=0.075_RP
    !
    a1=astrt-ainc
    do i=1,200
       a1=a1+ainc
       alu(i)=log(a1)
    enddo
    !GdCP is g/cp add for SCALE
    GdCP = - GRAV/CP ! inital set
    return
  end subroutine kf_lutab

end module scale_atmos_phy_cp_kf
