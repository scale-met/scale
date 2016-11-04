!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics
!!
!! @par Description
!!          Cloud Microphysics by Kessler-type parametarization
!!          Reference: Kessler(1969)
!!                     Klemp and Wilhelmson(1978)
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2012-01-14 (Y.Miyamoto) [new]
!! @li      2012-03-23 (H.Yashiro)  [mod] Explicit index parameter inclusion
!! @li      2012-12-23 (H.Yashiro)  [mod] Reconstruction
!! @li      2015-09-08 (Y.Sato)     [add] Add evaporated cloud number concentration
!!
!<
!-------------------------------------------------------------------------------
#include "inc_openmp.h"
module scale_atmos_phy_mp_kessler
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index

  use scale_atmos_hydrometer, only: &
       N_HYD, &
       I_QV, &
       I_QC, &
       I_QR, &
       I_HC, &
       I_HR
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_MP_kessler_config
  public :: ATMOS_PHY_MP_kessler_setup
  public :: ATMOS_PHY_MP_kessler
  public :: ATMOS_PHY_MP_kessler_CloudFraction
  public :: ATMOS_PHY_MP_kessler_EffectiveRadius
  public :: ATMOS_PHY_MP_kessler_Mixingratio

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, parameter :: QA_MP  = 3

  character(len=H_SHORT), public, target :: ATMOS_PHY_MP_kessler_NAME(QA_MP)
  character(len=H_MID)  , public, target :: ATMOS_PHY_MP_kessler_DESC(QA_MP)
  character(len=H_SHORT), public, target :: ATMOS_PHY_MP_kessler_UNIT(QA_MP)

  real(RP), public, target :: ATMOS_PHY_MP_kessler_DENS(N_HYD) ! hydrometeor density [kg/m3]=[g/L]

  data ATMOS_PHY_MP_kessler_NAME / &
                 'QV', &
                 'QC', &
                 'QR'  /

  data ATMOS_PHY_MP_kessler_DESC / &
                 'Ratio of Water Vapor mass to total mass (Specific humidity)',   &
                 'Ratio of Cloud Water mass to total mass', &
                 'Ratio of Rain Water mass to total mass'   /

  data ATMOS_PHY_MP_kessler_UNIT / &
                 'kg/kg', &
                 'kg/kg', &
                 'kg/kg'  /

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: MP_kessler
  private :: MP_kessler_vterm

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,  private, parameter   :: I_mp_QC = 1
  integer,  private, parameter   :: I_mp_QR = 2
  integer,  private              :: QS_MP
  integer,  private              :: QE_MP

  logical,  private              :: MP_donegative_fixer = .true.  ! apply negative fixer?
  logical,  private              :: MP_doprecipitation  = .true.  ! apply sedimentation (precipitation)?
  logical,  private              :: MP_couple_aerosol   = .false. ! Consider CCN effect ?
  real(RP), private              :: MP_limit_negative   = 1.0_RP  ! Abort if abs(fixed negative vaue) > abs(MP_limit_negative)

  real(RP), private, allocatable :: factor_vterm(:) ! collection factor for terminal velocity of QR

  integer,  private              :: MP_ntmax_sedimentation = 1 ! number of time step for sedimentation
  integer,  private              :: MP_NSTEP_SEDIMENTATION
  real(RP), private              :: MP_RNSTEP_SEDIMENTATION
  real(DP), private              :: MP_DTSEC_SEDIMENTATION

  logical,  private              :: first = .true.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_MP_kessler_config( &
       MP_TYPE, &
       QA, QS   )
    use scale_process, only: &
       PRC_MPIstop
    use scale_tracer, only: &
       TRACER_regist
    use scale_atmos_hydrometer, only: &
       ATMOS_HYDROMETER_regist
    implicit none
    character(len=*), intent(in) :: MP_TYPE
    integer, intent(out) :: QA
    integer, intent(out) :: QS
    !---------------------------------------------------------------------------

    if ( MP_TYPE /= 'KESSLER' ) then
       write(*,*) 'xxx ATMOS_PHY_MP_TYPE is not KESSLER. Check!'
       call PRC_MPIstop
    endif

    call ATMOS_HYDROMETER_regist( QS_MP, & ! (out)
                                  1, 2, 0,  & ! (in)
                                  ATMOS_PHY_MP_kessler_NAME, & ! (in)
                                  ATMOS_PHY_MP_kessler_DESC, & ! (in)
                                  ATMOS_PHY_MP_kessler_UNIT  ) ! (in)
    QA = QA_MP
    QS = QS_MP
    QE_MP = QS_MP + QA_MP - 1

    I_QV = QS
    I_QC = QS + I_mp_QC
    I_QR = QS + I_mp_QR

    return
  end subroutine ATMOS_PHY_MP_kessler_config

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_MP_kessler_setup
    use scale_process, only: &
       PRC_MPIstop
    use scale_const, only: &
       CONST_UNDEF, &
       CONST_DWATR
    use scale_time, only: &
       TIME_DTSEC_ATMOS_PHY_MP
    use scale_grid, only: &
       CDZ => GRID_CDZ
    implicit none

    real(RP), parameter :: max_term_vel = 10.0_RP  !-- terminal velocity for calculate dt of sedimentation
    integer :: nstep_max

    NAMELIST / PARAM_ATMOS_PHY_MP / &
       MP_doprecipitation,     &
       MP_donegative_fixer,    &
       MP_limit_negative,      &
       MP_ntmax_sedimentation, &
       MP_couple_aerosol

    integer :: ierr

    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[Cloud Microphysics] / Categ[ATMOS PHYSICS] / Origin[SCALElib]'
    if( IO_L ) write(IO_FID_LOG,*) '*** KESSLER-type 1-moment bulk 3 category'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_MP,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_ATMOS_PHY_MP. Check!'
       call PRC_MPIstop
    endif
    if( IO_LNML ) write(IO_FID_LOG,nml=PARAM_ATMOS_PHY_MP)

    if( MP_couple_aerosol ) then
       write(*,*) 'xxx MP_aerosol_couple should be .false. for KESSLER type MP!'
       call PRC_MPIstop
    endif

    if ( IO_L ) write(IO_FID_LOG,*)
    if ( IO_L ) write(IO_FID_LOG,*) '*** Enable negative fixer?                    : ', MP_donegative_fixer
    if ( IO_L ) write(IO_FID_LOG,*) '*** Value limit of negative fixer (abs)       : ', abs(MP_limit_negative)
    if ( IO_L ) write(IO_FID_LOG,*) '*** Enable sedimentation (precipitation)?     : ', MP_doprecipitation

    nstep_max = int( ( TIME_DTSEC_ATMOS_PHY_MP * max_term_vel ) / minval( CDZ ) )
    MP_ntmax_sedimentation = max( MP_ntmax_sedimentation, nstep_max )

    MP_NSTEP_SEDIMENTATION  = MP_ntmax_sedimentation
    MP_RNSTEP_SEDIMENTATION = 1.0_RP / real(MP_ntmax_sedimentation,kind=RP)
    MP_DTSEC_SEDIMENTATION  = TIME_DTSEC_ATMOS_PHY_MP * MP_RNSTEP_SEDIMENTATION

    if ( IO_L ) write(IO_FID_LOG,*) '*** Timestep of sedimentation is divided into : ', MP_ntmax_sedimentation, 'step'
    if ( IO_L ) write(IO_FID_LOG,*) '*** DT of sedimentation                       : ', MP_DTSEC_SEDIMENTATION, '[s]'

    ATMOS_PHY_MP_kessler_DENS(:)    = CONST_UNDEF
    ATMOS_PHY_MP_kessler_DENS(I_HC) = CONST_DWATR
    ATMOS_PHY_MP_kessler_DENS(I_HR) = CONST_DWATR

    allocate( factor_vterm(KA) )

    return
  end subroutine ATMOS_PHY_MP_kessler_setup

  !-----------------------------------------------------------------------------
  !> Cloud Microphysics
  subroutine ATMOS_PHY_MP_kessler( &
       DENS,      &
       MOMZ,      &
       MOMX,      &
       MOMY,      &
       RHOT,      &
       QTRC,      &
       CCN,       &
       EVAPORATE, &
       SFLX_rain, &
       SFLX_snow  )
    use scale_grid_index
    use scale_comm, only: &
       COMM_horizontal_mean
    use scale_history, only: &
       HIST_in
    use scale_tracer, only: &
       QA, &
       TRACER_R, &
       TRACER_CV, &
       TRACER_MASS
    use scale_atmos_thermodyn, only: &
       THERMODYN_rhoe        => ATMOS_THERMODYN_rhoe,       &
       THERMODYN_rhot        => ATMOS_THERMODYN_rhot,       &
       THERMODYN_temp_pres_E => ATMOS_THERMODYN_temp_pres_E
    use scale_atmos_phy_mp_common, only: &
       MP_negative_fixer        => ATMOS_PHY_MP_negative_fixer,       &
       MP_precipitation         => ATMOS_PHY_MP_precipitation,        &
       MP_saturation_adjustment => ATMOS_PHY_MP_saturation_adjustment
    implicit none

    real(RP), intent(inout) :: DENS(KA,IA,JA)
    real(RP), intent(inout) :: MOMZ(KA,IA,JA)
    real(RP), intent(inout) :: MOMX(KA,IA,JA)
    real(RP), intent(inout) :: MOMY(KA,IA,JA)
    real(RP), intent(inout) :: RHOT(KA,IA,JA)
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA)
    real(RP), intent(in)    :: CCN(KA,IA,JA)
    real(RP), intent(out)   :: EVAPORATE(KA,IA,JA)   !--- number of evaporated cloud [/m3]
    real(RP), intent(out)   :: SFLX_rain(IA,JA)
    real(RP), intent(out)   :: SFLX_snow(IA,JA)

    real(RP) :: RHOE_t(KA,IA,JA)
    real(RP) :: QTRC_t(KA,IA,JA,QA)
    real(RP) :: RHOE  (KA,IA,JA)
    real(RP) :: temp  (KA,IA,JA)
    real(RP) :: pres  (KA,IA,JA)

    real(RP) :: vterm   (KA,IA,JA,QA_MP-1) ! terminal velocity of each tracer [m/s]
    real(RP) :: FLX_rain(KA,IA,JA)
    real(RP) :: FLX_snow(KA,IA,JA)
    real(RP) :: wflux_rain(KA,IA,JA)
    real(RP) :: wflux_snow(KA,IA,JA)
    integer  :: step

    real(RP) :: rho_prof(KA) ! averaged profile of rho
    integer  :: k, i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*) '*** Physics step: Cloud microphysics(kessler)'

    if ( first ) then
       ! Calculate collection factor for terminal velocity of QR
       call COMM_horizontal_mean( rho_prof(:), DENS(:,:,:) )
       rho_prof(:) = rho_prof(:) * 1.E-3_RP ! [kg/m3]->[g/cc]

       do k = KS, KE
          factor_vterm(k) = sqrt( rho_prof(KS)/rho_prof(k) )
       enddo
       first = .false.
    endif

    EVAPORATE(:,:,:)   = 0.0_RP

    if ( MP_donegative_fixer ) then
       call MP_negative_fixer( DENS(:,:,:),      & ! [INOUT]
                               RHOT(:,:,:),      & ! [INOUT]
                               QTRC(:,:,:,:),    & ! [INOUT]
                               I_QV,             & ! [IN]
                               MP_limit_negative ) ! [IN]
    endif

    call THERMODYN_rhoe( RHOE(:,:,:),   & ! [OUT]
                         RHOT(:,:,:),   & ! [IN]
                         QTRC(:,:,:,:), & ! [IN]
                         TRACER_CV(:),  & ! [IN]
                         TRACER_R(:),   & ! [IN]
                         TRACER_MASS(:) ) ! [IN]

    !##### MP Main #####
    RHOE_t(:,:,:)   = 0.0_RP
    QTRC_t(:,:,:,:) = 0.0_RP

    call MP_kessler( RHOE_t(:,:,:),   & ! [INOUT]
                     QTRC_t(:,:,:,:), & ! [INOUT]
                     RHOE  (:,:,:),   & ! [INOUT]
                     QTRC  (:,:,:,:), & ! [INOUT]
                     DENS  (:,:,:)    ) ! [IN]

    if ( MP_doprecipitation ) then

       FLX_rain(:,:,:) = 0.0_RP
       FLX_snow(:,:,:) = 0.0_RP

       vterm(:,:,:,:) = 0.0_RP

       do step = 1, MP_NSTEP_SEDIMENTATION

         call MP_kessler_vterm( vterm(:,:,:,:), & ! [OUT]
                                DENS (:,:,:),   & ! [IN]
                                QTRC (:,:,:,:)  ) ! [IN]

         call THERMODYN_temp_pres_E( temp(:,:,:),   & ! [OUT]
                                     pres(:,:,:),   & ! [OUT]
                                     DENS(:,:,:),   & ! [IN]
                                     RHOE(:,:,:),   & ! [IN]
                                     QTRC(:,:,:,:), & ! [IN]
                                     TRACER_CV(:),  & ! [IN]
                                     TRACER_R(:),   & ! [IN]
                                     TRACER_MASS(:) ) ! [IN]

         call MP_precipitation( wflux_rain(:,:,:),     & ! [OUT]
                                wflux_snow(:,:,:),     & ! [OUT]
                                DENS    (:,:,:),       & ! [INOUT]
                                MOMZ    (:,:,:),       & ! [INOUT]
                                MOMX    (:,:,:),       & ! [INOUT]
                                MOMY    (:,:,:),       & ! [INOUT]
                                RHOE    (:,:,:),       & ! [INOUT]
                                QTRC    (:,:,:,:),     & ! [INOUT]
                                QA_MP,                 & ! [IN]
                                QS_MP,                 & ! [IN]
                                vterm   (:,:,:,:),     & ! [IN]
                                temp    (:,:,:),       & ! [IN]
                                TRACER_CV(:),          & ! [IN]
                                MP_DTSEC_SEDIMENTATION ) ! [IN]

         do j = JS, JE
         do i = IS, IE
         do k = KS-1, KE
            FLX_rain(k,i,j) = FLX_rain(k,i,j) + wflux_rain(k,i,j) * MP_RNSTEP_SEDIMENTATION
            FLX_snow(k,i,j) = FLX_snow(k,i,j) + wflux_snow(k,i,j) * MP_RNSTEP_SEDIMENTATION
         enddo
         enddo
         enddo

       enddo

    else
       vterm(:,:,:,:) = 0.0_RP

       FLX_rain(:,:,:) = 0.0_RP
       FLX_snow(:,:,:) = 0.0_RP
    endif

    call MP_saturation_adjustment( RHOE_t(:,:,:),       & ! [INOUT]
                                   QTRC_t(:,:,:,:),     & ! [INOUT]
                                   RHOE  (:,:,:),       & ! [INOUT]
                                   QTRC  (:,:,:,:),     & ! [INOUT]
                                   DENS  (:,:,:),       & ! [IN]
                                   I_QV,                & ! [IN]
                                   I_QC,                & ! [IN]
                                   -1,                  & ! [IN]
                                   flag_liquid = .true. ) ! [IN]

    call HIST_in( vterm(:,:,:,I_mp_QR), 'Vterm_QR', 'terminal velocity of QR', 'm/s' )

    !##### END MP Main #####

    call THERMODYN_rhot( RHOT(:,:,:),   & ! [OUT]
                         RHOE(:,:,:),   & ! [IN]
                         QTRC(:,:,:,:), & ! [IN]
                         TRACER_CV(:),  & ! [IN]
                         TRACER_R(:),   & ! [IN]
                         TRACER_MASS(:) ) ! [IN]

    if ( MP_donegative_fixer ) then
       call MP_negative_fixer( DENS(:,:,:),      & ! [INOUT]
                               RHOT(:,:,:),      & ! [INOUT]
                               QTRC(:,:,:,:),    & ! [INOUT]
                               I_QV,             & ! [IN]
                               MP_limit_negative ) ! [IN]
    endif

    SFLX_rain(:,:) = FLX_rain(KS-1,:,:)
    SFLX_snow(:,:) = FLX_snow(KS-1,:,:)

    return
  end subroutine ATMOS_PHY_MP_kessler

  !-----------------------------------------------------------------------------
  !> Kessler-type warm rain microphysics
  subroutine MP_kessler( &
       RHOE_t, &
       QTRC_t, &
       RHOE0,  &
       QTRC0,  &
       DENS0   )
    use scale_tracer, only: &
       QA, &
       TRACER_R, &
       TRACER_CV, &
       TRACER_MASS
    use scale_time, only: &
       dt => TIME_DTSEC_ATMOS_PHY_MP
    use scale_atmos_thermodyn, only: &
       THERMODYN_temp_pres_E => ATMOS_THERMODYN_temp_pres_E
    use scale_atmos_hydrometer, only: &
       LHV
    use scale_atmos_saturation, only: &
       SATURATION_dens2qsat_liq => ATMOS_SATURATION_dens2qsat_liq
    implicit none

    real(RP), intent(inout) :: RHOE_t(KA,IA,JA)    ! tendency rhoe             [J/m3/s]
    real(RP), intent(inout) :: QTRC_t(KA,IA,JA,QA) ! tendency tracer           [kg/kg/s]
    real(RP), intent(inout) :: RHOE0 (KA,IA,JA)    ! density * internal energy [J/m3]
    real(RP), intent(inout) :: QTRC0 (KA,IA,JA,QA) ! mass concentration        [kg/kg]
    real(RP), intent(in)    :: DENS0 (KA,IA,JA)    ! density                   [kg/m3]

    ! working
    real(RP) :: QSAT (KA,IA,JA) ! saturated water vapor [kg/kg]
    real(RP) :: TEMP0(KA,IA,JA) ! temperature           [K]
    real(RP) :: PRES0(KA,IA,JA) ! pressure              [Pa]

    real(RP) :: dens (KA)
    real(RP) :: rhoe (KA)
    real(RP) :: temp (KA)
    real(RP) :: pres (KA)
    real(RP) :: qv   (KA)
    real(RP) :: qc   (KA)
    real(RP) :: qr   (KA)
    real(RP) :: qsatl(KA)
    real(RP) :: Sliq (KA)

    ! tendency
    real(RP) :: dq_evap(KA) ! tendency q (evaporation)
    real(RP) :: dq_auto(KA) ! tendency q (autoconversion)
    real(RP) :: dq_accr(KA) ! tendency q (accretion)

    real(RP) :: dqv, dqc, dqr
    real(RP) :: vent_factor
    real(RP) :: rdt

    integer :: k, i, j
    !---------------------------------------------------------------------------

    call PROF_rapstart('MP_kessler', 3)

    rdt = 1.0_RP / dt

    call THERMODYN_temp_pres_E( TEMP0(:,:,:),   & ! [OUT]
                                PRES0(:,:,:),   & ! [OUT]
                                DENS0(:,:,:),   & ! [IN]
                                RHOE0(:,:,:),   & ! [IN]
                                QTRC0(:,:,:,:), & ! [IN]
                                TRACER_CV(:),   & ! [IN]
                                TRACER_R(:),    & ! [IN]
                                TRACER_MASS(:)  ) ! [IN]

    call SATURATION_dens2qsat_liq( QSAT (:,:,:), & ! [OUT]
                                   TEMP0(:,:,:), & ! [IN]
                                   DENS0(:,:,:)  ) ! [IN]

    !$omp parallel do private(i,j,k,dens,rhoe,temp,pres,qv,qc,qr,qsatl,Sliq,dq_evap,dq_auto,dq_accr,vent_factor,dqc,dqr,dqv) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE

       ! store to work
       dens (KS:KE) = DENS0(KS:KE,i,j)
       rhoe (KS:KE) = RHOE0(KS:KE,i,j)
       temp (KS:KE) = TEMP0(KS:KE,i,j)
       pres (KS:KE) = PRES0(KS:KE,i,j)
       qv   (KS:KE) = QTRC0(KS:KE,i,j,I_QV)
       qc   (KS:KE) = QTRC0(KS:KE,i,j,I_QC)
       qr   (KS:KE) = QTRC0(KS:KE,i,j,I_QR)
       qsatl(KS:KE) = QSAT (KS:KE,i,j)
       Sliq (KS:KE) = qv(KS:KE) / qsatl(KS:KE)

       ! Auto-conversion (QC->QR)
       do k = KS, KE
          dq_auto(k) = 1.E-3_RP * max( qc(k)-1.E-3_RP, 0.0_RP )
       enddo

       ! Accretion (QC->QR)
       do k = KS, KE
          dq_accr(k) = 2.2_RP * qc(k) * qr(k)**0.875_RP
       enddo

       ! Evaporation (QR->QV)
       do k = KS, KE
          vent_factor = 1.6_RP + 124.9_RP * ( dens(k)*qr(k) )**0.2046_RP

          dq_evap(k) = ( 1.0_RP-min(Sliq(k),1.0_RP) ) / dens(k) * vent_factor  &
                     * ( dens(k)*qr(k) )**0.525_RP / ( 5.4E5_RP + 2.55E8_RP / ( pres(k)*qsatl(k) ) )
       enddo

       ! limiter
       do k = KS, KE
          dqc = (-dq_auto(k)-dq_accr(k)            )

          QTRC_t(k,i,j,I_QC) = QTRC_t(k,i,j,I_QC) + max( dqc, -qc(k)*rdt )
       enddo

       do k = KS, KE
          dqr = ( dq_auto(k)+dq_accr(k)-dq_evap(k) )

          QTRC_t(k,i,j,I_QR) = QTRC_t(k,i,j,I_QR) + max( dqr, -qr(k)*rdt )
       enddo

       do k = KS, KE
          dqv = - ( QTRC_t(k,i,j,I_QC) &
                  + QTRC_t(k,i,j,I_QR) )

          QTRC_t(k,i,j,I_QV) = QTRC_t(k,i,j,I_QV) + max( dqv, -qv(k)*rdt )
       enddo

       do k = KS, KE
          RHOE_t(k,i,j) = RHOE_t(k,i,j) - dens(k) * ( LHV * QTRC_t(k,i,j,I_QV) )
       enddo
    enddo
    enddo

    ! mass & energy update
    !$omp parallel do private(i,j, k) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
       QTRC0(KS:KE,i,j,I_QV) = QTRC0(KS:KE,i,j,I_QV) + QTRC_t(KS:KE,i,j,I_QV) * dt
       QTRC0(KS:KE,i,j,I_QC) = QTRC0(KS:KE,i,j,I_QC) + QTRC_t(KS:KE,i,j,I_QC) * dt
       QTRC0(KS:KE,i,j,I_QR) = QTRC0(KS:KE,i,j,I_QR) + QTRC_t(KS:KE,i,j,I_QR) * dt

       RHOE0(KS:KE,i,j) = RHOE0(KS:KE,i,j) + RHOE_t(KS:KE,i,j) * dt
    enddo
    enddo

    call PROF_rapend  ('MP_kessler', 3)

    return
  end subroutine MP_kessler

  !-----------------------------------------------------------------------------
  !> Kessler-type warm rain microphysics (terminal velocity)
  subroutine MP_kessler_vterm( &
       vterm, &
       DENS0, &
       QTRC0  )
    use scale_tracer, only: &
       QA
    use scale_atmos_refstate, only: &
       REFSTATE_dens => ATMOS_REFSTATE_dens
    implicit none

    real(RP), intent(out) :: vterm(KA,IA,JA,QA_MP-1)
    real(RP), intent(in)  :: DENS0(KA,IA,JA)
    real(RP), intent(in)  :: QTRC0(KA,IA,JA,QA)

    real(RP) :: zerosw

    integer :: k, i, j
    !---------------------------------------------------------------------------

    !$omp parallel do private(i,j,k,zerosw) OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       vterm(k,i,j,I_mp_QC) = 0.0_RP
    enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       zerosw = 0.5_RP - sign(0.5_RP, QTRC0(k,i,j,I_QR) - 1.E-12_RP )

       vterm(k,i,j,I_mp_QR) = - 36.34_RP * ( DENS0(k,i,j) * ( QTRC0(k,i,j,I_QR) + zerosw ) )**0.1364_RP &
                            * REFSTATE_dens(KS,i,j) / REFSTATE_dens(k,i,j) * ( 1.0_RP - zerosw )
    enddo
    enddo
    enddo

    return
  end subroutine MP_kessler_vterm

  !-----------------------------------------------------------------------------
  !> Calculate Cloud Fraction
  subroutine ATMOS_PHY_MP_kessler_CloudFraction( &
       cldfrac, &
       QTRC     )
    use scale_grid_index
    use scale_tracer, only: &
       QA
    use scale_const, only: &
       EPS => CONST_EPS
    implicit none

    real(RP), intent(out) :: cldfrac(KA,IA,JA)
    real(RP), intent(in)  :: QTRC   (KA,IA,JA,QA)

    real(RP) :: qhydro
    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    do j  = JS, JE
    do i  = IS, IE
    do k  = KS, KE
       qhydro = 0.0_RP
       do iq = QS_MP+1, QE_MP
          qhydro = qhydro + QTRC(k,i,j,iq)
       enddo
       cldfrac(k,i,j) = 0.5_RP + sign(0.5_RP,qhydro-EPS)
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_MP_kessler_CloudFraction

  !-----------------------------------------------------------------------------
  !> Calculate Effective Radius
  subroutine ATMOS_PHY_MP_kessler_EffectiveRadius( &
       Re,    &
       QTRC0, &
       DENS0, &
       TEMP0  )
    use scale_grid_index
    use scale_tracer, only: &
       QA
    use scale_atmos_hydrometer, only: &
       N_HYD
    implicit none
    real(RP), intent(out) :: Re   (KA,IA,JA,N_HYD) ! effective radius          [cm]
    real(RP), intent(in)  :: QTRC0(KA,IA,JA,QA)    ! tracer mass concentration [kg/kg]
    real(RP), intent(in)  :: DENS0(KA,IA,JA)       ! density                   [kg/m3]
    real(RP), intent(in)  :: TEMP0(KA,IA,JA)       ! temperature               [K]

    real(RP), parameter :: um2cm = 100.0_RP

    integer :: ih
    !---------------------------------------------------------------------------

    do ih = 1, N_HYD
       if ( ih == I_HC ) then
          Re(:,:,:,ih) =   8.E-6_RP * um2cm
       else if ( ih == I_HR ) then
          Re(:,:,:,ih) = 100.E-6_RP * um2cm
       else
          Re(:,:,:,ih) = 0.0_RP
       end if
    end do

    return
  end subroutine ATMOS_PHY_MP_kessler_EffectiveRadius

  !-----------------------------------------------------------------------------
  !> Calculate mixing ratio of each category
  subroutine ATMOS_PHY_MP_kessler_Mixingratio( &
       Qe,   &
       QTRC0 )
    use scale_grid_index
    use scale_tracer, only: &
       QA
    use scale_atmos_hydrometer, only: &
       N_HYD
    implicit none
    real(RP), intent(out) :: Qe   (KA,IA,JA,N_HYD) ! mixing ratio of each cateory [kg/kg]
    real(RP), intent(in)  :: QTRC0(KA,IA,JA,QA)    ! tracer mass concentration [kg/kg]

    integer :: ih
    !---------------------------------------------------------------------------

    do ih = 1, N_HYD
       if ( ih == I_HC ) then
          Qe(:,:,:,ih) = QTRC0(:,:,:,I_QC)
       else if ( ih == I_HR ) then
          Qe(:,:,:,ih) = QTRC0(:,:,:,I_QR)
       else
          Qe(:,:,:,ih) = 0.0_RP
       end if
    end do

    return
  end subroutine ATMOS_PHY_MP_kessler_Mixingratio

end module scale_atmos_phy_mp_kessler
