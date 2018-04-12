!-------------------------------------------------------------------------------
!> module atmosphere / physics / microphysics / Kessler
!!
!! @par Description
!!          Cloud Microphysics by Kessler-type parametarization
!!          Reference: Kessler(1969)
!!                     Klemp and Wilhelmson(1978)
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_phy_mp_kessler
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_MP_kessler_setup
  public :: ATMOS_PHY_MP_kessler_adjustment
  public :: ATMOS_PHY_MP_kessler_terminal_velocity
  public :: ATMOS_PHY_MP_kessler_effective_radius
  public :: ATMOS_PHY_MP_kessler_cloud_fraction
  public :: ATMOS_PHY_MP_kessler_qtrc2qhyd
  public :: ATMOS_PHY_MP_kessler_qhyd2qtrc

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, private, parameter :: QA_MP  = 3

  integer,                parameter, public :: ATMOS_PHY_MP_kessler_ntracers = QA_MP
  integer,                parameter, public :: ATMOS_PHY_MP_kessler_nwaters = 2
  integer,                parameter, public :: ATMOS_PHY_MP_kessler_nices = 0
  character(len=H_SHORT), parameter, public :: ATMOS_PHY_MP_kessler_tracer_names(QA_MP) = (/ &
       'QV', &
       'QC', &
       'QR'  /)
  character(len=H_MID)  , parameter, public :: ATMOS_PHY_MP_kessler_tracer_descriptions(QA_MP) = (/ &
       'Ratio of Water Vapor mass to total mass (Specific humidity)', &
       'Ratio of Cloud Water mass to total mass                    ', &
       'Ratio of Rain Water mass to total mass                     '/)
  character(len=H_SHORT), parameter, public :: ATMOS_PHY_MP_kessler_tracer_units(QA_MP) = (/ &
       'kg/kg',  &
       'kg/kg',  &
       'kg/kg'   /)


  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: MP_kessler

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,  private, parameter   :: I_QV = 1
  integer,  private, parameter   :: I_QC = 2
  integer,  private, parameter   :: I_QR = 3

  integer,  private, parameter   :: I_hyd_QC =  1
  integer,  private, parameter   :: I_hyd_QR =  2

  logical,  private              :: flag_liquid = .true.     ! warm rain
  logical,  private              :: couple_aerosol = .false. ! Consider CCN effect ?

  real(RP), private, parameter   :: re_qc =  8.E-6_RP        ! effective radius for cloud water

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> ATMOS_PHY_MP_kessler_setup
  !! Setup
  !<
  subroutine ATMOS_PHY_MP_kessler_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_MP_kessler_setup",*) 'Setup'
    LOG_INFO("ATMOS_PHY_MP_kessler_setup",*) 'KESSLER-type 1-moment bulk 3 category'

    if( couple_aerosol ) then
       LOG_ERROR("ATMOS_PHY_MP_kessler_setup",*) 'MP_aerosol_couple should be .false. for KESSLER type MP!'
       call PRC_abort
    endif

    return
  end subroutine ATMOS_PHY_MP_kessler_setup

  !-----------------------------------------------------------------------------
  !> ATMOS_PHY_MP_kessler_adjustment
  !! calculate state after saturation process
  !<
  subroutine ATMOS_PHY_MP_kessler_adjustment( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       DENS, PRES,               &
       dt,                       &
       TEMP, QTRC, CPtot, CVtot, &
       RHOE_t, EVAPORATE         )
    use scale_const, only: &
       DWATR => CONST_DWATR, &
       PI    => CONST_PI
    use scale_file_history, only: &
       FILE_HISTORY_in
    use scale_atmos_phy_mp_common, only: &
       MP_saturation_adjustment => ATMOS_PHY_MP_saturation_adjustment
    implicit none
    integer,  intent(in) :: KA, KS, KE
    integer,  intent(in) :: IA, IS, IE
    integer,  intent(in) :: JA, JS, JE

    real(RP), intent(in) :: DENS (KA,IA,JA)         ! density [kg/m3]
    real(RP), intent(in) :: PRES (KA,IA,JA)         ! pressure [Pa]
    real(DP), intent(in) :: dt                      ! time interval of microphysics [s]

    real(RP), intent(inout) :: TEMP(KA,IA,JA)       ! temperature [K]
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA_MP) ! tracer mass concentration [kg/kg]
    real(RP), intent(inout) :: CPtot(KA,IA,JA)      ! total specific heat capacity at constant pressure [J/kg/K]
    real(RP), intent(inout) :: CVtot(KA,IA,JA)      ! total specific heat capacity at constant volume [J/kg/K]

    real(RP), intent(out) :: RHOE_t   (KA,IA,JA)    ! tendency of rhoe [J/m3/s]
    real(RP), intent(out) :: EVAPORATE(KA,IA,JA)    ! number of evaporated cloud [#/m3/s]

    real(RP) :: QTRC_dummy(KA,IA,JA)

    real(RP) :: RHOE_d_sat(KA,IA,JA)

    real(RP) :: QC_t_sat (KA,IA,JA)
    real(RP) :: dens_prof(KA)       ! averaged profile of rho

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    LOG_PROGRESS(*) 'atmosphere / physics / microphysics / Kessler'

    !##### MP Main #####
    call MP_kessler( &
         KA, KS, KE, IA, IS, IE, JA, JS, JE, &
         DENS(:,:,:), PRES(:,:,:),   & ! [IN]
         dt,                         & ! [IN]
         TEMP(:,:,:), QTRC(:,:,:,:), & ! [INOUT]
         CPtot(:,:,:), CVtot(:,:,:), & ! [INOUT]
         RHOE_t(:,:,:)               ) ! [OUT]

    ! save value before saturation adjustment
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       QC_t_sat(k,i,j) = QTRC(k,i,j,I_QC)
    enddo
    enddo
    enddo

!OCL XFILL
    QTRC_dummy(:,:,:) = -1.0_RP

    call MP_saturation_adjustment( &
         KA, KS, KE, IA, IS, IE, JA, JS, JE,  &
         DENS(:,:,:),                         & ! [IN]
         flag_liquid,                         & ! [IN]
         TEMP(:,:,:),                         & ! [INOUT]
         QTRC(:,:,:,I_QV),                    & ! [INOUT]
         QTRC(:,:,:,I_QC), QTRC_dummy(:,:,:), & ! [INOUT]
         CPtot(:,:,:), CVtot(:,:,:),          & ! [INOUT]
         RHOE_d_sat(:,:,:)                    ) ! [OUT]

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       RHOE_t(k,i,j) = RHOE_t(k,i,j) + RHOE_d_sat(k,i,j) / dt
    enddo
    enddo
    enddo
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       QC_t_sat(k,i,j) = ( QTRC(k,i,j,I_QC) - QC_t_sat(k,i,j) ) / dt
    enddo
    enddo
    enddo

    call FILE_HISTORY_in( QC_t_sat(:,:,:), 'Pcsat', 'QC production term by satadjust', 'kg/kg/s' )

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       EVAPORATE(k,i,j) = max( -QC_t_sat(k,i,j), 0.0_RP ) & ! if negative, condensation
                        * DENS(k,i,j) / (4.0_RP/3.0_RP*PI*DWATR*re_qc**3) ! mass -> number (assuming constant particle radius as re_qc)
    enddo
    enddo
    enddo

    !##### END MP Main #####

    return
  end subroutine ATMOS_PHY_MP_kessler_adjustment

  !-----------------------------------------------------------------------------
  !> Kessler-type warm rain microphysics
  subroutine MP_kessler( &
         KA, KS, KE, IA, IS, IE, JA, JS, JE, &
         DENS0, PRES0,   &
         dt,             &
         TEMP0, QTRC0,   &
         CPtot0, CVtot0, &
         RHOE_t          )
    use scale_const, only: &
       EPS => CONST_EPS
    use scale_atmos_hydrometeor, only: &
       LHV, &
       CP_VAPOR, &
       CP_WATER, &
       CV_VAPOR, &
       CV_WATER
    use scale_atmos_saturation, only: &
       SATURATION_dens2qsat_liq => ATMOS_SATURATION_dens2qsat_liq
    implicit none
    integer,  intent(in) :: KA, KS, KE
    integer,  intent(in) :: IA, IS, IE
    integer,  intent(in) :: JA, JS, JE

    real(RP), intent(in) :: DENS0(KA,IA,JA)           ! density [kg/m3]
    real(RP), intent(in) :: PRES0(KA,IA,JA)           ! pressure [Pa]
    real(DP), intent(in) :: dt                        ! time interval of microphysics [s]

    real(RP), intent(inout) :: TEMP0 (KA,IA,JA)       ! temperature [K]
    real(RP), intent(inout) :: QTRC0 (KA,IA,JA,QA_MP) ! tracer mass concentration [kg/kg]
    real(RP), intent(inout) :: CPtot0(KA,IA,JA)       ! total specific heat capacity at constant pressure [J/kg/K]
    real(RP), intent(inout) :: CVtot0(KA,IA,JA)       ! total specific heat capacity at constant volume [J/kg/K]

    real(RP), intent(out) :: RHOE_t(KA,IA,JA)         ! tendency of rhoe [J/m3/s]

    real(RP) :: dens
    real(RP) :: temp
    real(RP) :: pres
    real(RP) :: cptot, cvtot
    real(RP) :: qv, qc, qr
    real(RP) :: qv_t, qc_t, qr_t
    real(RP) :: e_t, cp_t, cv_t
    real(RP) :: QSATL ! saturated water vapor for liquid water [kg/kg]
    real(RP) :: Sliq  ! saturated ratio S for liquid water (0-1)

    ! tendency
    real(RP) :: dq_evap ! tendency q (evaporation)
    real(RP) :: dq_auto ! tendency q (autoconversion)
    real(RP) :: dq_accr ! tendency q (accretion)

    real(RP) :: vent_factor
    real(RP) :: rdt

    integer :: k, i, j
    !---------------------------------------------------------------------------

    call PROF_rapstart('MP_kessler', 3)

    rdt = 1.0_RP / dt

    !$omp parallel do default(none) &
    !$omp shared(KS,KE,IS,IE,JS,JE, &
    !$omp        DENS0,TEMP0,PRES0,QTRC0,CPtot0,CVtot0,dt, &
    !$omp        RHOE_t, &
    !$omp        EPS,LHV,CP_VAPOR,CP_WATER,CV_VAPOR,CV_WATER,rdt) &
    !$omp private(k,i,j,dens,temp,pres,cptot,cvtot,qv,qc,qr,qv_t,qc_t,qr_t,e_t,cp_t,cv_t, &
    !$omp         QSATL,Sliq, &
    !$omp         dq_evap,dq_auto,dq_accr,vent_factor) &
    !$omp OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       ! store to work
       dens = DENS0(k,i,j)
       temp = TEMP0(k,i,j)
       pres = PRES0(k,i,j)
       qv   = QTRC0(k,i,j,I_QV)
       qc   = QTRC0(k,i,j,I_QC)
       qr   = QTRC0(k,i,j,I_QR)

       call SATURATION_dens2qsat_liq( &
            temp, dens, & ! [IN]
            QSATL       ) ! [OUT]

       Sliq = qv / max( QSATL, EPS)

       ! Auto-conversion (QC->QR)
       dq_auto = 1.E-3_RP * max( qc-1.E-3_RP, 0.0_RP )

       ! Accretion (QC->QR)
       dq_accr = 2.2_RP * qc * qr**0.875_RP

       ! Evaporation (QR->QV)
       vent_factor = 1.6_RP + 124.9_RP * ( dens*qr )**0.2046_RP
       dq_evap = ( 1.0_RP-min(Sliq,1.0_RP) ) / dens * vent_factor  &
               * ( dens*qr )**0.525_RP / ( 5.4E5_RP + 2.55E8_RP / ( pres*QSATL ) )

       ! limiter
       qc_t = (-dq_auto-dq_accr         )
       qc_t = max( qc_t, -qc*rdt )

       qr_t = ( dq_auto+dq_accr-dq_evap )
       qr_t = max( qr_t, -qr*rdt )

       qv_t = - ( qc_t + qr_t )
       qv_t = max( qv_t, -qv*rdt)

       ! mass & energy update
       QTRC0(k,i,j,I_QV) = QTRC0(k,i,j,I_QV) + qv_t * dt
       QTRC0(k,i,j,I_QC) = QTRC0(k,i,j,I_QC) + qc_t * dt
       QTRC0(k,i,j,I_QR) = QTRC0(k,i,j,I_QR) + qr_t * dt

       e_t = -LHV * qv_t ! internal energy change
       cp_t = CP_VAPOR * qv_t &
            + CP_WATER * ( qc_t + qr_t )
       cptot = CPtot0(k,i,j) + cp_t * dt

       cv_t = CV_VAPOR * qv_t &
            + CV_WATER * ( qc_t + qr_t )
       cvtot = CVtot0(k,i,j) + cv_t * dt

       RHOE_t(k,i,j) = dens * e_t

       TEMP0(k,i,j) = ( temp * CVtot0(k,i,j) + e_t * dt ) / cvtot
       CPtot0(k,i,j) = cptot
       CVtot0(k,i,j) = cvtot

    enddo
    enddo
    enddo

    call PROF_rapend  ('MP_kessler', 3)

    return
  end subroutine MP_kessler

  !-----------------------------------------------------------------------------
  !> Kessler-type warm rain microphysics (terminal velocity)
  subroutine ATMOS_PHY_MP_kessler_terminal_velocity( &
       KA, KS, KE,            &
       DENS0, RHOQ0,          &
       REFSTATE_dens_profile, &
       vterm                  )
    implicit none
    integer,  intent(in)  :: KA, KS, KE

    real(RP), intent(in)  :: DENS0(KA)                  ! density [kg/m3]
    real(RP), intent(in)  :: RHOQ0(KA,QA_MP-1)          ! density of each hydrometeor [kg/m3]
    real(RP), intent(in)  :: REFSTATE_dens_profile(KA)  ! reference state density profile [kg/m3]

    real(RP), intent(out) :: vterm(KA,QA_MP-1) ! terminal velocity of each hydrometeor [m/s]

    real(RP) :: qr
    real(RP) :: zerosw

    integer :: k, i, j
    !---------------------------------------------------------------------------

    do k = KS, KE
       vterm(k,I_hyd_QC) = 0.0_RP
    enddo

    do k = KS, KE
       qr     = RHOQ0(k,I_hyd_QR) / DENS0(k)
       zerosw = 0.5_RP - sign(0.5_RP, qr - 1.E-12_RP )

       vterm(k,I_hyd_QR) = - 36.34_RP * ( DENS0(k) * ( qr + zerosw ) )**0.1364_RP &
                         * REFSTATE_dens_profile(KS) / REFSTATE_dens_profile(k) * ( 1.0_RP - zerosw )
    enddo

    return
  end subroutine ATMOS_PHY_MP_kessler_terminal_velocity

  !-----------------------------------------------------------------------------
  !> Calculate Cloud Fraction
  subroutine ATMOS_PHY_MP_kessler_cloud_fraction( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       QTRC,           &
       mask_criterion, &
       cldfrac         )
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: QTRC(KA,IA,JA,QA_MP-1) ! hydrometeor mass concentration [kg/kg]
    real(RP), intent(in)  :: mask_criterion         ! criterion of hydrometeor [kg/kg]

    real(RP), intent(out) :: cldfrac(KA,IA,JA)      ! cloud fraction (0 or 1)

    real(RP) :: qhydro
    integer  :: k, i, j
    !---------------------------------------------------------------------------

    do j  = JS, JE
    do i  = IS, IE
    do k  = KS, KE
       qhydro = QTRC(k,i,j,I_hyd_QC) + QTRC(k,i,j,I_hyd_QR)
       cldfrac(k,i,j) = 0.5_RP + sign(0.5_RP,qhydro-mask_criterion)
    enddo
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_MP_kessler_cloud_fraction

  !-----------------------------------------------------------------------------
  !> Calculate Effective Radius
  subroutine ATMOS_PHY_MP_kessler_effective_radius( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       DENS0, TEMP0, QTRC0, &
       Re                   )
    use scale_atmos_hydrometeor, only: &
       N_HYD, &
       I_HC, &
       I_HR
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: DENS0(KA,IA,JA)         ! density [kg/m3]
    real(RP), intent(in)  :: TEMP0(KA,IA,JA)         ! temperature [K]
    real(RP), intent(in)  :: QTRC0(KA,IA,JA,QA_MP-1) ! hydrometeor mass concentration [kg/kg]

    real(RP), intent(out) :: Re   (KA,IA,JA,N_HYD)   ! effective radius [cm]

    real(RP), parameter :: um2cm = 100.0_RP
    !---------------------------------------------------------------------------

!OCL XFILL
    Re(:,:,:,I_HC) =   8.E-6_RP * um2cm
!OCL XFILL
    Re(:,:,:,I_HR) = 100.E-6_RP * um2cm
!OCL XFILL
    Re(:,:,:,I_HR+1:) = 0.0_RP

    return
  end subroutine ATMOS_PHY_MP_kessler_effective_radius

  !-----------------------------------------------------------------------------
  !> Calculate mass ratio of each category
  subroutine ATMOS_PHY_MP_kessler_qtrc2qhyd( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       QTRC, &
       Qe    )
    use scale_atmos_hydrometeor, only: &
       N_HYD, &
       I_HC, &
       I_HR
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: QTRC(KA,IA,JA,QA_MP-1) ! hydrometeor mass concentration [kg/kg]

    real(RP), intent(out) :: Qe(KA,IA,JA,N_HYD)   ! mass ratio of each hydrometeor [kg/kg]
    !---------------------------------------------------------------------------

!OCL XFILL
    Qe(:,:,:,I_HC) = QTRC(:,:,:,I_hyd_QC)
!OCL XFILL
    Qe(:,:,:,I_HR) = QTRC(:,:,:,I_hyd_QR)
!OCL XFILL
    Qe(:,:,:,I_HR+1:) = 0.0_RP

    return
  end subroutine ATMOS_PHY_MP_kessler_qtrc2qhyd

  !-----------------------------------------------------------------------------
  !> Calculate mass ratio of each category
  subroutine ATMOS_PHY_MP_kessler_qhyd2qtrc( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       Qe,  &
       QTRC )
    use scale_atmos_hydrometeor, only: &
       N_HYD, &
       I_HC, &
       I_HR
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in) :: Qe(KA,IA,JA,N_HYD)   ! mass ratio of each hydrometeor [kg/kg]

    real(RP), intent(out) :: QTRC(KA,IA,JA,QA_MP-1) ! hydrometeor mass concentration [kg/kg]

    !---------------------------------------------------------------------------

!OCL XFILL
    QTRC(:,:,:,I_hyd_QC) = Qe(:,:,:,I_HC)
!OCL XFILL
    QTRC(:,:,:,I_hyd_QR) = Qe(:,:,:,I_HR)

    ! ignore ice water

    return
  end subroutine ATMOS_PHY_MP_kessler_qhyd2qtrc

end module scale_atmos_phy_mp_kessler
