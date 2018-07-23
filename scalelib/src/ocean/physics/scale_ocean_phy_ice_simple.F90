!-------------------------------------------------------------------------------
!> module ocean / physics / common
!!
!! @par Description
!!          ocean common module
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_ocean_phy_ice_simple
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
  public :: OCEAN_PHY_ICE_setup
  public :: OCEAN_PHY_ICE_fraction
  public :: OCEAN_PHY_ICE_adjustment
  public :: OCEAN_PHY_ICE_simple

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), public :: OCEAN_PHY_ICE_freezetemp      =  271.35_RP ! freezing temperature of sea ice [K]
  real(RP), public :: OCEAN_PHY_ICE_density         =  1000.0_RP ! density of sea ice              [kg/m3]
  real(RP), public :: OCEAN_PHY_ICE_mass_critical   =  1600.0_RP ! ice amount for fraction = 1     [kg/m2]
  real(RP), public :: OCEAN_PHY_ICE_mass_limit      = 50000.0_RP ! maximum ice amount              [kg/m2]
  real(RP), public :: OCEAN_PHY_ICE_fraction_limit  =     1.0_RP ! maximum ice fraction            [1]

  logical,  private :: OCEAN_PHY_ICE_nudging        = .false.    !< Sea ice nudging is used?
  real(DP), private :: OCEAN_PHY_ICE_nudging_tausec              !< Relaxation time [sec]
  logical,  private :: OCEAN_PHY_ICE_offline_mode   = .false.    !< Use offline mode?

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine OCEAN_PHY_ICE_setup
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_calendar, only: &
       CALENDAR_unit2sec
    use scale_file_external_input, only: &
       FILE_EXTERNAL_INPUT_file_limit, &
       FILE_EXTERNAL_INPUT_regist
    implicit none

    real(DP)               :: OCEAN_PHY_ICE_nudging_tau                                      = 0.0_DP  ! Relaxation time
    character(len=H_SHORT) :: OCEAN_PHY_ICE_nudging_tau_unit                                 = "SEC"
    character(len=H_LONG)  :: OCEAN_PHY_ICE_nudging_basename(FILE_EXTERNAL_INPUT_file_limit) = ''
    logical                :: OCEAN_PHY_ICE_nudging_enable_periodic_year                     = .false.
    logical                :: OCEAN_PHY_ICE_nudging_enable_periodic_month                    = .false.
    logical                :: OCEAN_PHY_ICE_nudging_enable_periodic_day                      = .false.
    integer                :: OCEAN_PHY_ICE_nudging_step_fixed                               = 0
    real(RP)               :: OCEAN_PHY_ICE_nudging_offset                                   = 0.0_RP
    real(RP)               :: OCEAN_PHY_ICE_nudging_defval                                  != UNDEF
    logical                :: OCEAN_PHY_ICE_nudging_check_coordinates                        = .true.
    integer                :: OCEAN_PHY_ICE_nudging_step_limit                               = 0

    namelist / PARAM_OCEAN_PHY_ICE / &
       OCEAN_PHY_ICE_density,                       &
       OCEAN_PHY_ICE_mass_critical,                 &
       OCEAN_PHY_ICE_mass_limit,                    &
       OCEAN_PHY_ICE_fraction_limit,                &
       OCEAN_PHY_ICE_nudging,                       &
       OCEAN_PHY_ICE_nudging_tau,                   &
       OCEAN_PHY_ICE_nudging_tau_unit,              &
       OCEAN_PHY_ICE_nudging_basename,              &
       OCEAN_PHY_ICE_nudging_enable_periodic_year,  &
       OCEAN_PHY_ICE_nudging_enable_periodic_month, &
       OCEAN_PHY_ICE_nudging_enable_periodic_day,   &
       OCEAN_PHY_ICE_nudging_step_fixed,            &
       OCEAN_PHY_ICE_nudging_offset,                &
       OCEAN_PHY_ICE_nudging_defval,                &
       OCEAN_PHY_ICE_nudging_check_coordinates,     &
       OCEAN_PHY_ICE_nudging_step_limit

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("OCEAN_PHY_ICE_setup",*) 'Setup'

    OCEAN_PHY_ICE_nudging_defval = UNDEF

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_PHY_ICE,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("OCEAN_PHY_ICE_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("OCEAN_PHY_ICE_setup",*) 'Not appropriate names in namelist PARAM_OCEAN_PHY_ICE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_OCEAN_PHY_ICE)

    LOG_NEWLINE
    LOG_INFO("OCEAN_PHY_ICE_setup",*) 'Ice amount for frac. = 1 [kg/m2] : ', OCEAN_PHY_ICE_mass_critical
    if ( OCEAN_PHY_ICE_nudging ) then
       call CALENDAR_unit2sec( OCEAN_PHY_ICE_nudging_tausec, OCEAN_PHY_ICE_nudging_tau, OCEAN_PHY_ICE_nudging_tau_unit )

       LOG_INFO("OCEAN_PHY_ICE_setup",*) 'Use nudging for sea ice fraction : ON'
       LOG_INFO("OCEAN_PHY_ICE_setup",*) 'Relaxation time Tau [sec]        : ', OCEAN_PHY_ICE_nudging_tausec

       if ( OCEAN_PHY_ICE_nudging_tausec == 0.0_RP ) then
          OCEAN_PHY_ICE_offline_mode = .true.
          LOG_INFO("OCEAN_PHY_ICE_setup",*) 'Tau=0 means that sea ice is completely replaced by the external data.'
       endif

       if ( OCEAN_PHY_ICE_nudging_basename(1) == '' ) then
          LOG_ERROR("OCEAN_PHY_ICE_setup",*) 'OCEAN_PHY_ICE_nudging_basename is necessary. STOP'
          call PRC_abort
       endif
    else
       LOG_INFO("OCEAN_PHY_ICE_setup",*) 'Use nudging for sea ice fraction : OFF'
    endif

    if ( OCEAN_PHY_ICE_nudging ) then
       call FILE_EXTERNAL_INPUT_regist( OCEAN_PHY_ICE_nudging_basename(:),           & ! [IN]
                                        'OCEAN_ICE_FRAC',                            & ! [IN]
                                        'XY',                                        & ! [IN]
                                        OCEAN_PHY_ICE_nudging_enable_periodic_year,  & ! [IN]
                                        OCEAN_PHY_ICE_nudging_enable_periodic_month, & ! [IN]
                                        OCEAN_PHY_ICE_nudging_enable_periodic_day,   & ! [IN]
                                        OCEAN_PHY_ICE_nudging_step_fixed,            & ! [IN]
                                        OCEAN_PHY_ICE_nudging_offset,                & ! [IN]
                                        OCEAN_PHY_ICE_nudging_defval,                & ! [IN]
                                        OCEAN_PHY_ICE_nudging_check_coordinates,     & ! [IN]
                                        OCEAN_PHY_ICE_nudging_step_limit             ) ! [IN]
    endif

    return
  end subroutine OCEAN_PHY_ICE_setup

  !-----------------------------------------------------------------------------
  subroutine OCEAN_PHY_ICE_fraction( &
       OIA, OIS, OIE, &
       OJA, OJS, OJE, &
       ICE_MASS,      &
       ICE_FRAC       )
    implicit none

    integer,  intent(in)  :: OIA, OIS, OIE
    integer,  intent(in)  :: OJA, OJS, OJE
    real(RP), intent(in)  :: ICE_MASS(OIA,OJA) ! sea ice amount        [kg/m2]
    real(RP), intent(out) :: ICE_FRAC(OIA,OJA) ! sea ice area fraction [1]

    integer :: i, j
    !---------------------------------------------------------------------------

    do j = OJS, OJE
    do i = OIS, OIE
       ICE_FRAC(i,j) = ICE_MASS(i,j) / OCEAN_PHY_ICE_mass_critical

       ICE_FRAC(i,j) = min( max( ICE_FRAC(i,j), 0.0_RP ), OCEAN_PHY_ICE_fraction_limit )
    enddo
    enddo

    return
  end subroutine OCEAN_PHY_ICE_fraction

  !-----------------------------------------------------------------------------
  subroutine OCEAN_PHY_ICE_adjustment( &
       OIA, OIS, OIE, &
       OJA, OJS, OJE, &
       calc_flag,     &
       dt,            &
       OCEAN_TEMP,    &
       ICE_TEMP,      &
       ICE_MASS       )
    use scale_const, only: &
       CONST_CL,    &
       CONST_EMELT, &
       CONST_DWATR
    use scale_ocean_dyn_slab, only: &
       OCEAN_DYN_SLAB_DEPTH
    implicit none

    integer,  intent(in)    :: OIA, OIS, OIE
    integer,  intent(in)    :: OJA, OJS, OJE
    logical,  intent(in)    :: calc_flag (OIA,OJA) ! to decide calculate or not
    real(DP), intent(in)    :: dt
    real(RP), intent(inout) :: OCEAN_TEMP(OIA,OJA) ! ocean temperature   [K]
    real(RP), intent(inout) :: ICE_TEMP  (OIA,OJA) ! sea ice temperature [K]
    real(RP), intent(inout) :: ICE_MASS  (OIA,OJA) ! sea ice amount      [kg/m2]

    real(RP) :: ICE_MASS_t
    real(RP) :: ICE_MASS_prev
    real(RP) :: factor
    real(RP) :: dt_RP

    integer  :: i, j
    !---------------------------------------------------------------------------

    dt_RP  = real(dt,kind=RP)
    factor = CONST_CL * CONST_DWATR * OCEAN_DYN_SLAB_DEPTH / dt_RP / CONST_EMELT

    do j = OJS, OJE
    do i = OIS, OIE
       if ( calc_flag(i,j) ) then
          ! update ice mass
          ICE_MASS_t    = ( OCEAN_PHY_ICE_freezetemp - OCEAN_TEMP(i,j) ) * factor ! [kg/m2/s], positive is freezing
          ICE_MASS_prev = ICE_MASS(i,j)
          ICE_MASS(i,j) = min( max( ICE_MASS(i,j) + ICE_MASS_t * dt_RP, 0.0_RP ), OCEAN_PHY_ICE_mass_limit ) ! update mass w/ limiter
          ICE_MASS_t    = ( ICE_MASS(i,j) - ICE_MASS_prev ) / dt_RP

          ! update ocean temperature
          OCEAN_TEMP(i,j) = OCEAN_TEMP(i,j) + ICE_MASS_t / factor

          ! update ice temperature
          if ( ICE_MASS(i,j) > 0.0_RP ) then
             ICE_TEMP(i,j) = ( ICE_TEMP(i,j) * ICE_MASS_prev + OCEAN_PHY_ICE_freezetemp * ICE_MASS_t * dt_RP ) &
                           / ICE_MASS(i,j)
          else
             ICE_TEMP(i,j) = OCEAN_PHY_ICE_freezetemp
          endif
       endif
    enddo
    enddo

    return
  end subroutine OCEAN_PHY_ICE_adjustment

  !-----------------------------------------------------------------------------
  !> Slab ocean model
  subroutine OCEAN_PHY_ICE_simple( &
       OIA, OIS, OIE, &
       OJA, OJS, OJE, &
       LHS,           &
       sflx_QV,       &
       sflx_rain,     &
       sflx_snow,     &
       sflx_hbalance, &
       subsfc_temp,   &
       TC_dz,         &
       ICE_TEMP,      &
       ICE_MASS,      &
       calc_flag,     &
       dt,            &
       ICE_TEMP_t,    &
       ICE_MASS_t,    &
       sflx_G,        &
       sflx_water,    &
       sflx_ice       )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       CONST_CI,    &
       CONST_EMELT
    use scale_file_external_input, only: &
       FILE_EXTERNAL_INPUT_update
    implicit none

    integer,  intent(in)  :: OIA,   OIS, OIE
    integer,  intent(in)  :: OJA,   OJS, OJE
    real(RP), intent(in)  :: LHS          (OIA,OJA) ! latent heat for sublimation [J/kg]
    real(RP), intent(in)  :: sflx_QV      (OIA,OJA) ! water vapor flux [kg/m2/s]
    real(RP), intent(in)  :: sflx_rain    (OIA,OJA) ! rain        flux [kg/m2/s]
    real(RP), intent(in)  :: sflx_snow    (OIA,OJA) ! snow        flux [kg/m2/s]
    real(RP), intent(in)  :: sflx_hbalance(OIA,OJA) ! surface heat flux balance [J/m2/s]
    real(RP), intent(in)  :: subsfc_temp  (OIA,OJA) ! subsurface temperature [K]
    real(RP), intent(in)  :: TC_dz        (OIA,OJA) ! Thermal conductance [K/m]
    real(RP), intent(in)  :: ICE_TEMP     (OIA,OJA) ! sea ice temperature [K]
    real(RP), intent(in)  :: ICE_MASS     (OIA,OJA) ! sea ice amount      [kg/m2]
    logical,  intent(in)  :: calc_flag    (OIA,OJA) ! to decide calculate or not
    real(DP), intent(in)  :: dt
    real(RP), intent(out) :: ICE_TEMP_t   (OIA,OJA) ! tendency of sea ice temperature [K/s]
    real(RP), intent(out) :: ICE_MASS_t   (OIA,OJA) ! tendency of sea ice amount      [kg/m2/s]
    real(RP), intent(out) :: SFLX_G       (OIA,OJA) ! heat         flux from sea ice to subsurface
    real(RP), intent(out) :: SFLX_water   (OIA,OJA) ! liquid water flux from sea ice to subsurface
    real(RP), intent(out) :: SFLX_ice     (OIA,OJA) ! ice    water flux from sea ice to subsurface

    real(RP) :: mass_budget
    real(RP) :: ICE_MASS_new
    real(RP) :: ICE_MASS_new2
    real(RP) :: heat_budget
    real(RP) :: heat_budget_new
    real(RP) :: ICE_TEMP_new
    real(RP) :: heating_limit
    real(RP) :: sflx_melt
    real(RP) :: dt_RP

    integer  :: i, j
    !---------------------------------------------------------------------------

    LOG_PROGRESS(*) 'ocean / physics / seaice'

    dt_RP  = real(dt,kind=RP)

    do j = OJS, OJE
    do i = OIS, OIE
       if ( calc_flag(i,j) ) then
          ! mass balance
          mass_budget     = -sflx_QV(i,j) + sflx_rain(i,j) + sflx_snow(i,j)
          ICE_MASS_new    = max( ICE_MASS(i,j) + mass_budget * dt_RP, 0.0_RP )
          mass_budget     = mass_budget - ( ICE_MASS_new - ICE_MASS(i,j) ) / dt_RP ! residual

          ! heat balance
          heat_budget     = - sflx_hbalance(i,j)                                & ! -SWD + SWU - LWD + LWU + SH + LH
                            + ( subsfc_temp(i,j) - ICE_TEMP(i,j) ) * TC_dZ(i,j) & ! heat flux from ocean
                            + sflx_rain(i,j) * CONST_EMELT                        ! rain freezing

          ! ice cooling/warming
          if ( ICE_MASS_new > 0.0_RP ) then
             heating_limit   = max( OCEAN_PHY_ICE_freezetemp - ICE_TEMP(i,j), 0.0_RP ) * CONST_CI * ICE_MASS_new / dt_RP
             heat_budget_new = max( heat_budget - heating_limit, 0.0_RP )
             ICE_TEMP_new    = ICE_TEMP(i,j) + ( heat_budget - heat_budget_new ) / CONST_CI / ICE_MASS_new * dt_RP
             heat_budget     = heat_budget_new
          else
             ICE_TEMP_new    = ICE_TEMP(i,j)
          endif

          ! ice melting
          sflx_melt       = max( heat_budget / CONST_EMELT, 0.0_RP ) ! only for positive heat flux
          ICE_MASS_new2   = max( ICE_MASS_new - sflx_melt * dt_RP, 0.0_RP )
          sflx_melt       = ( ICE_MASS_new2 - ICE_MASS_new ) / dt_RP
          heat_budget     = heat_budget - sflx_melt * CONST_EMELT
          mass_budget     = mass_budget + sflx_melt

          SFLX_G    (i,j) = heat_budget
          SFLX_water(i,j) = mass_budget
          SFLX_ice  (i,j) = 0.0_RP

          ICE_MASS_t(i,j) = ( ICE_MASS_new2 - ICE_MASS(i,j) ) / dt_RP
          ICE_TEMP_t(i,j) = ( ICE_TEMP_new  - ICE_TEMP(i,j) ) / dt_RP
       endif
    enddo
    enddo

    return
  end subroutine OCEAN_PHY_ICE_simple

end module scale_ocean_phy_ice_simple
