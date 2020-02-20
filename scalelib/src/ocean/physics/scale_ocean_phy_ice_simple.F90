!-------------------------------------------------------------------------------
!> module ocean / physics / ice / simple
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

  real(RP), private :: OCEAN_PHY_ICE_mass_critical  =  1600.0_RP ! ice amount for fraction = 1     [kg/m2]
  real(RP), private :: OCEAN_PHY_ICE_mass_limit     = 50000.0_RP ! maximum ice amount              [kg/m2]
  real(RP), private :: OCEAN_PHY_ICE_fraction_limit =     1.0_RP ! maximum ice fraction            [1]
  real(RP), private :: OCEAN_PHY_ICE_dT_max         =   5.E-2_RP ! maximum delta ice temperature   [K/s]

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
       FILE_EXTERNAL_INPUT_regist
    implicit none

    real(DP)               :: OCEAN_PHY_ICE_nudging_tau                   = 0.0_DP  ! Relaxation time
    character(len=H_SHORT) :: OCEAN_PHY_ICE_nudging_tau_unit              = "SEC"
    character(len=H_LONG)  :: OCEAN_PHY_ICE_nudging_basename              = ''
    logical                :: OCEAN_PHY_ICE_nudging_basename_add_num      = .false.
    integer                :: OCEAN_PHY_ICE_nudging_number_of_files       = 1
    logical                :: OCEAN_PHY_ICE_nudging_enable_periodic_year  = .false.
    logical                :: OCEAN_PHY_ICE_nudging_enable_periodic_month = .false.
    logical                :: OCEAN_PHY_ICE_nudging_enable_periodic_day   = .false.
    integer                :: OCEAN_PHY_ICE_nudging_step_fixed            = 0
    real(RP)               :: OCEAN_PHY_ICE_nudging_offset                = 0.0_RP
    real(RP)               :: OCEAN_PHY_ICE_nudging_defval                != UNDEF
    logical                :: OCEAN_PHY_ICE_nudging_check_coordinates     = .true.
    integer                :: OCEAN_PHY_ICE_nudging_step_limit            = 0

    namelist / PARAM_OCEAN_PHY_ICE / &
       OCEAN_PHY_ICE_density,                       &
       OCEAN_PHY_ICE_mass_critical,                 &
       OCEAN_PHY_ICE_mass_limit,                    &
       OCEAN_PHY_ICE_fraction_limit,                &
       OCEAN_PHY_ICE_dT_max!,                        &
!        OCEAN_PHY_ICE_nudging,                       &
!        OCEAN_PHY_ICE_nudging_tau,                   &
!        OCEAN_PHY_ICE_nudging_tau_unit,              &
!        OCEAN_PHY_ICE_nudging_basename,              &
!        OCEAN_PHY_ICE_nudging_basename_add_num,      &
!        OCEAN_PHY_ICE_nudging_number_of_files,       &
!        OCEAN_PHY_ICE_nudging_enable_periodic_year,  &
!        OCEAN_PHY_ICE_nudging_enable_periodic_month, &
!        OCEAN_PHY_ICE_nudging_enable_periodic_day,   &
!        OCEAN_PHY_ICE_nudging_step_fixed,            &
!        OCEAN_PHY_ICE_nudging_offset,                &
!        OCEAN_PHY_ICE_nudging_defval,                &
!        OCEAN_PHY_ICE_nudging_check_coordinates,     &
!        OCEAN_PHY_ICE_nudging_step_limit

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

       if ( OCEAN_PHY_ICE_nudging_basename == '' ) then
          LOG_ERROR("OCEAN_PHY_ICE_setup",*) 'OCEAN_PHY_ICE_nudging_basename is necessary. STOP'
          call PRC_abort
       endif
    else
       LOG_INFO("OCEAN_PHY_ICE_setup",*) 'Use nudging for sea ice fraction : OFF'
    endif

    if ( OCEAN_PHY_ICE_nudging ) then
       call FILE_EXTERNAL_INPUT_regist( OCEAN_PHY_ICE_nudging_basename,              & ! [IN]
                                        OCEAN_PHY_ICE_nudging_basename_add_num,      & ! [IN]
                                        OCEAN_PHY_ICE_nudging_number_of_files,       & ! [IN]
                                        'OCEAN_ICE_FRAC',                            & ! [IN]
                                        'XY',                                        & ! [IN]
                                        OCEAN_PHY_ICE_nudging_enable_periodic_year,  & ! [IN]
                                        OCEAN_PHY_ICE_nudging_enable_periodic_month, & ! [IN]
                                        OCEAN_PHY_ICE_nudging_enable_periodic_day,   & ! [IN]
                                        OCEAN_PHY_ICE_nudging_step_fixed,            & ! [IN]
                                        OCEAN_PHY_ICE_nudging_offset,                & ! [IN]
                                        OCEAN_PHY_ICE_nudging_defval,                & ! [IN]
                                        check_coordinates = OCEAN_PHY_ICE_nudging_check_coordinates, & ! [IN]
                                        step_limit        = OCEAN_PHY_ICE_nudging_step_limit         ) ! [IN]
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

       ICE_FRAC(i,j) = min( sqrt( max( ICE_FRAC(i,j), 0.0_RP ) ), OCEAN_PHY_ICE_fraction_limit )
    enddo
    enddo

    return
  end subroutine OCEAN_PHY_ICE_fraction

  !-----------------------------------------------------------------------------
  subroutine OCEAN_PHY_ICE_adjustment( &
       OIA, OIS, OIE, &
       OJA, OJS, OJE, &
       calc_flag,     &
       OCEAN_DEPTH,   &
       OCEAN_TEMP,    &
       ICE_TEMP,      &
       ICE_MASS,      &
       MASS_FLUX,     &
       ENGI_FLUX,     &
       MASS_SUPL,     &
       ENGI_SUPL      )
    use scale_const, only: &
       DWATR => CONST_DWATR
    use scale_atmos_hydrometeor, only: &
       CV_WATER, &
       CV_ICE,   &
       LHF
    implicit none

    integer,  intent(in)    :: OIA, OIS, OIE
    integer,  intent(in)    :: OJA, OJS, OJE
    logical,  intent(in)    :: calc_flag (OIA,OJA) ! to decide calculate or not
    real(RP), intent(in)    :: OCEAN_DEPTH         ! depth of the first layer of the ocean
    real(RP), intent(inout) :: OCEAN_TEMP(OIA,OJA) ! ocean temperature   [K]
    real(RP), intent(inout) :: ICE_TEMP  (OIA,OJA) ! sea ice temperature [K]
    real(RP), intent(inout) :: ICE_MASS  (OIA,OJA) ! sea ice amount      [kg/m2]
    real(RP), intent(out)   :: MASS_FLUX(OIA,OJA)
    real(RP), intent(out)   :: ENGI_FLUX(OIA,OJA)
    real(RP), intent(out)   :: MASS_SUPL (OIA,OJA)
    real(RP), intent(out)   :: ENGI_SUPL (OIA,OJA)

    real(RP) :: C_w
    real(RP) :: ICE_MASS_frz
    real(RP) :: ICE_MASS_prev

    integer  :: i, j
    !---------------------------------------------------------------------------

    C_w = CV_WATER * DWATR * OCEAN_DEPTH

    !$omp parallel do &
    !$omp private(ICE_MASS_frz,ICE_MASS_prev)
    do j = OJS, OJE
    do i = OIS, OIE
       if (       calc_flag(i,j) &
            .and. OCEAN_TEMP(i,j) < OCEAN_PHY_ICE_freezetemp .and. ICE_MASS(i,j) < OCEAN_PHY_ICE_mass_limit ) then
          ICE_MASS_frz = C_w * ( OCEAN_PHY_ICE_freezetemp - OCEAN_TEMP(i,j) ) &
                       / ( ( CV_WATER - CV_ICE ) * OCEAN_PHY_ICE_freezetemp + LHF )
          ICE_MASS_frz = min( ICE_MASS_frz, DWATR * OCEAN_DEPTH )

          ! update ice mass
          ICE_MASS_prev = ICE_MASS(i,j)
          ICE_MASS(i,j) = ICE_MASS(i,j) + ICE_MASS_frz
          ICE_MASS(i,j) = min( ICE_MASS(i,j), OCEAN_PHY_ICE_mass_limit ) ! apply limiter
          ICE_MASS_frz = ICE_MASS(i,j) - ICE_MASS_prev

          ! update ice temperature
          ICE_TEMP(i,j) = ICE_TEMP(i,j) &
                        + ( OCEAN_PHY_ICE_freezetemp - ICE_TEMP(i,j) ) * ICE_MASS_frz / ICE_MASS(i,j)

          ! update ocean temperature
          OCEAN_TEMP(i,j) = OCEAN_TEMP(i,j) &
                          + ( CV_WATER * OCEAN_TEMP(i,j) - CV_ICE * OCEAN_PHY_ICE_freezetemp + LHF ) * ICE_MASS_frz &
                          / ( C_w - CV_WATER * ICE_MASS_frz )

          MASS_FLUX(i,j) = ICE_MASS_frz
          ENGI_FLUX(i,j) = ( CV_ICE * OCEAN_PHY_ICE_freezetemp - LHF ) * ICE_MASS_frz
          MASS_SUPL(i,j) = ICE_MASS_frz
          ENGI_SUPL(i,j) = ICE_MASS_frz * CV_WATER * OCEAN_TEMP(i,j)
       else
          MASS_FLUX(i,j) = 0.0_RP
          ENGI_FLUX(i,j) = 0.0_RP
          MASS_SUPL(i,j) = 0.0_RP
          ENGI_SUPL(i,j) = 0.0_RP
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
       iflx_water,    &
       iflx_hbalance, &
       subsfc_temp,   &
       TC_dz,         &
       ICE_TEMP,      &
       ICE_MASS,      &
       ICE_FRAC,      &
       calc_flag,     &
       dt,            &
       ICE_TEMP_t,    &
       ICE_MASS_t,    &
       sflx_G,        &
       sflx_water,    &
       sflx_RHOE      )
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_hydrometeor, only: &
       CV_WATER, &
       CV_ICE,   &
       LHF
    use scale_file_external_input, only: &
       FILE_EXTERNAL_INPUT_update
    implicit none

    integer,  intent(in)  :: OIA, OIS, OIE
    integer,  intent(in)  :: OJA, OJS, OJE
    real(RP), intent(in)  :: iflx_water   (OIA,OJA) ! input mass flux [kg/m2/s] (downward)
    real(RP), intent(in)  :: iflx_hbalance(OIA,OJA) ! input heat flux [J/m2/s]  (downward)
    real(RP), intent(in)  :: subsfc_temp  (OIA,OJA) ! subsurface temperature [K]
    real(RP), intent(in)  :: TC_dz        (OIA,OJA) ! Thermal conductance [K/m]
    real(RP), intent(in)  :: ICE_TEMP     (OIA,OJA) ! sea ice temperature [K]
    real(RP), intent(in)  :: ICE_MASS     (OIA,OJA) ! sea ice amount      [kg/m2]
    real(RP), intent(in)  :: ICE_FRAC     (OIA,OJA) ! sea ice fraction    [0-1]
    logical,  intent(in)  :: calc_flag    (OIA,OJA) ! to decide calculate or not
    real(DP), intent(in)  :: dt
    real(RP), intent(out) :: ICE_TEMP_t   (OIA,OJA) ! tendency of sea ice temperature [K/s]
    real(RP), intent(out) :: ICE_MASS_t   (OIA,OJA) ! tendency of sea ice amount      [kg/m2/s]
    real(RP), intent(out) :: SFLX_G       (OIA,OJA) ! heat flux from sea ice to subsurface
    real(RP), intent(out) :: SFLX_water   (OIA,OJA) ! mass flux from sea ice to subsurface
    real(RP), intent(out) :: SFLX_RHOE    (OIA,OJA) ! internal energy flux from sea ice to subsurface

    real(RP) :: ICE_MASS_new    ! [kg/m2]
    real(RP) :: ICE_TEMP_new    ! [K]
    real(RP) :: mass_budget     ! [kg/m2/s]
    real(RP) :: heat_budget     ! [J/m2/s]
    real(RP) :: G               ! [J/m2/s]
    real(RP) :: dM              ! [kg/m2]
    real(RP) :: dE              ! [J/m2]
    real(RP) :: M_mlt           ! [kg/m2]
    real(RP) :: dt_RP

    integer  :: i, j
    !---------------------------------------------------------------------------

    LOG_PROGRESS(*) 'ocean / physics / seaice'

    dt_RP = real(dt,kind=RP)

    !$omp parallel do &
    !$omp private(mass_budget,heat_budget,dM,dE,G,M_mlt, &
    !$omp         ICE_TEMP_new,ICE_MASS_new)
    do j = OJS, OJE
    do i = OIS, OIE
       if ( calc_flag(i,j) ) then

          ! mass change
          dM  = iflx_water(i,j) * ICE_FRAC(i,j) * dt_RP
          ICE_MASS_new = ICE_MASS(i,j) + dM

          if ( ICE_MASS_new > 0.0_RP ) then
             ! internal energy change
             G = ( subsfc_temp(i,j) - ICE_TEMP(i,j) ) * TC_dZ(i,j) ! heat flux from ocean
             dE = ( iflx_hbalance(i,j) + G ) * ICE_FRAC(i,j) * dt_RP
             ICE_TEMP_new = ICE_TEMP(i,j) &
                          + ( dE - ( CV_ICE * ICE_TEMP(i,j) - LHF ) * dM ) / ( CV_ICE * ICE_MASS_new )

             ! melting ice
             M_mlt = CV_ICE * ( ICE_TEMP_new - OCEAN_PHY_ICE_freezetemp ) * ICE_MASS_new &
                   / ( ( CV_WATER - CV_ICE ) * OCEAN_PHY_ICE_freezetemp + LHF )
             M_mlt = min( max( M_mlt, 0.0_RP ), ICE_MASS_new )

             ICE_MASS_new = ICE_MASS_new - M_mlt
             ICE_TEMP_new = ICE_TEMP_new &
                          + ( CV_ICE * ICE_TEMP_new - LHF - CV_WATER * OCEAN_PHY_ICE_freezetemp ) * M_mlt &
                          / ( CV_ICE * ICE_MASS_new )

             ! ice to ocean flux
             mass_budget = M_mlt / dt_RP
             SFLX_RHOE (i,j) = CV_WATER * OCEAN_PHY_ICE_freezetemp * mass_budget
             SFLX_G    (i,j) = - G * ICE_FRAC(i,j)
             SFLX_water(i,j) = mass_budget

          else

             ICE_MASS_new = 0.0_RP
             ICE_TEMP_new = OCEAN_PHY_ICE_freezetemp ! dummy

             SFLX_RHOE (i,j) = CV_WATER * subsfc_temp(i,j) * ICE_MASS_new / dt_RP
             SFLX_G    (i,j) = ( CV_ICE * ICE_TEMP(i,j) - LHF ) * ICE_MASS(i,j) &
                             - SFLX_RHOE(i,j) + iflx_hbalance(i,j) * ICE_FRAC(i,j)
             SFLX_water(i,j) = ICE_MASS_new ! (negative)

          endif

          ICE_MASS_t(i,j) = ( ICE_MASS_new - ICE_MASS(i,j) ) / dt_RP
          ICE_TEMP_t(i,j) = ( ICE_TEMP_new - ICE_TEMP(i,j) ) / dt_RP

       else
          SFLX_G    (i,j) = 0.0_RP
          SFLX_water(i,j) = 0.0_RP
          SFLX_RHOE (i,j) = 0.0_RP
          ICE_MASS_t(i,j) = 0.0_RP
          ICE_TEMP_t(i,j) = 0.0_RP
       endif
    enddo
    enddo

    return
  end subroutine OCEAN_PHY_ICE_simple

end module scale_ocean_phy_ice_simple
