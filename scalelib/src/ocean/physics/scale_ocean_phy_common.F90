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
module scale_ocean_phy_common
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
  public :: OCEAN_PHY_setup
  public :: OCEAN_PHY_ICEFRACTION

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
  real(RP), public :: OCEAN_PHY_DENS_seaice     =  1000.0_RP ! density of sea ice          [kg/m3]
  real(RP), public :: OCEAN_PHY_seaice_critical =   300.0_RP ! ice amount for fraction = 1 [kg/m2]
  real(RP), public :: OCEAN_PHY_seaice_max      = 50000.0_RP ! maximum ice amount          [kg/m2]
  real(RP), public :: OCEAN_PHY_icefraction_max =     1.0_RP ! maximum ice fraction        [1]

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine OCEAN_PHY_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    namelist / PARAM_OCEAN_PHY / &
       OCEAN_PHY_DENS_seaice,     &
       OCEAN_PHY_seaice_critical, &
       OCEAN_PHY_seaice_max,      &
       OCEAN_PHY_icefraction_max

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("OCEAN_PHY_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_PHY,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("OCEAN_PHY_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("OCEAN_PHY_setup",*) 'Not appropriate names in namelist PARAM_OCEAN_PHY. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_OCEAN_PHY)

    return
  end subroutine OCEAN_PHY_setup

  !-----------------------------------------------------------------------------
  subroutine OCEAN_PHY_ICEFRACTION( &
       OIA, OIS, OIE, &
       OJA, OJS, OJE, &
       ICE_MASS,      &
       ICE_FRAC       )
    use scale_const, only: &
       EPS => CONST_EPS
    implicit none

    integer,  intent(in)  :: OIA, OIS, OIE
    integer,  intent(in)  :: OJA, OJS, OJE
    real(RP), intent(in)  :: ICE_MASS(OIA,OJA) ! sea ice amount        [kg/m2]
    real(RP), intent(out) :: ICE_FRAC(OIA,OJA) ! sea ice area fraction [1]

    integer  :: i, j
    !---------------------------------------------------------------------------

    do j = OJS, OJE
    do i = OIS, OIE
       ICE_FRAC(i,j) = ICE_MASS(i,j) / OCEAN_PHY_seaice_critical

       ICE_FRAC(i,j) = min( max( ICE_FRAC(i,j), 0.0_RP ), OCEAN_PHY_icefraction_max )
    enddo
    enddo

    return
  end subroutine OCEAN_PHY_ICEFRACTION

end module scale_ocean_phy_common
