!-------------------------------------------------------------------------------
!> module ocean / physics / surface thermal conductivity
!!
!! @par Description
!!          surface thermal conductivity common module
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_ocean_phy_tc
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
  public :: OCEAN_PHY_TC_seaice_setup
  public :: OCEAN_PHY_TC_seaice

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
  real(RP), public :: OCEAN_PHY_thermalcond_max    =   10.0_RP ! maximum thermal conductivity / depth [J/m2/s/K]
  real(RP), public :: OCEAN_PHY_thermalcond_seaice =    2.0_RP ! thermal conductivity of sea ice      [J/m/s/K]

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine OCEAN_PHY_TC_seaice_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    namelist / PARAM_OCEAN_PHY_TC_seaice / &
       OCEAN_PHY_thermalcond_max,   &
       OCEAN_PHY_thermalcond_seaice

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("OCEAN_PHY_TC_seaice_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_PHY_TC_seaice,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("OCEAN_PHY_TC_seaice_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("OCEAN_PHY_TC_seaice_setup",*) 'Not appropriate names in namelist PARAM_OCEAN_PHY_TC_seaice. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_OCEAN_PHY_TC_seaice)

    return
  end subroutine OCEAN_PHY_TC_seaice_setup

  !-----------------------------------------------------------------------------
  subroutine OCEAN_PHY_TC_seaice( &
       OIA, OIS, OIE, &
       OJA, OJS, OJE, &
       ICE_MASS,      &
       ICE_FRAC,      &
       mask,          &
       TC_dz          )
    use scale_const, only: &
       EPS => CONST_EPS
    use scale_ocean_phy_ice_simple, only: &
       OCEAN_PHY_ICE_density
    implicit none

    integer,  intent(in)  :: OIA, OIS, OIE
    integer,  intent(in)  :: OJA, OJS, OJE
    real(RP), intent(in)  :: ICE_MASS(OIA,OJA) ! sea ice amount               [kg/m2]
    real(RP), intent(in)  :: ICE_FRAC(OIA,OJA) ! sea ice area fraction        [1]
    logical,  intent(in)  :: mask    (OIA,OJA)
    real(RP), intent(out) :: TC_dz   (OIA,OJA) ! thermal conductivity / depth [J/m2/s/K]

    real(RP) :: ice_depth
    integer  :: i, j
    !---------------------------------------------------------------------------
    !$acc data copyin(ICE_MASS,ICE_FRAC,mask) copyout(TC_dz)

    !$acc kernels
    !$omp parallel do private(ice_depth)
    do j = OJS, OJE
    do i = OIS, OIE
       if ( mask(i,j) ) then
          ice_depth  = ICE_MASS(i,j) / OCEAN_PHY_ICE_density / max(ICE_FRAC(i,j),EPS)

          TC_dz(i,j) = OCEAN_PHY_thermalcond_seaice / max(ice_depth*0.5_RP,EPS) ! at the middle point of the layer

          TC_dz(i,j) = min( TC_dz(i,j), OCEAN_PHY_thermalcond_max )
       end if
    enddo
    enddo
    !$acc end kernels

    !$acc end data
    return
  end subroutine OCEAN_PHY_TC_seaice

end module scale_ocean_phy_tc
