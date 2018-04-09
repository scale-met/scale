!-------------------------------------------------------------------------------
!> module ocean / physics / surface roughness length
!!
!! @par Description
!!          surface roughness length Miller92 scheme
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_ocean_phy_roughness_miller92
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: OCEAN_PHY_ROUGHNESS_miller92_setup
  public :: OCEAN_PHY_ROUGHNESS_miller92

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
  real(RP), private :: OCEAN_PHY_ROUGHNESS_miller92_CM0   = 1.0E-3_RP ! bulk coef. for U*
  real(RP), private :: OCEAN_PHY_ROUGHNESS_miller92_Z0MI  = 0.0E-0_RP ! base roughness length for momentum
  real(RP), private :: OCEAN_PHY_ROUGHNESS_miller92_Z0MR  = 1.8E-2_RP ! rough factor          for momentum
  real(RP), private :: OCEAN_PHY_ROUGHNESS_miller92_Z0MS  = 1.1E-1_RP ! smooth factor         for momentum
  real(RP), private :: OCEAN_PHY_ROUGHNESS_miller92_Z0HI  = 1.4E-5_RP ! base roughness length for heat
  real(RP), private :: OCEAN_PHY_ROUGHNESS_miller92_Z0HR  = 0.0E-0_RP ! rough factor          for heat
  real(RP), private :: OCEAN_PHY_ROUGHNESS_miller92_Z0HS  = 4.0E-1_RP ! smooth factor         for heat
  real(RP), private :: OCEAN_PHY_ROUGHNESS_miller92_Z0EI  = 1.3E-4_RP ! base roughness length for moisture
  real(RP), private :: OCEAN_PHY_ROUGHNESS_miller92_Z0ER  = 0.0E-0_RP ! rough factor          for moisture
  real(RP), private :: OCEAN_PHY_ROUGHNESS_miller92_Z0ES  = 6.2E-1_RP ! smooth factor         for moisture

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine OCEAN_PHY_ROUGHNESS_miller92_setup
    use scale_prc, only: &
       PRC_abort
    use scale_ocean_phy_roughness, only: &
       OCEAN_PHY_ROUGHNESS_setup
    implicit none

    NAMELIST / PARAM_OCEAN_PHY_ROUGHNESS_MILLER92 / &
       OCEAN_PHY_ROUGHNESS_miller92_CM0,    &
       OCEAN_PHY_ROUGHNESS_miller92_Z0MI,   &
       OCEAN_PHY_ROUGHNESS_miller92_Z0MR,   &
       OCEAN_PHY_ROUGHNESS_miller92_Z0MS,   &
       OCEAN_PHY_ROUGHNESS_miller92_Z0HI,   &
       OCEAN_PHY_ROUGHNESS_miller92_Z0HR,   &
       OCEAN_PHY_ROUGHNESS_miller92_Z0HS,   &
       OCEAN_PHY_ROUGHNESS_miller92_Z0EI,   &
       OCEAN_PHY_ROUGHNESS_miller92_Z0ER,   &
       OCEAN_PHY_ROUGHNESS_miller92_Z0ES

    integer :: ierr
    !---------------------------------------------------------------------------

    ! common setup
    call OCEAN_PHY_ROUGHNESS_setup

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_PHY_ROUGHNESS_MILLER92,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("OCEAN_PHY_ROUGHNESS_miller92_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("OCEAN_PHY_ROUGHNESS_miller92_setup",*) 'Not appropriate names in namelist PARAM_OCEAN_PHY_ROUGHNESS_MILLER92. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_OCEAN_PHY_ROUGHNESS_MILLER92)

    return
  end subroutine OCEAN_PHY_ROUGHNESS_miller92_setup

  !-----------------------------------------------------------------------------
  subroutine OCEAN_PHY_ROUGHNESS_miller92( &
       OIA, OIS, OIE, OJA, OJS, OJE, &
       Uabs,          & ! [IN]
       Z0M, Z0H, Z0E  ) ! [OUT]
    use scale_const, only: &
       GRAV => CONST_GRAV
    use scale_ocean_phy_roughness, only: &
       OCEAN_PHY_ROUGHNESS_visck, &
       OCEAN_PHY_ROUGHNESS_Ustar_min, &
       OCEAN_PHY_ROUGHNESS_Z0M_min, &
       OCEAN_PHY_ROUGHNESS_Z0H_min, &
       OCEAN_PHY_ROUGHNESS_Z0E_min
    implicit none
    integer, intent(in) :: OIA, OIS, OIE
    integer, intent(in) :: OJA, OJS, OJE

    real(RP), intent(in) :: Uabs(OIA,OJA) ! velocity at the lowest atomspheric layer [m/s]

    real(RP), intent(inout) :: Z0M(OIA,OJA) ! roughness length for momentum [m]
    real(RP), intent(inout) :: Z0H(OIA,OJA) ! roughness length for heat [m]
    real(RP), intent(inout) :: Z0E(OIA,OJA) ! roughness length for vapor [m]

    real(RP) :: Ustar

    integer  :: i, j
    !---------------------------------------------------------------------------

    !$omp parallel do &
    !$omp private(Ustar)
    do j = OJS, OJE
    do i = OIS, OIE

       Ustar = max( sqrt( OCEAN_PHY_ROUGHNESS_miller92_CM0 ) * Uabs(i,j), OCEAN_PHY_ROUGHNESS_Ustar_min )

       Z0M(i,j) = max( OCEAN_PHY_ROUGHNESS_miller92_Z0MI &
                     + OCEAN_PHY_ROUGHNESS_miller92_Z0MR / GRAV * Ustar * Ustar &
                     + OCEAN_PHY_ROUGHNESS_miller92_Z0MS * OCEAN_PHY_ROUGHNESS_visck / Ustar, &
                       OCEAN_PHY_ROUGHNESS_Z0M_min )
       Z0H(i,j) = max( OCEAN_PHY_ROUGHNESS_miller92_Z0HI &
                     + OCEAN_PHY_ROUGHNESS_miller92_Z0HR / GRAV * Ustar * Ustar &
                     + OCEAN_PHY_ROUGHNESS_miller92_Z0HS * OCEAN_PHY_ROUGHNESS_visck / Ustar, &
                       OCEAN_PHY_ROUGHNESS_Z0H_min )
       Z0E(i,j) = max( OCEAN_PHY_ROUGHNESS_miller92_Z0EI &
                     + OCEAN_PHY_ROUGHNESS_miller92_Z0ER / GRAV * Ustar * Ustar &
                     + OCEAN_PHY_ROUGHNESS_miller92_Z0ES * OCEAN_PHY_ROUGHNESS_visck / Ustar, &
                       OCEAN_PHY_ROUGHNESS_Z0E_min )

    enddo
    enddo

    return
  end subroutine OCEAN_PHY_ROUGHNESS_miller92

end module scale_ocean_phy_roughness_miller92
