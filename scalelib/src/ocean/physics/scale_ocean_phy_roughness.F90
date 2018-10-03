!-------------------------------------------------------------------------------
!> module ocean / physics / surface roughness length
!!
!! @par Description
!!          surface roughness length common module
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_ocean_phy_roughness
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
  public :: OCEAN_PHY_ROUGHNESS_setup
  public :: OCEAN_PHY_ROUGHNESS_const_setup
  public :: OCEAN_PHY_ROUGHNESS_seaice_setup
  public :: OCEAN_PHY_ROUGHNESS_const
  public :: OCEAN_PHY_ROUGHNESS_seaice

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
  real(RP), public :: OCEAN_PHY_ROUGHNESS_visck      = 1.5E-5_RP ! kinematic viscosity [m2/s]
  real(RP), public :: OCEAN_PHY_ROUGHNESS_Ustar_min  = 1.0E-3_RP ! minimum friction velocity [m/s]
  real(RP), public :: OCEAN_PHY_ROUGHNESS_Z0M_min    = 1.0E-5_RP ! minimum roughness length for momentum [m]
  real(RP), public :: OCEAN_PHY_ROUGHNESS_Z0H_min    = 1.0E-5_RP ! minimum roughness length for heat     [m]
  real(RP), public :: OCEAN_PHY_ROUGHNESS_Z0E_min    = 1.0E-5_RP ! minimum roughness length for moisture [m]

  real(RP), public :: OCEAN_PHY_ROUGHNESS_Z0M        = 1.0E-5_RP ! constant roughness length for momentum [m]
  real(RP), public :: OCEAN_PHY_ROUGHNESS_Z0H        = 1.0E-5_RP ! constant roughness length for heat     [m]
  real(RP), public :: OCEAN_PHY_ROUGHNESS_Z0E        = 1.0E-5_RP ! constant roughness length for moisture [m]

  real(RP), public :: OCEAN_PHY_ROUGHNESS_seaice_Z0M = 2.0E-2_RP ! seaice roughness length for momentum [m]
  real(RP), public :: OCEAN_PHY_ROUGHNESS_seaice_Z0H = 2.0E-3_RP ! seaice roughness length for heat     [m]
  real(RP), public :: OCEAN_PHY_ROUGHNESS_seaice_Z0E = 2.0E-3_RP ! seaice roughness length for moisture [m]

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine OCEAN_PHY_ROUGHNESS_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    namelist / PARAM_OCEAN_PHY_ROUGHNESS / &
       OCEAN_PHY_ROUGHNESS_visck,     &
       OCEAN_PHY_ROUGHNESS_Ustar_min, &
       OCEAN_PHY_ROUGHNESS_Z0M_min,   &
       OCEAN_PHY_ROUGHNESS_Z0H_min,   &
       OCEAN_PHY_ROUGHNESS_Z0E_min

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("OCEAN_PHY_ROUGHNESS_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_PHY_ROUGHNESS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("OCEAN_PHY_ROUGHNESS_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("OCEAN_PHY_ROUGHNESS_setup",*) 'Not appropriate names in namelist PARAM_OCEAN_PHY_ROUGHNESS. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_OCEAN_PHY_ROUGHNESS)

    return
  end subroutine OCEAN_PHY_ROUGHNESS_setup

  !-----------------------------------------------------------------------------
  subroutine OCEAN_PHY_ROUGHNESS_const_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    namelist / PARAM_OCEAN_PHY_ROUGHNESS_const / &
       OCEAN_PHY_ROUGHNESS_Z0M, &
       OCEAN_PHY_ROUGHNESS_Z0H, &
       OCEAN_PHY_ROUGHNESS_Z0E

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("OCEAN_PHY_ROUGHNESS_const_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_PHY_ROUGHNESS_const,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("OCEAN_PHY_ROUGHNESS_const_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("OCEAN_PHY_ROUGHNESS_const_setup",*) 'Not appropriate names in namelist PARAM_OCEAN_PHY_ROUGHNESS_const. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_OCEAN_PHY_ROUGHNESS_const)

    return
  end subroutine OCEAN_PHY_ROUGHNESS_const_setup

  !-----------------------------------------------------------------------------
  subroutine OCEAN_PHY_ROUGHNESS_seaice_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    namelist / PARAM_OCEAN_PHY_ROUGHNESS_seaice / &
       OCEAN_PHY_ROUGHNESS_seaice_Z0M, &
       OCEAN_PHY_ROUGHNESS_seaice_Z0H, &
       OCEAN_PHY_ROUGHNESS_seaice_Z0E

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("OCEAN_PHY_ROUGHNESS_seaice_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_OCEAN_PHY_ROUGHNESS_seaice,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("OCEAN_PHY_ROUGHNESS_seaice_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("OCEAN_PHY_ROUGHNESS_seaice_setup",*) 'Not appropriate names in namelist PARAM_OCEAN_PHY_ROUGHNESS_seaice. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_OCEAN_PHY_ROUGHNESS_seaice)

    return
  end subroutine OCEAN_PHY_ROUGHNESS_seaice_setup

  !-----------------------------------------------------------------------------
  subroutine OCEAN_PHY_ROUGHNESS_const( &
       OIA, OIS, OIE, &
       OJA, OJS, OJE, &
       Z0M, Z0H, Z0E  )
    implicit none

    integer,  intent(in)  :: OIA, OIS, OIE
    integer,  intent(in)  :: OJA, OJS, OJE
    real(RP), intent(out) :: Z0M(OIA,OJA) ! roughness length for momentum [m]
    real(RP), intent(out) :: Z0H(OIA,OJA) ! roughness length for heat [m]
    real(RP), intent(out) :: Z0E(OIA,OJA) ! roughness length for vapor [m]

    integer :: i, j
    !---------------------------------------------------------------------------

    do j = OJS, OJE
    do i = OIS, OIE
       Z0M(i,j) = max( OCEAN_PHY_ROUGHNESS_Z0M, OCEAN_PHY_ROUGHNESS_Z0M_min )
       Z0H(i,j) = max( OCEAN_PHY_ROUGHNESS_Z0H, OCEAN_PHY_ROUGHNESS_Z0H_min )
       Z0E(i,j) = max( OCEAN_PHY_ROUGHNESS_Z0E, OCEAN_PHY_ROUGHNESS_Z0E_min )
    enddo
    enddo

    return
  end subroutine OCEAN_PHY_ROUGHNESS_const

  !-----------------------------------------------------------------------------
  subroutine OCEAN_PHY_ROUGHNESS_seaice( &
       OIA, OIS, OIE, &
       OJA, OJS, OJE, &
       Z0M, Z0H, Z0E  )
    implicit none

    integer,  intent(in)  :: OIA, OIS, OIE
    integer,  intent(in)  :: OJA, OJS, OJE
    real(RP), intent(out) :: Z0M(OIA,OJA) ! roughness length for momentum [m]
    real(RP), intent(out) :: Z0H(OIA,OJA) ! roughness length for heat [m]
    real(RP), intent(out) :: Z0E(OIA,OJA) ! roughness length for vapor [m]

    integer :: i, j
    !---------------------------------------------------------------------------

    do j = OJS, OJE
    do i = OIS, OIE
       Z0M(i,j) = max( OCEAN_PHY_ROUGHNESS_seaice_Z0M, OCEAN_PHY_ROUGHNESS_Z0M_min )
       Z0H(i,j) = max( OCEAN_PHY_ROUGHNESS_seaice_Z0H, OCEAN_PHY_ROUGHNESS_Z0H_min )
       Z0E(i,j) = max( OCEAN_PHY_ROUGHNESS_seaice_Z0E, OCEAN_PHY_ROUGHNESS_Z0E_min )
    enddo
    enddo

    return
  end subroutine OCEAN_PHY_ROUGHNESS_seaice

end module scale_ocean_phy_roughness
