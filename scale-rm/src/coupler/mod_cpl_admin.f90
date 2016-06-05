!-------------------------------------------------------------------------------
!> module Coupler admin
!!
!! @par Description
!!          Coupler submodel administrator
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
module mod_cpl_admin
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
  public :: CPL_ADMIN_setup
  public :: CPL_ADMIN_getscheme

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical, public :: CPL_sw ! do coupler calculation?

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine CPL_ADMIN_setup
    use mod_ocean_admin, only: &
       OCEAN_sw
    use mod_land_admin, only: &
       LAND_sw
    use mod_urban_admin, only: &
       URBAN_sw
    use scale_process, only: &
       PRC_MPIstop
    implicit none
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '++++++ Module[ADMIN] / Categ[CPL] / Origin[SCALE-RM]'

    !-----< module component check >-----

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '*** Coupler components ***'

    ! Atoms-Ocean/Land/Urban Switch
    if ( OCEAN_sw .OR. LAND_sw .OR. URBAN_sw ) then
       CPL_sw = .true.
    else
       CPL_sw = .false.
    endif

    if ( CPL_sw ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Coupler : ON'
    else
       if( IO_L ) write(IO_FID_LOG,*) '*** Coupler : OFF'
    endif

    return
  end subroutine CPL_ADMIN_setup

  !-----------------------------------------------------------------------------
  !> Get name of scheme for each component
  subroutine CPL_ADMIN_getscheme
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine CPL_ADMIN_getscheme

end module mod_cpl_admin
