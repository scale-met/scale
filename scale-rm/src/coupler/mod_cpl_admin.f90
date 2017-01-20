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
    use scale_process, only: &
       PRC_MPIstop
    use mod_atmos_admin, only: &
       ATMOS_PHY_SF_TYPE, &
       ATMOS_sw_phy_sf
    use mod_ocean_admin, only: &
       OCEAN_sw
    use mod_land_admin, only: &
       LAND_sw
    use mod_urban_admin, only: &
       URBAN_sw
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

    ! Check Atmos_Surface setting
    if ( CPL_sw ) then
       if( IO_L ) write(IO_FID_LOG,*) '*** Coupler : ON'

       if ( ATMOS_PHY_SF_TYPE == 'COUPLE' ) then
          ! do nothing
       elseif( ATMOS_PHY_SF_TYPE == 'NONE' ) then
          if( IO_L ) write(IO_FID_LOG,*) '*** -> Surface Flux Type is forced to change from NONE to COUPLE.'
          ! overwrite
          ATMOS_PHY_SF_TYPE = 'COUPLE'
          ATMOS_sw_phy_sf   = .true.
       else
          if( IO_L ) write(IO_FID_LOG,*) '*** Surface Flux : ', trim(ATMOS_PHY_SF_TYPE)
          if( IO_L ) write(IO_FID_LOG,*) 'xxx Setting conflicts between coupler and surface flux! STOP.'
          write(*,*)                     'xxx Setting conflicts between coupler and surface flux! STOP.'
          call PRC_MPIstop
       endif
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
