!-------------------------------------------------------------------------------
!> module CONVERT driver
!!
!! @par Description
!!          administrator of convert tools (GM)
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_gm_convert
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
  public :: CONVERT_setup
  public :: CONVERT

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
  logical :: CONVERT_TOPO    = .false.
  logical :: CONVERT_LANDUSE = .false.
  logical :: CONVERT_USER    = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine CONVERT_setup
    use scale_prc, only: &
       PRC_abort
    use mod_gm_cnvtopo, only: &
       CNVTOPO_setup
!     use mod_cnvlanduse_gm, only: &
!        CNVLANDUSE_setup
!     use mod_cnvuser_gm, only: &
!        CNVUSER_setup
    implicit none

    namelist / PARAM_CONVERT / &
       CONVERT_TOPO!,    &
!        CONVERT_LANDUSE, &
!        CONVERT_USER

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("CONVERT_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_CONVERT,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("CONVERT_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("CONVERT_setup",*) 'Not appropriate names in namelist PARAM_CONVERT. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_CONVERT)

    ! set up TOPO
    if ( CONVERT_TOPO ) then
       call CNVTOPO_setup
    endif

    ! set up LANDUSE
    if ( CONVERT_LANDUSE ) then
!        call CNVLANDUSE_setup
    endif

    ! set up LANDUSE
    if ( CONVERT_USER ) then
!        call CNVUSER_setup
    endif

    return
  end subroutine CONVERT_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine CONVERT
    use scale_prc, only: &
       PRC_abort
    use mod_gm_cnvtopo, only: &
       CNVTOPO
!     use mod_cnvlanduse_gm, only: &
!        CNVLANDUSE
!     use mod_cnvuser_gm, only: &
!        CNVUSER
    implicit none
    !---------------------------------------------------------------------------

    if ( CONVERT_TOPO .OR. CONVERT_LANDUSE .OR. CONVERT_USER ) then
       LOG_NEWLINE
       LOG_PROGRESS(*) 'start convert boundary data'

       if ( CONVERT_TOPO ) then
          call CNVTOPO
       endif

       if ( CONVERT_LANDUSE ) then
!           call CNVLANDUSE
       endif

       if ( CONVERT_USER ) then
!           call CNVUSER
       endif

       LOG_PROGRESS(*) 'end   convert boundary data'
    else
       LOG_NEWLINE
       LOG_PROGRESS(*) 'skip  convert boundary data'
    endif

    return
  end subroutine CONVERT

end module mod_gm_convert
