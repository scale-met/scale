!-------------------------------------------------------------------------------
!> module CONVERT driver
!!
!! @par Description
!!          administrator of convert tools
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_convert
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_tracer
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
  logical :: CONVERT_2D      = .false.

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine CONVERT_setup
    use scale_prc, only: &
       PRC_abort
    use mod_cnvtopo, only: &
       CNVTOPO_setup
    use mod_cnvlanduse, only: &
       CNVLANDUSE_setup
    use mod_cnv2d, only: &
       CNV2D_setup
    implicit none

    NAMELIST / PARAM_CONVERT / &
       CONVERT_TOPO,    &
       CONVERT_LANDUSE, &
       CONVERT_2D

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
    if( CONVERT_TOPO ) then
       call CNVTOPO_setup
    end if

    ! set up LANDUSE
    if( CONVERT_LANDUSE ) then
       call CNVLANDUSE_setup
    end if

    ! set up LANDUSE
    if( CONVERT_2D ) then
       call CNV2D_setup
    end if

    return
  end subroutine CONVERT_setup

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine CONVERT
    use scale_prc, only: &
       PRC_abort
    use mod_cnvtopo, only: &
       CNVTOPO
    use mod_cnvlanduse, only: &
       CNVLANDUSE
    use mod_cnv2d, only: &
       CNV2D
    implicit none
    !---------------------------------------------------------------------------

    if ( CONVERT_TOPO .OR. CONVERT_LANDUSE .OR. CONVERT_2D ) then
       LOG_NEWLINE
       LOG_PROGRESS(*) 'start convert boundary data'

       if( CONVERT_TOPO ) then
          call CNVTOPO
       end if

       if( CONVERT_LANDUSE ) then
          call CNVLANDUSE
       end if

       if( CONVERT_2D ) then
          call CNV2D
       end if

       LOG_PROGRESS(*) 'end   convert boundary data'
    else
       LOG_NEWLINE
       LOG_PROGRESS(*) 'skip  convert boundary data'
    endif

    return
  end subroutine CONVERT

end module mod_convert
