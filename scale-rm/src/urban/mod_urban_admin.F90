!-------------------------------------------------------------------------------
!> module Urban admin
!!
!! @par Description
!!          Urban submodel administrator
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_urban_admin
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
  public :: URBAN_ADMIN_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,                public :: URBAN_do   = .true.  ! main switch for the model
  logical,                public :: URBAN_land = .false. ! urban is handled as a land use type

  character(len=H_SHORT), public :: URBAN_DYN_TYPE = 'NONE'
                                                   ! 'OFF'
                                                   ! 'LAND'
                                                   ! 'KUSAKA01'
  character(len=H_SHORT), public :: URBAN_SFC_TYPE = 'NONE'
                                                   ! 'KUSAKA01'

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
  subroutine URBAN_ADMIN_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    namelist / PARAM_URBAN / &
       URBAN_DYN_TYPE
    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("URBAN_ADMIN_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_URBAN,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("URBAN_ADMIN_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("URBAN_ADMIN_setup",*) 'Not appropriate names in namelist PARAM_URBAN. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_URBAN)

    !-----< module component check >-----

    LOG_NEWLINE
    LOG_INFO("URBAN_ADMIN_setup",*) 'Urban model components '

    if ( URBAN_DYN_TYPE == 'LAND' ) then
       LOG_INFO_CONT(*) 'Urban model : OFF (Land model is used for urban)'
       URBAN_do   = .false.
       URBAN_land = .true.
    else if ( URBAN_DYN_TYPE /= 'OFF' .AND. URBAN_DYN_TYPE /= 'NONE' ) then
       LOG_INFO_CONT(*) 'Urban model : ON, ', trim(URBAN_DYN_TYPE)
       URBAN_do = .true.
    else
       LOG_INFO_CONT(*) 'Urban model : OFF'
       URBAN_do = .false.
    endif

    return
  end subroutine URBAN_ADMIN_setup

end module mod_urban_admin
