!-------------------------------------------------------------------------------
!> module Lake admin
!!
!! @par Description
!!          Lake model administrator
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_lake_admin
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
  public :: LAKE_ADMIN_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical,                public :: LAKE_do   = .true. ! main switch for the model

  character(len=H_SHORT), public :: LAKE_DYN_TYPE = 'NONE'
                                                  ! 'OFF'
  character(len=H_SHORT), public :: LAKE_SFC_TYPE = 'NONE'

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
  subroutine LAKE_ADMIN_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

!    namelist / PARAM_LAKE / &
!       LAKE_DYN_TYPE

!    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("LAKE_ADMIN_setup",*) 'Setup'

    !--- read namelist
!!$    rewind(IO_FID_CONF)
!!$    read(IO_FID_CONF,nml=PARAM_LAKE,iostat=ierr)
!!$    if( ierr < 0 ) then !--- missing
!!$       LOG_INFO("LAKE_ADMIN_setup",*) 'Not found namelist. Default used.'
!!$    elseif( ierr > 0 ) then !--- fatal error
!!$       LOG_ERROR("LAKE_ADMIN_setup",*) 'Not appropriate names in namelist PARAM_LAKE. Check!'
!!$       call PRC_abort
!!$    endif
!!$    LOG_NML(PARAM_LAKE)

    !-----< module component check >-----

    LOG_NEWLINE
    LOG_INFO("LAKE_ADMIN_setup",*) 'Lake model components '

    if ( LAKE_DYN_TYPE /= 'OFF' .AND. LAKE_DYN_TYPE /= 'NONE' ) then
       LOG_ERROR("LAKE_ADMIN_setup",*) 'Currently, no lake model is impremented'
!!$       if ( WKMAX < 0 ) then
!!$          LOG_ERROR("LAKE_ADMIN_setup",*) 'LAKE_DYN_TYPE is set but WKMAX < 0'
!!$          call PRC_abort
!!$       end if
       LOG_INFO_CONT(*) 'Lake model : ON, ', trim(LAKE_DYN_TYPE)
       LAKE_do = .true.
    else
       LOG_INFO_CONT(*) 'Lake model : OFF'
       LAKE_do = .false.
    endif

    return
  end subroutine LAKE_ADMIN_setup

end module mod_lake_admin
