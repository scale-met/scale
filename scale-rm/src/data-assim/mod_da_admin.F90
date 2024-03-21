!-------------------------------------------------------------------------------
!> module Data Assimilation admin
!!
!! @par Description
!!          Data Assimilation administrator
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_da_admin
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
  public :: DA_ADMIN_setup

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  character(len=H_SHORT), public :: DA_TYPE = 'OFF'

  logical, public :: DA_do ! do data assimilation?

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
  subroutine DA_ADMIN_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    namelist / PARAM_DA / &
       DA_TYPE

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("DA_ADMIN_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_DA,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("DA_ADMIN_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("DA_ADMIN_setup",*) 'Not appropriate names in namelist PARAM_DA. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_DA)

    !-----< module component check >-----

    LOG_NEWLINE
    LOG_INFO("DA_ADMIN_setup",*) 'Data Assimilation components '

    if ( DA_TYPE /= 'OFF' ) then
       LOG_INFO_CONT(*) 'Data Assimilation: ON, ', trim(DA_TYPE)
       DA_do = .true.
    else
       LOG_INFO_CONT(*) 'Data Assimilation: OFF'
       DA_do = .false.
    endif

    return
  end subroutine DA_ADMIN_setup

end module mod_da_admin
