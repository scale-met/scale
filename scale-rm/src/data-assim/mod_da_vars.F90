!-------------------------------------------------------------------------------
!> module Data Assimilation Variables
!!
!! @par Description
!!          Container for Data Assimilation variables
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_da_vars
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_debug
  use scale_atmos_grid_cartesC_index
  use scale_tracer

  use scale_const, only: &
    UNDEF => CONST_UNDEF
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: DA_vars_setup
  public :: DA_vars_history
  public :: DA_vars_finalize
  public :: DA_vars_monitor

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  logical, public :: DA_COMPUTE_ENS_HISTORY = .false.

  integer, public, parameter :: OBS_IN_BASENAME_MAXSIZE = 1

  character(len=H_LONG), public :: OBS_IN_BASENAME(OBS_IN_BASENAME_MAXSIZE)
  character(len=H_LONG), public :: OBS_IN_FORMAT  (OBS_IN_BASENAME_MAXSIZE)
  character(len=H_LONG), public :: OBS_IN_MASKFILE = ''

  logical, public :: POSITIVE_DEFINITE_Q    = .false.
  logical, public :: POSITIVE_DEFINITE_QHYD = .false.

  integer, public :: OBS_IN_NUM = 1

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
  subroutine DA_vars_setup
    use scale_prc, only: &
      PRC_abort
    implicit none

    integer :: ierr

    namelist / PARAM_DA_VARS / &
       DA_COMPUTE_ENS_HISTORY, &
       OBS_IN_BASENAME,        &
       OBS_IN_FORMAT,          &
       OBS_IN_MASKFILE,        &
       POSITIVE_DEFINITE_Q,    &
       POSITIVE_DEFINITE_QHYD
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("DA_vars_setup",*) 'Setup'

    OBS_IN_BASENAME(:) = ''
    OBS_IN_FORMAT  (:) = ''

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_DA_VARS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("DA_vars_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("DA_vars_setup",*) 'Not appropriate names in namelist PARAM_DA_VARS. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_DA_VARS)

    return
  end subroutine DA_vars_setup

  !-----------------------------------------------------------------------------
  !> History output set for data-assimilation variables
  subroutine DA_vars_history
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine DA_vars_history

  !-----------------------------------------------------------------------------
  !> Finalize
  subroutine DA_vars_finalize
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine DA_vars_finalize

  !-----------------------------------------------------------------------------
  !> monitor output
  subroutine DA_vars_monitor
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine DA_vars_monitor

end module mod_da_vars
