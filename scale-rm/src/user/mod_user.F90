!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_user
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
  public :: USER_tracer_setup
  public :: USER_setup
  public :: USER_mkinit
  public :: USER_calc_tendency
  public :: USER_update

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
  logical, private :: USER_do = .false. !< do user step?

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Config before setup of tracers
  subroutine USER_tracer_setup
    use scale_tracer, only: &
       TRACER_regist
    implicit none

    ! if you want to add tracers, call the TRACER_regist subroutine.
    ! e.g.,
!    integer, parameter     :: NQ = 1
!    integer                :: QS
!    character(len=H_SHORT) :: NAME(NQ)
!    character(len=H_MID)   :: DESC(NQ)
!    character(len=H_SHORT) :: UNIT(NQ)
!
!    data NAME (/ 'name' /)
!    data DESC (/ 'tracer name' /)
!    data UNIT (/ 'kg/kg' /)
    !---------------------------------------------------------------------------

!    call TRACER_regist( QS,   & ! [OUT]
!                        NQ,   & ! [IN]
!                        NAME, & ! [IN]
!                        DESC, & ! [IN]
!                        UNIT  ) ! [IN]

    return
  end subroutine USER_tracer_setup

  !-----------------------------------------------------------------------------
  !> Setup before setup of other components
  subroutine USER_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    namelist / PARAM_USER / &
       USER_do

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_PROGRESS(*) 'Module[DRIVER] / Categ[USER] / Origin[SCALE-RM]'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_USER,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("USER_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("USER_setup",*) 'Not appropriate names in namelist PARAM_USER. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_USER)

    LOG_NEWLINE
    LOG_INFO("USER_setup",*) 'This module is dummy.'

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Make initial state
  subroutine USER_mkinit
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_mkinit

  !-----------------------------------------------------------------------------
  !> Calculation tendency
  subroutine USER_calc_tendency
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_calc_tendency

  !-----------------------------------------------------------------------------
  !> User step
  subroutine USER_update
    use scale_prc, only: &
       PRC_abort
    implicit none
    !---------------------------------------------------------------------------

    if ( USER_do ) then
       call PRC_abort
    endif

    return
  end subroutine USER_update

end module mod_user
