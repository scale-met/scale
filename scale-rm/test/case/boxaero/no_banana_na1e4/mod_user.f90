!-------------------------------------------------------------------------------
!> module USER
!!
!! @par Description
!!          User defined module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2012-12-26 (H.Yashiro)   [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_user
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
  use scale_tracer
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: USER_config
  public :: USER_setup
  public :: USER_resume0
  public :: USER_resume
  public :: USER_step

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
  !> Config before setup of other components
  subroutine USER_config
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_regist
    implicit none

    integer, parameter :: NQ = 1
    integer :: QS
    character(len=H_SHORT) :: NAME(NQ)
    character(len=H_MID)   :: DESC(NQ)
    character(len=H_SHORT) :: UNIT(NQ)

    data NAME / 'QV' /
    data DESC / 'Specific humidity' /
    data UNIT / 'kg/kg' /

    call ATMOS_HYDROMETEOR_regist( QS,              & ! (out)
                                   1, 0, 0,         & ! (in)
                                   NAME, DESC, UNIT ) ! (in)

    return
  end subroutine USER_config

  !-----------------------------------------------------------------------------
  !> Setup before setup of other components
  subroutine USER_setup
    use scale_process, only: &
       PRC_MPIstop
    implicit none

    return
  end subroutine USER_setup

  !-----------------------------------------------------------------------------
  !> Resuming operation, before calculating tendency
  subroutine USER_resume0
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_resume0

  !-----------------------------------------------------------------------------
  !> Resuming operation
  subroutine USER_resume
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_resume

  !-----------------------------------------------------------------------------
  !> User step
  subroutine USER_step
    implicit none
    !---------------------------------------------------------------------------

    return
  end subroutine USER_step

end module mod_user
