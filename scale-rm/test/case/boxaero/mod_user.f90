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
module mod_user
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_stdio
  use scale_prof
  use scale_grid_index
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
  integer,                private, parameter :: QA = 1
  character(len=H_SHORT), private            :: QNAME(QA)
  character(len=H_MID),   private            :: QDESC(QA)
  character(len=H_SHORT), private            :: QUNIT(QA)

  data QNAME / 'QV' /

  data QDESC / 'Ratio of Water Vapor mass to total mass (Specific humidity)' /

  data QUNIT / 'kg/kg' /

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Config before setup of other components
  subroutine USER_config
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_regist, &
       I_QV
    implicit none
    !---------------------------------------------------------------------------

    call ATMOS_HYDROMETEOR_regist( I_QV,               & ! (out)
                                   1, 0, 0,            & ! (in)
                                   QNAME, QDESC, QUNIT ) ! (in)

    return
  end subroutine USER_config

  !-----------------------------------------------------------------------------
  !> Setup before setup of other components
  subroutine USER_setup
    implicit none
    !---------------------------------------------------------------------------

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
