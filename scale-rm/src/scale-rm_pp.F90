!-------------------------------------------------------------------------------
!> Program SCALE-RM_pp (a launcher of main routine)
!!
!! @par Description
!!          SCALE: Scalable Computing by Advanced Library and Environment
!!          Limited area model for regional weather, regional climate, and large-Eddy Simulation (LES)
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
program scalerm_pp
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_launcher
  !-----------------------------------------------------------------------------
  implicit none

  logical :: EXECUTE_PREPROCESS = .true.  ! execute preprocess tools?
  logical :: EXECUTE_MODEL      = .false. ! execute main model?

  call launcher( EXECUTE_PREPROCESS, & ! (in)
                 EXECUTE_MODEL       ) ! (in)

  stop
end program scalerm_pp
