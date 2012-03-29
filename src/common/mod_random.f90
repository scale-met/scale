!-------------------------------------------------------------------------------
!> module random
!!
!! @par Description
!!          random number generation module
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2012-03-28 (H.Yashiro)  [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_random
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_stdio, only: &
     IO_FID_LOG, &
     IO_L
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: RANDOM_setup
  public :: RANDOM_reset
  public :: RANDOM_get

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, private, allocatable, save :: RANDOM_seedvar(:)
  integer, private,              save :: RANDOM_count
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
  subroutine RANDOM_setup
    implicit none

    integer :: nseeds
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[RANDOM]/Categ[COMMON]'

    call random_seed
    call random_seed(size=nseeds)

    allocate( RANDOM_seedvar(nseeds))

    if( IO_L ) write(IO_FID_LOG,*) '*** Array size for random seed:', nseeds

    RANDOM_count = 0

    call RANDOM_reset

    return
  end subroutine RANDOM_setup

  !-----------------------------------------------------------------------------
  subroutine RANDOM_reset
    use mod_process, only: &
       PRC_myrank
    implicit none

    integer :: time1
    real(8) :: time2
    !---------------------------------------------------------------------------

    call time(time1)
    call cpu_time(time2)
    RANDOM_count = RANDOM_count + 1

    RANDOM_seedvar(:) = time1 + int(time2*1.D6) + PRC_myrank

    call random_seed(put=RANDOM_seedvar)

    return
  end subroutine RANDOM_reset

  !-----------------------------------------------------------------------------
  subroutine RANDOM_get( var )
    implicit none

    real(8), intent(out) :: var(:,:,:)
    !---------------------------------------------------------------------------

    call random_number(var)

    return
  end subroutine RANDOM_get

end module mod_random
