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
  public :: RANDOM_get
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  include 'inc_precision.h'
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: RANDOM_reset
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private, allocatable, save :: RANDOM_seedvar(:)
  integer, private,              save :: RANDOM_count
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

    return
  end subroutine RANDOM_setup

  !-----------------------------------------------------------------------------
  subroutine RANDOM_reset
    use mod_process, only: &
       PRC_myrank
    implicit none

    integer :: time1(8)
    real(RP) :: time2
    !---------------------------------------------------------------------------

    call date_and_time(values=time1)
    call cpu_time(time2)
    RANDOM_count = RANDOM_count + 1

    RANDOM_seedvar(:) = &
         + ( time1(1) - 1970 ) * 32140800 &
         + time1(2) * 2678400 &
         + time1(3) * 86400 &
         + time1(4) * 60 &
         + time1(5) * 3600 &
         + time1(6) * 60 &
         + time1(7) &
         + int(time2*1.D6) + PRC_myrank

    call random_seed(put=RANDOM_seedvar)

    return
  end subroutine RANDOM_reset

  !-----------------------------------------------------------------------------
  subroutine RANDOM_get( var )
    implicit none

    real(RP), intent(out) :: var(:,:,:)
    !---------------------------------------------------------------------------

    call RANDOM_reset
    call random_number(var)

    return
  end subroutine RANDOM_get

end module mod_random
