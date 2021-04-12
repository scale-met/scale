!-------------------------------------------------------------------------------
!> module RANDOM
!!
!! @par Description
!!          random number generation module
!!
!! @author Team SCALE
!!
!<
#include "scalelib.h"
module scale_random
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: RANDOM_setup
  public :: RANDOM_uniform
  public :: RANDOM_normal
  public :: RANDOM_finalize

  interface RANDOM_uniform
     module procedure RANDOM_uniform_1D
     module procedure RANDOM_uniform_2D
     module procedure RANDOM_uniform_3D
  end interface

  interface RANDOM_normal
     module procedure RANDOM_normal_1D
     module procedure RANDOM_normal_2D
     module procedure RANDOM_normal_3D
  end interface RANDOM_normal

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
  logical, private :: RANDOM_FIX = .false.

  integer, private, allocatable :: RANDOM_seedvar(:)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine RANDOM_setup
    use scale_prc, only: &
       PRC_abort
    implicit none

    namelist / PARAM_RANDOM / &
       RANDOM_FIX

    integer :: nseeds, ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("RANDOM_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_RANDOM,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("RANDOM_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("RANDOM_setup",*) 'Not appropriate names in namelist PARAM_RANDOM. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_RANDOM)

    call random_seed
    call random_seed(size=nseeds)

    allocate( RANDOM_seedvar(nseeds))

    LOG_NEWLINE
    LOG_INFO("RANDOM_setup",*) 'Array size for random seed:', nseeds
    if ( RANDOM_FIX ) then
       LOG_INFO("RANDOM_setup",*) 'random seed is fixed.'
    endif

    call RANDOM_reset

    return
  end subroutine RANDOM_setup

  !-----------------------------------------------------------------------------
  !> Reset random seed
  subroutine RANDOM_reset
    use scale_prc, only: &
       PRC_myrank
    implicit none

    integer  :: time1(8) ! date and time information
    real(RP) :: time2    ! CPU time
    !---------------------------------------------------------------------------

    if ( RANDOM_FIX ) then
       ! birthday of SCALE
       time1(1) = 2011 ! The year
       time1(2) = 12   ! The month
       time1(3) = 5    ! The day of the month
       time1(4) = 9    ! Time difference with UTC in minutes
       time1(5) = 10   ! The hour of the day
       time1(6) = 20   ! The minutes of the hour
       time1(7) = 41   ! The seconds of the minute
       time1(8) = 0    ! The milliseconds of the second
       time2    = 0.0_RP
    else
       call date_and_time(values=time1)
       call cpu_time(time2)
    endif

    RANDOM_seedvar(:) = PRC_myrank                     &
                      + ( time1(1) - 1970 ) * 32140800 &
                      + time1(2) * 2678400             &
                      + time1(3) * 86400               &
                      + time1(5) * 60                  &
                      + time1(6) * 3600                &
                      + time1(7) * 60                  &
                      + time1(8)                       &
                      + int(time2*1.E6_RP)

    call random_seed(put=RANDOM_seedvar)

    return
  end subroutine RANDOM_reset

  !-----------------------------------------------------------------------------
  !> Get uniform random number ( 1D )
  subroutine RANDOM_uniform_1D( var )
    implicit none
    real(RP), intent(out) :: var(:)
    !---------------------------------------------------------------------------

    call random_number(var)

    return
  end subroutine RANDOM_uniform_1D

  !-----------------------------------------------------------------------------
  !> Get uniform random number ( 2D )
  subroutine RANDOM_uniform_2D( var )
    implicit none
    real(RP), intent(out) :: var(:,:)
    !---------------------------------------------------------------------------

    call random_number(var)

    return
  end subroutine RANDOM_uniform_2D

  !-----------------------------------------------------------------------------
  !> Get uniform random number ( 3D )
  subroutine RANDOM_uniform_3D( var )
    implicit none
    real(RP), intent(out) :: var(:,:,:)
    !---------------------------------------------------------------------------

    call random_number(var)

    return
  end subroutine RANDOM_uniform_3D

  !-----------------------------------------------------------------------------
  !> Get normal random number (1D)
  subroutine RANDOM_normal_1D( var )
    implicit none
    real(RP), intent(out) :: var(:)
    integer :: n

    n = size(var)
    call get_normal( n, var(:) )

    return
  end subroutine RANDOM_normal_1D

  !> Get normal random number (2D)
  subroutine RANDOM_normal_2D( var )
    implicit none
    real(RP), intent(out) :: var(:,:)
    integer :: n

    n = size(var)
    call get_normal( n, var(:,:) )

    return
  end subroutine RANDOM_normal_2D

  !> Get normal random number (3D)
  subroutine RANDOM_normal_3D( var )
    implicit none
    real(RP), intent(out) :: var(:,:,:)
    integer :: n

    n = size(var)
    call get_normal( n, var(:,:,:) )

    return
  end subroutine RANDOM_normal_3D

  !> finalize
  subroutine RANDOM_finalize
    implicit none
    !---------------------------------------------------------------------------

    RANDOM_FIX = .false.

    deallocate( RANDOM_seedvar )

    return
  end subroutine RANDOM_finalize

  ! private

  subroutine get_normal( n, var )
    use scale_const, only: &
       PI => CONST_PI
    implicit none
    integer,  intent(in)  :: n
    real(RP), intent(out) :: var(n)

    real(RP) :: rnd(n+1)
    real(RP) :: fact
    real(RP) :: theta
    integer :: i
    !---------------------------------------------------------------------------

    call random_number(rnd)

    !$omp parallel do &
    !$omp private(fact,theta)
    do i = 1, n/2
       fact = sqrt(-2.0_RP * log( 1.0_RP - rnd(i*2-1) ) ) ! 0 <= rnd < 1
       theta = 2.0_RP * PI * rnd(i*2)
       var(i*2-1) = fact * cos(theta)
       var(i*2  ) = fact * sin(theta)
    end do
    if ( mod(n,2) == 1 ) then
       fact = sqrt(-2.0_RP * log( 1.0_RP - rnd(n) ) )
       theta = 2.0_RP * PI * rnd(n+1)
       var(n) = fact * cos(theta)
    end if

    return
  end subroutine get_normal

end module scale_random
