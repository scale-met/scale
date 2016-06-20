module mod_prof
  !
  !-----------------------------------------------------------
  ! 2014/11/26: Original (Ryuji Yoshida)
  !
  !-----------------------------------------------------------

  use mpi !external
  use mod_vars

  implicit none

  !--- Public
  !-----------------------------------------------------------
  public :: prof_initialize
  public :: prof_start
  public :: prof_end
  public :: prof_finalize

  !--- Private
  !-----------------------------------------------------------
  integer, private, parameter :: ml = 7
  integer, private, parameter :: mp = 10

  integer, private :: ncall(mp)

  real(DP), private :: stime(mp)
  real(DP), private :: etime(mp)
  real(DP), private :: ttime(mp)

  character(len=ml) :: pname(mp)

  integer, private :: p, pid, lastpid
 
  !-----------------------------------------------------------
  contains
  !-----------------------------------------------------------

  !-----------------------------------------------------------
  subroutine prof_initialize ()
    ! -----

    lastpid  = 0
    ncall(:) = 0
    stime(:) = 0.0D0
    etime(:) = 0.0D0
    ttime(:) = -999.9D0

    do p = 1, mp
       pname(p) = ''
    enddo

   return
  end subroutine prof_initialize
  !-----------------------------------------------------------

  !-----------------------------------------------------------
  subroutine prof_start ( inname )
    character(len=*), intent(in) :: inname

    logical :: firsttime
    ! -----

    firsttime = .true.
    do p = 1, mp
       if ( trim(inname) == trim(pname(p)) ) then
          pid       = p
          firsttime = .false.
          exit
       endif
    enddo

    if ( firsttime ) then
       lastpid = lastpid + 1
       pid = lastpid
       pname(pid) = trim(inname)
       ttime(pid) = 0.0D0
       stime(pid) = MPI_WTIME()
    else
       stime(pid) = MPI_WTIME()
    endif

   return
  end subroutine prof_start
  !-----------------------------------------------------------

  !-----------------------------------------------------------
  subroutine prof_end ( inname )
    character(len=*), intent(in) :: inname
    ! -----

    do p = 1, mp
       if ( trim(inname) == trim(pname(p)) ) then
          pid = p
          exit
       endif
    enddo

    etime(pid) = MPI_WTIME()
    ttime(pid) = ttime(pid) + (etime(pid) - stime(pid))
    ncall(pid) = ncall(pid) + 1

   return
  end subroutine prof_end
  !-----------------------------------------------------------

  !-----------------------------------------------------------
  subroutine prof_finalize ()
    ! -----

    do p = 1, lastpid
       write ( LFID, '(1X,A,A,A,F9.4,A,I4)' ) &
        "PROF: ", pname(p), "  time = ", ttime(p), "[sec]  call = ", ncall(p)
    enddo

   return
  end subroutine prof_finalize
  !-----------------------------------------------------------

  !-----------------------------------------------------------
end module mod_prof
