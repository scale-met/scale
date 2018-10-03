!-------------------------------------------------------------------------------
!> module sigvars
!!
!! @par Description
!!          get SIGNAL values via C Language.
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2016-06-30 (R.Yoshida)  [new]
!!
!<
module scale_sigvars
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use iso_c_binding
  use scale_io
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: SIGVARS_Get_all

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer(c_int), public :: SIGINT  = -1
  integer(c_int), public :: SIGQUIT = -1
  integer(c_int), public :: SIGABRT = -1
  integer(c_int), public :: SIGFPE  = -1
  integer(c_int), public :: SIGSEGV = -1
  integer(c_int), public :: SIGTERM = -1

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  interface
    function get_sigint(ierr) bind(c)
      import
      integer(c_int),intent(out) :: ierr
      integer(c_int) get_sigint
    end function
  end interface

  interface
    function get_sigquit(ierr) bind(c)
      import
      integer(c_int),intent(out) :: ierr
      integer(c_int) get_sigquit
    end function
  end interface

  interface
    function get_sigabrt(ierr) bind(c)
      import
      integer(c_int),intent(out) :: ierr
      integer(c_int) get_sigabrt
    end function
  end interface

  interface
    function get_sigfpe(ierr) bind(c)
      import
      integer(c_int),intent(out) :: ierr
      integer(c_int) get_sigfpe
    end function
  end interface

  interface
    function get_sigsegv(ierr) bind(c)
      import
      integer(c_int),intent(out) :: ierr
      integer(c_int) get_sigsegv
    end function
  end interface

  interface
    function get_sigterm(ierr) bind(c)
      import
      integer(c_int),intent(out) :: ierr
      integer(c_int) get_sigterm
    end function
  end interface

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !> Get signal values
  subroutine SIGVARS_Get_all( &
       master    )
    implicit none
    logical, intent(in) :: master   ! master flag
    integer(c_int) :: ierr
    !---------------------------------------------------------------------------

    if( master ) write(*,*) ''
    if( master ) write(*,*) 'Get system signals'

    SIGINT = get_sigint(ierr)
    if ( ierr == 0 ) then
       if( master ) write(*,*) '*** signal: SIGINT  = ', SIGINT
    else
       if( master ) write(*,*) 'xxx [WARNING] Not Exist: SIGINT'
    endif

    SIGQUIT = get_sigquit(ierr)
    if ( ierr == 0 ) then
       if( master ) write(*,*) '*** signal: SIGQUIT = ', SIGQUIT
    else
       if( master ) write(*,*) 'xxx [WARNING] Not Exist: SIGQUIT'
    endif

    SIGABRT = get_sigabrt(ierr)
    if ( ierr == 0 ) then
       if( master ) write(*,*) '*** signal: SIGABRT = ', SIGABRT
    else
       if( master ) write(*,*) 'xxx [WARNING] Not Exist: SIGABRT'
    endif

    SIGFPE = get_sigfpe(ierr)
    if ( ierr == 0 ) then
       if( master ) write(*,*) '*** signal: SIGFPE  = ', SIGFPE
    else
       if( master ) write(*,*) 'xxx [WARNING] Not Exist: SIGFPE'
    endif

    SIGSEGV = get_sigsegv(ierr)
    if ( ierr == 0 ) then
       if( master ) write(*,*) '*** signal: SIGSEGV = ', SIGSEGV
    else
       if( master ) write(*,*) 'xxx [WARNING] Not Exist: SIGSEGV'
    endif

    SIGTERM = get_sigterm(ierr)
    if ( ierr == 0 ) then
       if( master ) write(*,*) '*** signal: SIGTERM = ', SIGTERM
    else
       if( master ) write(*,*) 'xxx [WARNING] Not Exist: SIGTERM'
    endif

    return
  end subroutine SIGVARS_Get_all

end module scale_sigvars
!-------------------------------------------------------------------------------
