!-------------------------------------------------------------------------------
!> Program SPDDUMP
!!
!! @par Description
!!          SCALE: Scalable Computing by Advanced Library and Environment
!!          header/data veiwer for formatted data
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-11-11 (H.Yashiro) [new]
!!
!<
!-------------------------------------------------------------------------------
program prg_fio_dump
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_fileio_h, only : &
    FIO_HLONG,         &
    FIO_LITTLE_ENDIAN, &
    FIO_BIG_ENDIAN,    &
    FIO_DUMP_HEADER,   &
    FIO_DUMP_ALL
  !-----------------------------------------------------------------------------
  implicit none
  !-----------------------------------------------------------------------------
  !
  !++ parameters & variables
  !

  character(LEN=FIO_HLONG) :: fname     = ""
  integer                  :: mode      = FIO_DUMP_HEADER
  integer                  :: endian    = FIO_BIG_ENDIAN
  logical                  :: filelok   = .false.
  logical                  :: modelok   = .false.
  logical                  :: endianlok = .false.

  character(LEN=FIO_HLONG) :: argstr

  integer :: command_argument_count
  integer :: n, narg
  integer :: fid ! return from C program
  !=============================================================================

  narg = command_argument_count()

  if ( narg == 0 ) then
     write(*,*) 'Usage : fio_dump [option] [file]'
     write(*,*) '  -h show header only'
     write(*,*) '  -d dump all data   '
     write(*,*) '  -b force dump with big-endian'
     write(*,*) '  -l force dump with little-endian'
     stop
  endif

  do n = 1, narg
     call get_command_argument(n,argstr)

     if ( argstr(1:1) == '-' ) then
        select case(argstr(2:2)) 
        case('h')
           if(.not. modelok) mode = FIO_DUMP_HEADER
           modelok = .true.
        case('d')
           if(.not. modelok) mode = FIO_DUMP_ALL
           modelok = .true.
        case('b')
           if(.not. endianlok) endian = FIO_BIG_ENDIAN
           endianlok = .true.
        case('l')
           if(.not. endianlok) endian = FIO_LITTLE_ENDIAN
           endianlok = .true.
        endselect
     else
        if(.not. filelok) fname = trim(argstr)
        filelok = .true.
     endif
  enddo

  call fio_syscheck()

  call fio_register_file(fid,trim(fname))

  call fio_dump_finfo(fid,endian,mode)

end program prg_fio_dump
!-------------------------------------------------------------------------------
