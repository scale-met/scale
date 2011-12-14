!-------------------------------------------------------------------------------
!> module STDIO
!!
!! @par Description
!!          Standard, common I/O module
!!
!! @author H.Tomita and SCALE developpers
!!
!! @par History
!! @li      2011-11-11 (H.Yashiro) [new]
!!
!<
!-------------------------------------------------------------------------------
module mod_stdio
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: IO_setup
  public :: IO_get_available_fid
  public :: IO_make_idstr
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

  !> Character length of system control
  integer, public, parameter :: IO_SYSCHR   = 32

  !> Character length of file name
  integer, public, parameter :: IO_FILECHR  = 128

  !> Config file ID
  integer, public, save      :: IO_FID_CONF = 10

  !> Log file ID
  integer, public, save      :: IO_FID_LOG  = 11

  !> output log or not?
  logical, public, save      :: IO_L        = .true.

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

  !> minimum available fid
  integer, private, parameter      :: IO_MINFID = 7

  !> maximum available fid
  integer, private, parameter      :: IO_MAXFID = 99  

  !> maximum digit for region/preccess specific filename
  integer, private, save           :: IO_MAXDIGIT = 6

  !> offset number of region file
  integer, private, parameter      :: IO_RGNOFFSET = 0

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> initialize file ID
  !-----------------------------------------------------------------------------
  subroutine IO_setup()
    implicit none

    logical :: IO_OUTPUT_LOGMSG     !< output log or not?

    namelist / PARAM_IO / &
       IO_MAXDIGIT,      &
       IO_OUTPUT_LOGMSG

    character(len=IO_FILECHR) :: param_fname

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- Read from argument
    if ( COMMAND_ARGUMENT_COUNT() /= 1 ) then
       write(*,*) ' xxx Program needs config file! STOP.'
       stop
    else
       call get_command_argument(1,param_fname)
    endif

    !--- Open config file till end
    IO_FID_CONF = IO_get_available_fid()
    open( IO_FID_CONF,                &
          file   = trim(param_fname), &
          form   = 'formatted',       &
          status = 'old',             &
          iostat = ierr               )

    !--- copy default value
    IO_OUTPUT_LOGMSG = IO_L

    !--- read PARAM
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_IO,iostat=ierr)

    if( ierr > 0 ) then !--- fatal error
       write(*,*) ' xxx Not appropriate names in namelist PARAM_IO . Check!'
       stop
    endif

    IO_L = IO_OUTPUT_LOGMSG

    return
  end subroutine IO_setup

  !-----------------------------------------------------------------------------
  !> search & get available file id
  !> @return fid
  !-----------------------------------------------------------------------------
  function IO_get_available_fid() result(fid)
    implicit none

    integer :: fid      !< file ID
    logical :: i_opened !< file ID is already used or not?
    !---------------------------------------------------------------------------

    do fid = IO_MINFID, IO_MAXFID
       inquire(fid,OPENED=i_opened)
       if ( .NOT. i_opened ) return
    enddo

  end function IO_get_available_fid

  !-----------------------------------------------------------------------------
  !> generate region/process specific filename
  !-----------------------------------------------------------------------------
  subroutine IO_make_idstr( &
      outstr, &
      instr,  &
      ext,    &
      rank    )
    implicit none

    character(len=*), intent(out) :: outstr !< generated string
    character(len=*), intent(in)  :: instr  !< strings
    character(len=*), intent(in)  :: ext    !< extention
    integer,          intent(in)  :: rank   !< number

    character(len=IO_SYSCHR) :: srank
    !---------------------------------------------------------------------------

    write(srank,'(I6.6)') rank + IO_RGNOFFSET

    outstr = trim(instr)//'.'//trim(ext)//trim(srank)

    return
  end subroutine IO_make_idstr

end module mod_stdio
