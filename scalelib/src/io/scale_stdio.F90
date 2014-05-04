!-------------------------------------------------------------------------------
!> module STDIO
!!
!! @par Description
!!          Standard, common I/O module
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2011-11-11 (H.Yashiro)  [new]
!!
!<
module scale_stdio
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use gtool_file_h, only: &
     File_HSHORT, &
     File_HMID,   &
     File_HLONG
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
#include "scalelib.h"
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
  integer,               public, parameter :: H_SHORT         = File_HSHORT !< Character length (short=16)
  integer,               public, parameter :: H_MID           = File_HMID   !< Character length (short=64)
  integer,               public, parameter :: H_LONG          = File_HLONG  !< Character length (short=256)

  character(len=H_MID),  public            :: H_SOURCE        = 'SCALE-LES ver. '//VERSION !< for header
  character(len=H_MID),  public            :: H_INSTITUTE     = 'AICS/RIKEN'               !< for header

  integer,               public            :: IO_FID_CONF     = 7           !< Config file ID
  integer,               public            :: IO_FID_LOG      = 8           !< Log file ID

  character(len=H_LONG), public            :: IO_LOG_BASENAME = 'LOG'       !< basename of logfile
  logical,               public            :: IO_L            = .false.     !< output log or not? (this process)
  logical,               public            :: IO_LOG_SUPPRESS = .false.     !< suppress all of log output?
  logical,               public            :: IO_LOG_ALLNODE  = .false.     !< output log for each node?


  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private, parameter :: IO_MINFID    = 7  !< minimum available fid
  integer, private, parameter :: IO_MAXFID    = 99 !< maximum available fid
  integer, private, parameter :: IO_RGNOFFSET = 0  !< offset number of process file

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  !> read command argument, open config file
  subroutine IO_setup
    implicit none

    namelist / PARAM_IO / &
       H_SOURCE,        &
       H_INSTITUTE,     &
       IO_LOG_BASENAME, &
       IO_LOG_SUPPRESS, &
       IO_LOG_ALLNODE

    character(len=H_LONG) :: fname !< name of config file for each process

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- Read from argument
    if ( COMMAND_ARGUMENT_COUNT() /= 1 ) then
       write(*,*) ' xxx Program needs config file! STOP.'
       stop
    else
       call get_command_argument(1,fname)
    endif

    !--- Open config file till end
    IO_FID_CONF = IO_get_available_fid()
    open( IO_FID_CONF,          &
          file   = trim(fname), &
          form   = 'formatted', &
          status = 'old',       &
          iostat = ierr         )

    !--- read PARAM
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_IO,iostat=ierr)

    if( ierr > 0 ) then !--- fatal error
       write(*,*) ' xxx Not appropriate names in namelist PARAM_IO . Check!'
       stop
    endif

    return
  end subroutine IO_setup

  !-----------------------------------------------------------------------------
  !> search & get available file ID
  !> @return fid
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
  !> generate process specific filename
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

    character(len=H_SHORT) :: srank
    !---------------------------------------------------------------------------

    write(srank,'(I6.6)') rank + IO_RGNOFFSET

    outstr = trim(instr)//'.'//trim(ext)//trim(srank)

    return
  end subroutine IO_make_idstr

end module scale_stdio
!-------------------------------------------------------------------------------
