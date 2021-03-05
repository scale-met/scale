!-------------------------------------------------------------------------------
!> module STDIO
!!
!! @par Description
!!          Standard, common I/O module
!!
!! @author Team SCALE
!<
#include "scalelib.h"
module scale_io
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_file_h, only: &
     FILE_HSHORT, &
     FILE_HMID,   &
     FILE_HLONG
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ included parameters
  !
  character(len=*), parameter :: LIBVERSION = VERSION_MACRO

  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: IO_setup
  public :: IO_LOG_setup
  public :: IO_finalize
  public :: IO_set_universalrank
  public :: IO_get_available_fid
  public :: IO_get_fname
  public :: IO_ARG_getfname
  public :: IO_CNF_open

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer,               public, parameter :: H_SHORT     = File_HSHORT     !< Character length (short=16)
  integer,               public, parameter :: H_MID       = File_HMID       !< Character length (short=64)
  integer,               public, parameter :: H_LONG      = File_HLONG      !< Character length (short=256)

  character(len=H_MID),  public            :: H_APPNAME                     !< name of the application
  character(len=H_MID),  public            :: H_LIBNAME                     !< name and version of the library
  character(len=H_MID),  public            :: H_SOURCE                      !< for file header
  character(len=H_MID),  public            :: H_INSTITUTE = 'RIKEN'         !< for file header

  character(len=9),      public, parameter :: IO_NULLFILE   = "/dev/null"
  character(len=6),      public, parameter :: IO_STDOUT     = "STDOUT"
  integer,               public, parameter :: IO_FID_STDOUT = 6
  integer,               public            :: IO_FID_CONF   = -1            !< Config file ID
  integer,               public            :: IO_FID_LOG    = -1            !< Log file ID
  integer,               public            :: IO_FID_NML    = -1            !< Log file ID (only for output namelist)

  character(len=H_LONG), public            :: IO_LOG_BASENAME     = 'LOG'   !< basename of logfile
  character(len=H_LONG), public            :: IO_NML_FILENAME     = ''      !< filename of logfile (only for output namelist)
  logical,               public            :: IO_L                = .false. !< output log or not? (this process)
  logical,               public            :: IO_NML              = .false. !< output log or not? (for namelist, this process)
  logical,               public            :: IO_LOG_SUPPRESS     = .false. !< suppress all of log output?
  logical,               public            :: IO_LOG_NML_SUPPRESS = .false. !< suppress all of log output? (for namelist)
  logical,               public            :: IO_LOG_ALLNODE      = .false. !< output log for each node?
  integer,               public            :: IO_STEP_TO_STDOUT   = -1      !< interval for output current step to STDOUT (negative is off)

  character(len=6),      public            :: IO_UNIVERSALRANK    = "UNKNWN"!< universal rank    for error log
  character(len=6),      public            :: IO_JOBID            = "UNKNWN"!< bulk job id       for error log
  character(len=6),      public            :: IO_DOMAINID         = "UNKNWN"!< nesting domain id for error log
  character(len=6),      public            :: IO_LOCALRANK        = "UNKNWN"!< local     rank    for error log

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private, parameter :: IO_MINFID    =  10 !< minimum available fid
  integer, private, parameter :: IO_MAXFID    = 256 !< maximum available fid

  character(len=H_LONG), private            :: IO_prefix

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine IO_setup( &
       APPNAME,     &
       conf_name,   &
       prefix,      &
       allow_noconf )

    implicit none

    namelist / PARAM_IO / &
       H_SOURCE,            &
       H_INSTITUTE,         &
       IO_LOG_BASENAME,     &
       IO_NML_FILENAME,     &
       IO_LOG_SUPPRESS,     &
       IO_LOG_NML_SUPPRESS, &
       IO_LOG_ALLNODE,      &
       IO_STEP_TO_STDOUT

    character(len=*), intent(in) :: APPNAME                !< name of the application
    character(len=*), intent(in), optional :: conf_name    !< name of config file for each process
    character(len=*), intent(in), optional :: prefix       !< prefix for files
    logical,          intent(in), optional :: allow_noconf !< if true, allow program to run without configure file

    character(len=H_LONG) :: fname

    integer :: ierr
    !---------------------------------------------------------------------------

    if ( present(conf_name) ) then
       fname = conf_name
    else
       fname = IO_ARG_getfname( is_master=.true., allow_noconf=allow_noconf )
       IO_LOG_BASENAME = IO_STDOUT
    endif

    if ( present(prefix) ) then
       IO_prefix = prefix
    else
       IO_prefix = ""
    end if

    !--- Open config file till end
    IO_FID_CONF = IO_CNF_open( fname,           & ! [IN]
                               is_master=.true. ) ! [IN]

    H_APPNAME = trim(APPNAME)
    H_LIBNAME = 'SCALE Library ver. '//trim(LIBVERSION)
    H_SOURCE  = trim(APPNAME)

    !--- read PARAM
    if ( IO_FID_CONF > 0 ) then
       rewind(IO_FID_CONF)
       read(IO_FID_CONF,nml=PARAM_IO,iostat=ierr)
       if ( ierr > 0 ) then !--- fatal error
          LOG_ERROR('IO_setup',*) 'Not appropriate names in namelist PARAM_IO . Check!'
          stop
       endif
    end if

    return
  end subroutine IO_setup

  !-----------------------------------------------------------------------------
  !> Setup LOG
  subroutine IO_LOG_setup( &
       myrank,   &
       is_master )
    implicit none

    integer, intent(in) :: myrank    !< my rank ID
    logical, intent(in) :: is_master !< master process?

    namelist / PARAM_IO / &
       H_SOURCE,            &
       H_INSTITUTE,         &
       IO_LOG_BASENAME,     &
       IO_NML_FILENAME,     &
       IO_LOG_SUPPRESS,     &
       IO_LOG_NML_SUPPRESS, &
       IO_LOG_ALLNODE

    character(len=H_LONG) :: fname

    integer :: ierr
    !---------------------------------------------------------------------------

    if ( IO_LOG_SUPPRESS ) then
       IO_L = .false.
    else
       if ( is_master ) then ! master node
          IO_L = .true.
       else
          IO_L = IO_LOG_ALLNODE
       endif
    endif

    if ( IO_LOG_NML_SUPPRESS ) then
       IO_NML = .false.
    else
       IO_NML = IO_L
    endif

    if ( IO_L ) then

       !--- Open logfile
       if ( IO_LOG_BASENAME == IO_STDOUT ) then
          IO_FID_LOG = IO_FID_STDOUT
       else
          IO_FID_LOG = IO_get_available_fid()
          call IO_get_fname(fname, IO_LOG_BASENAME, rank=myrank)
          open( unit   = IO_FID_LOG,  &
                file   = trim(fname), &
                form   = 'formatted', &
                iostat = ierr         )
          if ( ierr /= 0 ) then
             LOG_ERROR('IO_LOG_setup',*) 'File open error! :', trim(fname)
             stop
          endif
       endif

       if ( IO_FID_LOG .ne. IO_FID_STDOUT ) then
          write(IO_FID_LOG,*) ''
          write(IO_FID_LOG,*) '                                               -+++++++++;              '
          write(IO_FID_LOG,*) '                                             ++++++++++++++=            '
          write(IO_FID_LOG,*) '                                           ++++++++++++++++++-          '
          write(IO_FID_LOG,*) '                                          +++++++++++++++++++++         '
          write(IO_FID_LOG,*) '                                        .+++++++++++++++++++++++        '
          write(IO_FID_LOG,*) '                                        +++++++++++++++++++++++++       '
          write(IO_FID_LOG,*) '                                       +++++++++++++++++++++++++++      '
          write(IO_FID_LOG,*) '                                      =++++++=x######+++++++++++++;     '
          write(IO_FID_LOG,*) '                                     .++++++X####XX####++++++++++++     '
          write(IO_FID_LOG,*) '         =+xxx=,               ++++  +++++=##+       .###++++++++++-    '
          write(IO_FID_LOG,*) '      ,xxxxxxxxxx-            +++++.+++++=##           .##++++++++++    '
          write(IO_FID_LOG,*) '     xxxxxxxxxxxxx           -+++x#;+++++#+              ##+++++++++.   '
          write(IO_FID_LOG,*) '    xxxxxxxx##xxxx,          ++++# +++++XX                #+++++++++-   '
          write(IO_FID_LOG,*) '   xxxxxxx####+xx+x         ++++#.++++++#                  #+++++++++   '
          write(IO_FID_LOG,*) '  +xxxxxX#X   #Xx#=        =+++#x=++++=#.                  x=++++++++   '
          write(IO_FID_LOG,*) '  xxxxxx#,    x###        .++++#,+++++#=                    x++++++++   '
          write(IO_FID_LOG,*) ' xxxxxx#.                 ++++# +++++x#                     #++++++++   '
          write(IO_FID_LOG,*) ' xxxxxx+                 ++++#-+++++=#                      #++++++++   '
          write(IO_FID_LOG,*) ',xxxxxX                 -+++XX-+++++#,                      +++++++++   '
          write(IO_FID_LOG,*) '=xxxxxX                .++++#.+++++#x                       -++++++++   '
          write(IO_FID_LOG,*) '+xxxxx=                ++++#.++++++#                        ++++++++#   '
          write(IO_FID_LOG,*) 'xxxxxx;               ++++#+=++++=#                         ++++++++#   '
          write(IO_FID_LOG,*) 'xxxxxxx              ,+++x#,+++++#-                        ;++++++++-   '
          write(IO_FID_LOG,*) '#xxxxxx              +++=# +++++xX                         ++++++++#    '
          write(IO_FID_LOG,*) 'xxxxxxxx            ++++#-+++++=#                         +++++++++X    '
          write(IO_FID_LOG,*) '-+xxxxxx+          ++++X#-++++=#.            -++;        =++++++++#     '
          write(IO_FID_LOG,*) ' #xxxxxxxx.      .+++++# +++++#x            =++++-      +++++++++XX     '
          write(IO_FID_LOG,*) ' #xxxxxxxxxx=--=++++++#.++++++#             ++++++    -+++++++++x#      '
          write(IO_FID_LOG,*) '  #+xxxxxxxxxx+++++++#x=++++=#              ++++++;=+++++++++++x#       '
          write(IO_FID_LOG,*) '  =#+xxxxxxxx+++++++##,+++++#=             =++++++++++++++++++##.       '
          write(IO_FID_LOG,*) '   X#xxxxxxxx++++++## +++++x#              ;x++++++++++++++++##.        '
          write(IO_FID_LOG,*) '    x##+xxxx+++++x## +++++=#                ##++++++++++++x##X          '
          write(IO_FID_LOG,*) '     ,###Xx+++x###x -++++=#,                .####x+++++X####.           '
          write(IO_FID_LOG,*) '       -########+   -#####x                   .X#########+.             '
          write(IO_FID_LOG,*) '           .,.      ......                         .,.                  '
          write(IO_FID_LOG,*) '                                                                        '
          write(IO_FID_LOG,*) '    .X#######     +###-        =###+           ###x         x########   '
          write(IO_FID_LOG,*) '   .#########    ######X      #######         ####        .#########x   '
          write(IO_FID_LOG,*) '   ####+++++=   X#######.    -#######x       .###;        ####x+++++.   '
          write(IO_FID_LOG,*) '   ###          ###= ####    #### x###       ####        -###.          '
          write(IO_FID_LOG,*) '  .###         ####   ###+  X###   ###X     =###.        ####           '
          write(IO_FID_LOG,*) '   ###-       ;###,        .###+   -###     ####        x##########+    '
          write(IO_FID_LOG,*) '   +####x     ####         ####     ####   ####         ###########.    '
          write(IO_FID_LOG,*) '    x######. =###          ###,     .###-  ###+        x###--------     '
          write(IO_FID_LOG,*) '      =##### X###         -###       #### ,###         ####             '
          write(IO_FID_LOG,*) '        .###=x###;        .###+       ###X ###X        ####.            '
          write(IO_FID_LOG,*) ' ###########; ###########+ ###########     ########### ,##########.     '
          write(IO_FID_LOG,*) '-###########  ,##########,  #########X      ##########  +#########      '
          write(IO_FID_LOG,*) ',,,,,,,,,,.     ,,,,,,,,,    .,,,,,,,.       .,,,,,,,,    ,,,,,,,,      '
          write(IO_FID_LOG,*) '                                                                        '
          write(IO_FID_LOG,*) '    SCALE : Scalable Computing by Advanced Library and Environment      '
          write(IO_FID_LOG,*) ''
          write(IO_FID_LOG,*) trim(H_LIBNAME)
          write(IO_FID_LOG,*) trim(H_APPNAME)


          LOG_NEWLINE
          LOG_INFO("IO_LOG_setup",*) 'Setup'

          LOG_INFO('IO_LOG_setup','(1x,A,I3)') 'Open config file, FID = ', IO_FID_CONF
          LOG_INFO('IO_LOG_setup','(1x,A,I3)') 'Open log    file, FID = ', IO_FID_LOG
          LOG_INFO('IO_LOG_setup','(1x,2A)')   'basename of log file  = ', trim(IO_LOG_BASENAME)
          LOG_NEWLINE
       end if

    else
       if( is_master ) write(*,*) '*** Log report is suppressed.'
    endif

    if ( IO_NML_FILENAME /= '' ) then
       call IO_get_fname(fname, IO_NML_FILENAME)
       LOG_INFO("IO_LOG_setup",*) 'The used configurations are output to the file.'
       LOG_INFO("IO_LOG_setup",*) 'filename of used config file   = ', trim(fname)

       if ( is_master ) then ! write from master node only
          IO_NML     = .true. ! force on
          IO_FID_NML = IO_get_available_fid()
          open( unit   = IO_FID_NML,  &
                file   = fname,       &
                form   = 'formatted', &
                iostat = ierr         )
          if ( ierr /= 0 ) then
             LOG_ERROR('IO_LOG_setup',*) 'File open error! :', trim(fname)
             stop 1
          endif

          LOG_INFO("IO_LOG_setup",'(1x,A,I3)') 'Open file to output used config, FID = ', IO_FID_NML

          write(IO_FID_NML,'(A)')  '################################################################################'
          write(IO_FID_NML,'(A)')  '#! configulation'
          write(IO_FID_NML,'(2A)') '#! ', trim(H_LIBNAME)
          write(IO_FID_NML,'(2A)') '#! ', trim(H_APPNAME)
          write(IO_FID_NML,'(A)')  '################################################################################'
          LOG_NML(PARAM_IO)
       else
          IO_NML     = .false. ! force off
          IO_FID_NML = -1

          LOG_INFO("IO_LOG_setup",*) 'The file for used config is open by the master rank'
       endif
    else
       if ( IO_NML ) then
          IO_FID_NML = IO_FID_LOG

          LOG_INFO("IO_LOG_setup",*) 'The used config is output to the log.'
       else
          LOG_INFO("IO_LOG_setup",*) 'The used config is not output.'
       endif
    endif

    write(IO_LOCALRANK,'(I6.6)') myrank

    return
  end subroutine IO_LOG_setup

  subroutine IO_finalize

    if ( IO_FID_CONF > 0 ) then
       close( IO_FID_CONF )
       IO_FID_CONF = -1
    end if

    if ( IO_FID_LOG /= IO_FID_STDOUT ) then
       close( IO_FID_LOG )
       IO_FID_LOG = -1
    end if

    if ( IO_FID_NML > 0 ) then
       close( IO_FID_NML )
       IO_FID_NML = -1
    end if

    return
  end subroutine IO_finalize

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

    if ( fid >= IO_MAXFID ) then ! reach limit
       LOG_ERROR("IO_get_available_fid",*) 'Used I/O unit number reached to the limit! STOP'
       stop
    endif

  end function IO_get_available_fid

  !-----------------------------------------------------------------------------
  !> Put for error log
  subroutine IO_set_universalrank( &
       myrank,  &
       jobid,   &
       domainid )
    implicit none

    integer, intent(in) :: myrank   !< my rank ID (universal)
    integer, intent(in) :: jobid    !< bulk job ID
    integer, intent(in) :: domainid !< nesting domain ID
    !---------------------------------------------------------------------------

    write(IO_UNIVERSALRANK,'(I6.6)') myrank
    write(IO_JOBID        ,'(I6.6)') jobid
    write(IO_DOMAINID     ,'(I6.6)') domainid

    return
  end subroutine IO_set_universalrank

  !-----------------------------------------------------------------------------
  !> generate process specific filename
  subroutine IO_get_fname( &
       outstr, &
       instr,  &
       rank,   &
       ext,    &
       len     )
    implicit none

    character(len=*), intent(out) :: outstr !< generated string
    character(len=*), intent(in)  :: instr  !< strings
    integer,          intent(in), optional :: rank
    character(len=*), intent(in), optional :: ext    !< extention
    integer,          intent(in), optional :: len

    integer :: len_
    character(len=H_SHORT) :: ext_
    character(len=8) :: srank
    character(len=6) :: fmt
    !---------------------------------------------------------------------------

    if ( present(rank) ) then
       if ( present(ext) ) then
          ext_ = ext
       else
          ext_ = "pe"
       end if
       if ( rank >= 0 ) then
          if ( present(len) ) then
             len_ = len
          else
             len_ = 6
          end if
          write(fmt,'(A,I1,A,I1,A)') '(I', len_, '.', len_, ')'
          write(srank, fmt) rank
          outstr = trim(instr)//'.'//trim(ext_)//trim(srank)
       else
          outstr = trim(instr)//'.'//trim(ext_)//"all"
       end if
    else
       outstr = instr
    end if

    if ( IO_prefix /= "" .and. outstr(1:1) /= "/" ) then
       outstr = trim(IO_prefix) // trim(outstr)
    end if

    return
  end subroutine IO_get_fname

  !-----------------------------------------------------------------------------
  !> get config filename from argument
  !> @return fname
  function IO_ARG_getfname( is_master, allow_noconf ) result(fname)
    implicit none
    logical, intent(in)           :: is_master    !> master process?
    logical, intent(in), optional :: allow_noconf !> allow no configuration file

    character(len=H_LONG) :: fname   !< filename
    logical :: allow_noconf_
    !---------------------------------------------------------------------------

    if ( COMMAND_ARGUMENT_COUNT() < 1 ) then
       allow_noconf_ = .false.
       if ( present(allow_noconf) ) allow_noconf_ = allow_noconf
       if ( .not. allow_noconf_ ) then
          if(is_master) then
             LOG_ERROR("IO_ARG_getfname",*) 'Program needs config file from argument! STOP.'
          end if
          stop
       else
          fname = IO_NULLFILE
       end if
    else
       call get_command_argument(1,fname)
    endif

  end function IO_ARG_getfname

  !-----------------------------------------------------------------------------
  !> open config file
  !> @return fid
  function IO_CNF_open( &
       fname,    &
       is_master ) &
       result(fid)
    implicit none

    character(len=*), intent(in) :: fname     !< filename
    logical,          intent(in) :: is_master !< master process?
    integer                      :: fid       !< file ID

    integer :: ierr
    !---------------------------------------------------------------------------

    fid = IO_get_available_fid()

    open( unit   = fid,         &
          file   = trim(fname), &
          form   = 'formatted', &
          status = 'old',       &
          iostat = ierr         )

    if ( ierr /= 0 ) then
       if(is_master) then
          LOG_ERROR("IO_CNF_open",*) 'Failed to open config file! STOP.'
          LOG_ERROR("IO_CNF_open",*) 'filename : ', trim(fname)
       end if
       stop
    endif

  end function IO_CNF_open

end module scale_io
