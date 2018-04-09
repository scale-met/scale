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
!<
module scale_stdio
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
#include "scalelib.h"
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: IO_setup
  public :: IO_LOG_setup
  public :: IO_get_available_fid
  public :: IO_make_idstr
  public :: IO_ARG_getfname
  public :: IO_CNF_open
  public :: IO_filename_replace_setup
  public :: IO_filename_replace

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
  character(len=H_MID),  public            :: H_INSTITUTE = 'AICS/RIKEN'    !< for file header

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

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: IO_str_replace

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private, parameter :: IO_MINFID    = 10 !< minimum available fid
  integer, private, parameter :: IO_MAXFID    = 99 !< maximum available fid

  integer, private              :: IO_FILENAME_KEYWORD_NUM    = 0                       !< actual number of keywords for filename replacement
  integer, private, parameter   :: IO_FILENAME_KEYWORD_MAXNUM = 10                      !< maximum number of keywords for filename replacement
  character(len=H_MID), private :: IO_FILENAME_KEYWORD(IO_FILENAME_KEYWORD_MAXNUM) = '' !< array of keywords for filename replacement
  character(len=H_MID), private :: IO_FILENAME_VALUE(IO_FILENAME_KEYWORD_MAXNUM)   = '' !< array or values for filename replacement

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine IO_setup( &
       APPNAME,     &
       conf_name,   &
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
          write(*,*) 'xxx Not appropriate names in namelist PARAM_IO . Check!'
          stop
       endif
    end if

    call IO_filename_replace( IO_LOG_BASENAME, 'IO_LOG_BASENAME' )

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
          call IO_make_idstr(fname,trim(IO_LOG_BASENAME),'pe',myrank)
          open( unit   = IO_FID_LOG,  &
                file   = trim(fname), &
                form   = 'formatted', &
                iostat = ierr         )
          if ( ierr /= 0 ) then
             write(*,*) 'xxx File open error! :', trim(fname)
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
          write(IO_FID_LOG,*) ''
          write(IO_FID_LOG,*) '++++++ Module[STDIO] / Categ[IO] / Origin[SCALElib]'
          write(IO_FID_LOG,*) ''
          write(IO_FID_LOG,'(1x,A,I3)') '*** Open config file, FID = ', IO_FID_CONF
          write(IO_FID_LOG,'(1x,A,I3)') '*** Open log    file, FID = ', IO_FID_LOG
          write(IO_FID_LOG,'(1x,2A)')   '*** basename of log file  = ', trim(IO_LOG_BASENAME)
          write(IO_FID_LOG,*) ''
       end if

    else
       if( is_master ) write(*,*) '*** Log report is suppressed.'
    endif

    if ( IO_NML_FILENAME /= '' ) then
       if( IO_L ) write(IO_FID_LOG,*)         '*** The used config is output to the file.'
       if( IO_L ) write(IO_FID_LOG,'(1x,2A)') '*** filename of used config file   = ', trim(IO_NML_FILENAME)

       if ( is_master ) then ! write from master node only
          IO_NML     = .true. ! force on
          IO_FID_NML = IO_get_available_fid()
          open( unit   = IO_FID_NML,            &
                file   = trim(IO_NML_FILENAME), &
                form   = 'formatted',           &
                iostat = ierr                   )
          if ( ierr /= 0 ) then
             write(*,*) 'xxx File open error! :', trim(IO_NML_FILENAME)
             stop 1
          endif

          if( IO_L ) write(IO_FID_LOG,'(1x,A,I3)') '*** Open file to output used config, FID = ', IO_FID_NML

          write(IO_FID_NML,'(A)')  '################################################################################'
          write(IO_FID_NML,'(A)')  '#! configulation'
          write(IO_FID_NML,'(2A)') '#! ', trim(H_LIBNAME)
          write(IO_FID_NML,'(2A)') '#! ', trim(H_APPNAME)
          write(IO_FID_NML,'(A)')  '################################################################################'
          write(IO_FID_NML,nml=PARAM_IO)
       else
          IO_NML     = .false. ! force off
          IO_FID_NML = -1

          if( IO_L ) write(IO_FID_LOG,*) '*** The file for used config is open by the master rank'
       endif
    else
       if ( IO_NML ) then
          IO_FID_NML = IO_FID_LOG

          if( IO_L ) write(IO_FID_LOG,*) '*** The used config is output to the log.'
       else
          if( IO_L ) write(IO_FID_LOG,*) '*** The used config is not output.'
       endif
    endif

    return
  end subroutine IO_LOG_setup

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
       rank,   &
       isrgn   )
    implicit none

    character(len=H_LONG), intent(out) :: outstr !< generated string
    character(len=*),      intent(in)  :: instr  !< strings
    character(len=*),      intent(in)  :: ext    !< extention
    integer,               intent(in)  :: rank   !< number
    logical,               intent(in), optional :: isrgn !< for region? (8 digits)

    character(len=H_SHORT) :: srank
    !---------------------------------------------------------------------------

    write(srank,'(I6.6)') rank

    if ( present(isrgn) ) then
       if(isrgn) write(srank,'(I8.8)') rank-1
    endif

    outstr = trim(instr)//'.'//trim(ext)//trim(srank)

    return
  end subroutine IO_make_idstr

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
          if(is_master) write(*,*) 'xxx Program needs config file from argument! STOP.'
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
       if(is_master) write(*,*) 'xxx Failed to open config file! STOP.'
       if(is_master) write(*,*) 'xxx filename : ', trim(fname)
       stop 1
    endif

  end function IO_CNF_open

  !-----------------------------------------------------------------------------
  !> Replace the first occurrence of 'oldsub' in 'str' with 'newsub';
  !> note that 'str' will be left-adjusted no matter whether 'oldsub' is found
  subroutine IO_str_replace( &
       str,    &
       oldsub, &
       newsub, &
       pos     )
    implicit none

    character(len=*), intent(inout) :: str !< input string; output string with substring replaced
    character(len=*), intent(in) :: oldsub !< old substring to be replaced
    character(len=*), intent(in) :: newsub !< new substring
    integer, intent(out) :: pos            !< the start position of the replaced substring; if not found, return 0

    integer :: str_lent, oldsub_len, newsub_len, shift
    !---------------------------------------------------------------------------

    str = adjustl(str)
    str_lent = len_trim(str)
    oldsub_len = len(oldsub)
    newsub_len = len(newsub)

    pos = index(str, oldsub)
    if ( pos >= 1 ) then
       shift = newsub_len - oldsub_len
       if ( shift > 0 ) then
          if ( str_lent+shift > len(str) ) then
             write(*,*) "xxx The length of 'str' string is not enough for substitution."
             write(*,*) 'xxx str : ', trim(str)
             stop 1
          endif
          str(pos+oldsub_len:str_lent+shift) = adjustr(str(pos+oldsub_len:str_lent+shift))
       elseif (shift < 0) then
          str(pos+newsub_len:pos+oldsub_len-1) = repeat(' ', 0-shift)
          str(pos+newsub_len:str_lent) = adjustl(str(pos+newsub_len:str_lent))
       endif
       str(pos:pos+newsub_len-1) = newsub
    endif

    return
  end subroutine IO_str_replace

  !-----------------------------------------------------------------------------
  !> Setup keywords and values for filename replacement
  subroutine IO_filename_replace_setup( &
       keyword, &
       value    )
    implicit none

    character(len=*), intent(in) :: keyword !< key word to be replaced
    character(len=*), intent(in) :: value   !< new value

    integer :: i
    logical :: found_keyword
    !---------------------------------------------------------------------------

    found_keyword = .false.
    do i = 1, IO_FILENAME_KEYWORD_NUM
       if ( trim(keyword) == trim(IO_FILENAME_KEYWORD(i)) ) then
          found_keyword = .true.
          IO_FILENAME_VALUE(i) = trim(value)
          if( IO_L ) write(IO_FID_LOG,'(1x,4A)') '*** Filename replacement setup: ', trim(keyword), ' = ', trim(value)
          exit
       endif
    enddo

    if ( .not. found_keyword ) then
       if ( IO_FILENAME_KEYWORD_NUM >= IO_FILENAME_KEYWORD_MAXNUM ) then
          write(*,*) 'xxx The number of keywords for filename replacement reached maximum.'
          stop 1
       endif
       IO_FILENAME_KEYWORD_NUM = IO_FILENAME_KEYWORD_NUM + 1
       IO_FILENAME_KEYWORD(IO_FILENAME_KEYWORD_NUM) = trim(keyword)
       IO_FILENAME_VALUE(IO_FILENAME_KEYWORD_NUM) = trim(value)
       if( IO_L ) write(IO_FID_LOG,'(1x,4A)') '*** Filename replacement setup: ', trim(keyword), ' = ', trim(value)
    endif

    return
  end subroutine IO_filename_replace_setup

  !-----------------------------------------------------------------------------
  !> Filename replacement
  subroutine IO_filename_replace( &
       str,    &
       varname )
    implicit none

    character(len=*), intent(inout)        :: str     !< input/output filename string
    character(len=*), intent(in), optional :: varname !< variable name for printing message

    character(len=H_LONG) :: str_orig
    integer :: i, pos
    !---------------------------------------------------------------------------

    if ( IO_FILENAME_KEYWORD_NUM > 0 ) then
       str_orig = trim(str)
       if ( present(varname) ) then
          if( IO_L ) write(IO_FID_LOG,'(1x,3A)') '*** Filename replacement for ', trim(varname), '...'
       else
          if( IO_L ) write(IO_FID_LOG,'(1x,1A)') '*** Filename replacement...'
       endif

       do i = 1, IO_FILENAME_KEYWORD_NUM
          call IO_str_replace(str, trim(IO_FILENAME_KEYWORD(i)), trim(IO_FILENAME_VALUE(i)), pos)
          if ( pos == 0 ) then
             if( IO_L ) write(IO_FID_LOG,'(1x,5A)') "***   Keyword '", trim(IO_FILENAME_KEYWORD(i)), &
                        "' is not found in '", trim(str_orig), "'."
          endif
       enddo

       if ( trim(str) /= trim(str_orig) ) then
          if( IO_L ) write(IO_FID_LOG,'(1x,2A)') '***   Original value : ', trim(str_orig)
          if( IO_L ) write(IO_FID_LOG,'(1x,2A)') '***   New      value : ', trim(str)
       endif
    endif

    return
  end subroutine IO_filename_replace

end module scale_stdio
