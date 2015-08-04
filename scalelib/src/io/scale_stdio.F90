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
  public :: IO_LOG_setup
  public :: IO_get_available_fid
  public :: IO_make_idstr
  public :: IO_ARG_getfname
  public :: IO_CNF_open

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer,               public, parameter :: H_SHORT     = File_HSHORT  !< Character length (short=16)
  integer,               public, parameter :: H_MID       = File_HMID    !< Character length (short=64)
  integer,               public, parameter :: H_LONG      = File_HLONG   !< Character length (short=256)

  character(len=H_MID),  public            :: H_MODELNAME                !< name and version of the model
  character(len=H_MID),  public            :: H_LIBNAME                  !< name and version of the library
  character(len=H_MID),  public            :: H_SOURCE                   !< for file header
  character(len=H_MID),  public            :: H_INSTITUTE = 'AICS/RIKEN' !< for file header

  character(len=6),      public, parameter :: IO_STDOUT     = "STDOUT"
  integer,               public, parameter :: IO_FID_STDOUT = 6
  integer,               public            :: IO_FID_CONF   = 7             !< Config file ID
  integer,               public            :: IO_FID_LOG    = 8             !< Log file ID

  character(len=H_LONG), public            :: IO_LOG_BASENAME     = 'LOG'   !< basename of logfile
  logical,               public            :: IO_L                = .false. !< output log or not? (this process)
  logical,               public            :: IO_LNML             = .false. !< output log or not? (for namelist, this process)
  logical,               public            :: IO_LOG_SUPPRESS     = .false. !< suppress all of log output?
  logical,               public            :: IO_LOG_ALLNODE      = .false. !< output log for each node?
  logical,               public            :: IO_LOG_NML_SUPPRESS = .false. !< suppress all of log output?

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private, parameter :: IO_MINFID    = 10 !< minimum available fid
  integer, private, parameter :: IO_MAXFID    = 99 !< maximum available fid
  integer, private, parameter :: IO_RGNOFFSET = 0  !< offset number of process file

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine IO_setup( &
       MODELNAME,          &
       call_from_launcher, &
       fname_in            )
    implicit none

    namelist / PARAM_IO / &
       H_SOURCE,            &
       H_INSTITUTE,         &
       IO_LOG_BASENAME,     &
       IO_LOG_SUPPRESS,     &
       IO_LOG_ALLNODE,      &
       IO_LOG_NML_SUPPRESS

    character(len=H_MID),  intent(in) :: MODELNAME !< name of the model
    logical,               intent(in) :: call_from_launcher  !< flag to get command argument
    character(len=H_LONG), intent(in), optional :: fname_in !< name of config file for each process

    character(len=H_LONG) :: fname
    integer :: ierr
    !---------------------------------------------------------------------------

    if ( call_from_launcher ) then
       if ( present(fname_in) ) then
          fname = fname_in
       else
          write(*,*) ' xxx Not imported name of config file! STOP.'
          stop
       endif
    else
       fname = IO_ARG_getfname( is_master=.true. )
    endif

    !--- Open config file till end
    IO_FID_CONF = IO_CNF_open( fname,           & ! [IN]
                               is_master=.true. ) ! [IN]

    H_MODELNAME = trim(MODELNAME)
    H_LIBNAME   = 'SCALE Library ver. '//trim(LIBVERSION)
    H_SOURCE    = trim(MODELNAME)

    !--- read PARAM
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_IO,iostat=ierr)
    if ( ierr > 0 ) then !--- fatal error
       write(*,*) ' xxx Not appropriate names in namelist PARAM_IO . Check!'
       stop
    endif

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

    character(len=H_LONG) :: fname
    integer :: ierr
    !---------------------------------------------------------------------------

    if ( .not. IO_LOG_SUPPRESS ) then
       if ( is_master ) then ! master node
          IO_L = .true.
       else
          IO_L = IO_LOG_ALLNODE
       endif
    endif

    if ( IO_LOG_NML_SUPPRESS ) then
       IO_LNML = .false.
    else
       IO_LNML = IO_L
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
       write(IO_FID_LOG,*) trim(H_MODELNAME)
       write(IO_FID_LOG,*) ''
       write(IO_FID_LOG,*) '++++++ Module[STDIO] / Categ[IO] / Origin[SCALElib]'
       write(IO_FID_LOG,*) ''
       write(IO_FID_LOG,*) '*** Open config file, FID = ', IO_FID_CONF
       write(IO_FID_LOG,*) '*** Open log    file, FID = ', IO_FID_LOG
       write(IO_FID_LOG,*) '*** basename of log file  = ', trim(IO_LOG_BASENAME)
       write(IO_FID_LOG,*) '*** detailed log output   = ', IO_LNML

    else
       if( is_master ) write(*,*) '*** Log report is suppressed.'
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

  !-----------------------------------------------------------------------------
  !> get config filename from argument
  !> @return fname
  function IO_ARG_getfname( is_master ) result(fname)
    implicit none

    logical, intent(in)   :: is_master !< master process?
    character(len=H_LONG) :: fname     !< filename
    !---------------------------------------------------------------------------

    if ( COMMAND_ARGUMENT_COUNT() < 1 ) then
       if(is_master) write(*,*) 'xxx Program needs config file from argument! STOP.'
       stop
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
       stop
    endif

  end function IO_CNF_open

end module scale_stdio
!-------------------------------------------------------------------------------
