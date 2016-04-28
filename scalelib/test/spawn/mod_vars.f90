module mod_vars
  !
  !-----------------------------------------------------------
  ! 2014/11/26: Original (Ryuji Yoshida)
  !
  !-----------------------------------------------------------

  implicit none

  !--- Public
  !-----------------------------------------------------------
  integer,  public, parameter :: SP    = 4
  integer,  public, parameter :: DP    = 8
  real(SP), public, parameter :: UNDEF = -999.99999999
  integer,  public, parameter :: FID_NML = 10   ! namelist file
  integer,  public            :: LFID           ! File ID for a Logfile

  integer, public :: nx     = 3
  integer, public :: ny     = 3
  integer, public :: nz     = 3
  integer, public :: nt     = 3
  integer, public :: DOMAIN_NUM
  integer, public :: LAUNCH_PRC
  integer, public :: PARENT_SIZE
  integer, public :: INTERCOMM_PARENT
  integer, public :: INTERCOMM_CHILD
  logical, public :: IAM_PARENT
  logical, public :: IAM_CHILD

  integer, public,  parameter :: parent = 1
  integer, public,  parameter :: child  = 2
  integer, public,  parameter :: master = 0
  integer, public,  parameter :: ypmax  = 10

  integer,  public              :: myproc
  integer,  public              :: nprocs
  integer,  public              :: dist
  integer,  public              :: source
  integer,  public              :: tag
  integer,  public, allocatable :: ireq_p(:)
  integer,  public, allocatable :: ireq_c(:)
  real(SP), public, allocatable :: var(:,:,:)
  real(SP), public, allocatable :: dummy(:,:,:)
  real(SP), public, allocatable :: send_buffer(:,:,:)
  real(SP), public, allocatable :: recv_buffer(:,:,:,:)

  character(len=128), public :: fnamelist = "run.conf"
  character(len=128), public :: flogfile  = "LOG.d00.pe000.txt"
  character(len=10),  public :: cmd       = './driver'
  character(len=3),   public :: prcnum    = '000'
  character(len=2),   public :: domnum    = '00'

  integer, public            :: sleep_loop


  !--- Private
  !-----------------------------------------------------------
  !
  !-----------------------------------------------------------
end module mod_vars
