  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ parameters of post-process for scale high-performance
  !
  !-----------------------------------------------------------------------------
  integer,  parameter :: CSHT       = 16
  integer,  parameter :: CMID       = 64
  integer,  parameter :: CLNG       = 128
  integer,  parameter :: SP         = 4
  integer,  parameter :: DP         = 8
  integer,  parameter :: max_vcount = 500
  integer,  parameter :: max_tcount = 1000
  integer,  parameter :: max_zcount = 200
  integer,  parameter :: master     = 0
  integer,  parameter :: FID_STD    = 6
  integer,  parameter :: FID_CONF   = 20
  integer,  parameter :: FID_RCNF   = 21
  integer,  parameter :: FID_LOGF   = 22
  integer,  parameter :: FID_CTL    = 23
  integer,  parameter :: FID_DAT    = 24

  integer,  parameter :: err_internal = 0
  integer,  parameter :: err_netcdf   = -1
  real(SP), parameter :: UNDEF_SP     = -9.9999D7

  integer,  parameter :: vt_2d     = 0  ! vtype index
  integer,  parameter :: vt_3d     = 1  ! vtype index
  integer,  parameter :: vt_height = 2  ! vtype index
  integer,  parameter :: vt_land   = 3  ! vtype index
  integer,  parameter :: vt_urban  = 4  ! vtype index
  integer,  parameter :: vt_tpmsk  = 5  ! vtype index

  integer,  parameter :: a_slice   = 0  ! atype index
  integer,  parameter :: a_max     = 1  ! atype index
  integer,  parameter :: a_min     = 2  ! atype index
  integer,  parameter :: a_sum     = 3  ! atype index
  integer,  parameter :: a_ave     = 4  ! atype index

  character(3) :: cmm(12)
  data cmm / "JAN", "FEB", "MAR", "APL", "MAY", "JUN", &
             "JUL", "AUG", "SEP", "OCT", "NOV", "DEC" /

