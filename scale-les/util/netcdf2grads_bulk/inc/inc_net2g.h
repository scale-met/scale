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
  integer,  parameter :: max_vcount = 100
  integer,  parameter :: max_tcount = 1000
  integer,  parameter :: max_zcount = 192
  integer,  parameter :: master     = 0
  integer,  parameter :: FID_STD    = 6
  integer,  parameter :: FID_CONF   = 20
  integer,  parameter :: FID_RCNF   = 21
  integer,  parameter :: FID_LOGF   = 22
  integer,  parameter :: FID_CTL    = 23
  integer,  parameter :: FID_DAT    = 24

  real(SP), parameter :: UNDEF_SP     = -9.9999D7
  real(DP), parameter :: UNDEF_DP     = -0.9999900E+31
  real(SP), parameter :: EPS_SP       = 1.0D-6

  integer,  parameter :: err_internal = 0
  integer,  parameter :: err_known    = 1
  integer,  parameter :: err_netcdf   = -1

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
  integer,  parameter :: a_conv    = 5  ! atype index

  integer,  parameter :: c_height  = 0  ! ctype index
  integer,  parameter :: c_pres    = 1  ! ctype index

  integer,  parameter :: loc_main   = 0
  integer,  parameter :: loc_anal   = 1
  integer,  parameter :: loc_cal    = 2
  integer,  parameter :: loc_comm   = 3
  integer,  parameter :: loc_netcdf = 4
  integer,  parameter :: loc_setup  = 5
  integer,  parameter :: loc_vars   = 6

  character(3) :: cmm(12)
  data cmm / "JAN", "FEB", "MAR", "APR", "MAY", "JUN", &
             "JUL", "AUG", "SEP", "OCT", "NOV", "DEC" /

  integer,  parameter :: num_std_plev = 14
  real(SP) :: std_plev(num_std_plev)
  data std_plev / 1000, 975, 950, 925, 900, 850, 800, &
                   700, 600, 500, 400, 300, 250, 200 /

  integer,  parameter :: num_std_zlev = 14
  real(SP) :: std_zlev(num_std_zlev)
  data std_zlev /  500, 1000, 1500, 2000,  3000,  4000,  5000, &
                  6000, 7000, 8000, 9000, 10000, 11000, 12000 /

  integer,  parameter :: num_std_vname = 6
  character(CSHT) :: std_vname(num_std_vname)
  data std_vname / "PT", "PRES", "U", "V", "W", "QHYD" /
