!-------------------------------------------------------------------------------------------
!> module NET2G vars
!!
!! @par Description
!!          Variables module for post-process of scale
!!
!! @author Team SCALE
!!
!! @par History
!! @li  2015-02-03 (R.Yoshida)  original
!!
!<
!-------------------------------------------------------------------------------------------
module mod_net2g_vars
  !-----------------------------------------------------------------------------------------
  !
  !++ used modules
  !
  !-----------------------------------------------------------------------------------------
  implicit none
  private
  !++ included parameters
#include "inc_net2g.h"
  !-----------------------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  !-----------------------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer,         public :: START_TSTEP    = 1
  integer,         public :: END_TSTEP      = 1
  integer,         public :: INC_TSTEP      = 1
  integer,         public :: DOMAIN_NUM     = 1
  integer,         public :: VCOUNT         = 1
  integer,         public :: ZCOUNT         = 0
  integer,         public :: ZSTART         = 1
  integer,         public :: TARGET_ZLEV(max_zcount) = -1
  real(DP),        public :: EXTRA_TINTERVAL = -9.999
  character(5),    public :: EXTRA_TUNIT    = ""
  character(CLNG), public :: IDIR           = "./data"
  character(CLNG), public :: ODIR           = "."
  character(CLNG), public :: CONFFILE       = "./run.conf"
  character(CSHT), public :: VNAME(max_vcount) = ""
  character(CSHT), public :: Z_LEV_TYPE     = "plev"
  character(5),    public :: DELT           = "1mn"
  character(15),   public :: STIME          = "00:00Z01JAN2000"
  character(15),   public :: FTIME          = "200001010000"
  character(15),   public :: FTIME_SAVE     = "200001010000"
  character(CLNG), public :: LOG_BASENAME   = "LOG"
  logical,         public :: LOG_ALL_OUTPUT = .false.
  logical,         public :: LOG_DBUG       = .false.
  integer,         public :: LOG_LEVEL      = 1
  logical,         public :: Z_LEV_LIST     = .false.
  logical,         public :: Z_MERGE_OUT    = .true.  ! only for slice and conv
  logical,         public :: T_MERGE_OUT    = .true.
  logical,         public :: MAPPROJ_ctl    = .true.

  character(CSHT), public :: ANALYSIS       = "ave" ! max, min, sum, ave

  integer,         public :: PRC_NUM_X
  integer,         public :: PRC_NUM_Y
  integer,         public :: TIME_STARTDATE(6)
  real(DP),        public :: HISTORY_DEFAULT_TINTERVAL
  character(CMID), public :: HISTORY_DEFAULT_BASENAME
  character(5),    public :: HISTORY_DEFAULT_TUNIT
  logical,         public :: HISTORY_DEFAULT_ZINTERP = .true.
  logical,         public :: HIST_BND = .true.

  real(DP),        public :: MPRJ_basepoint_lon
  real(DP),        public :: MPRJ_basepoint_lat
  character(CSHT), public :: MPRJ_type
  real(DP),        public :: MPRJ_LC_lat1
  real(DP),        public :: MPRJ_LC_lat2

  integer,         public :: FID_LOG        = 22
  logical,         public :: LOUT           = .false.

  character(CLNG), public :: fname_bank(max_fcount)

end module mod_net2g_vars
