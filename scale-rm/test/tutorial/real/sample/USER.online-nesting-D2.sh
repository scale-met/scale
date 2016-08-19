#!/bin/bash

#################################################
#
# &PARAM_NEST            (ALL configs)
# &PARAM_TIME            (ALL configs)
# &PARAM_HISTORY         (run config)
# &PARAM_ATMOS_BOUNDARY  (run config)
# &PARAM_ATMOS_REFSTATE  (run config)
#
#################################################

NUM_DOMAIN=2 # set number of domains

RUN_DATE_YEAR=2007
RUN_DATE_MON=7
RUN_DATE_DAY=14
RUN_DATE_HOUR=18
RUN_DATE_MIN=0
RUN_DATE_SEC=0
RUN_DATE_MSEC=0

TIME_DURATION="6.0"
TIME_DURATION_UNIT="HOUR" # unit: DAY / HOUR / MIN / SEC

TIME_DT_RESTART="10800.0"   # unit: SEC only
TIME_DT_BOUNDARY="21600.0"  # unit: SEC only
TIME_DT_REFSTATE="10800.0"  # unit: SEC only
TIME_DT_HISTORY_2D="3600.0" # unit: SEC only
TIME_DT_HISTORY_3D="3600.0" # unit: SEC only

TIME_DT=(               "90.0"  "18.0" ) # required parameters for each domain - unit: SEC only
TIME_DT_ATMOS_DYN=(     "45.0"   "9.0" ) # required parameters for each domain - unit: SEC only
TIME_DT_ATMOS_PHY_MP=(  "90.0"  "18.0" ) # required parameters for each domain - unit: SEC only
TIME_DT_ATMOS_PHY_RD=( "900.0" "180.0" ) # required parameters for each domain - unit: SEC only
TIME_DT_ATMOS_PHY_SF=(  "90.0"  "18.0" ) # required parameters for each domain - unit: SEC only
TIME_DT_ATMOS_PHY_TB=(  "90.0"  "18.0" ) # required parameters for each domain - unit: SEC only
TIME_DT_OCEAN=(        "450.0"  "90.0" ) # required parameters for each domain - unit: SEC only
TIME_DT_LAND=(         "450.0"  "90.0" ) # required parameters for each domain - unit: SEC only
TIME_DT_URBAN=(        "450.0"  "90.0" ) # required parameters for each domain - unit: SEC only

#################################################
#
# &PARAM_PRC          (ALL configs)
# &PARAM_INDEX        (ALL configs)
# &PARAM_LAND_INDEX   (ALL configs)
# &PARAM_LAND_GRID    (ALL configs)
# &PARAM_URBAN_INDEX  (ALL configs)
# &PARAM_URBAN_GRID   (ALL configs)
# &PARAM_GRID         (ALL configs)
# &PARAM_MAPPROJ      (ALL configs)
#
#################################################

PRC_NUM_X=( 2 4 ) # required parameters for each domain
PRC_NUM_Y=( 2 4 ) # required parameters for each domain
 
KMAX=( 36 36 ) # required parameters for each domain
IMAX=( 45 32 ) # required parameters for each domain
JMAX=( 45 32 ) # required parameters for each domain

LKMAX=5
LDZ=( "0.05" "0.15" "0.30" "0.50" "1.00" ) # required parameters for LKMAX

UKMAX=5
UDZ=( "0.01" "0.01" "0.03" "0.05" "0.10" ) # required parameters for UKMAX

DX=( "20000.0" "4000.0" ) # required parameters for each domain
DY=( "20000.0" "4000.0" ) # required parameters for each domain

DEF_Z=( 
"FZ(:) =    80.841,   248.821,   429.882,   625.045,   835.409,  1062.158,
          1306.565,  1570.008,  1853.969,  2160.047,  2489.963,  2845.575,
          3228.883,  3642.044,  4087.384,  4567.409,  5084.820,  5642.530,
          6243.676,  6891.642,  7590.074,  8342.904,  9154.367, 10029.028,
         10971.815, 11988.030, 13083.390, 14264.060, 15536.685, 16908.430,
         18387.010, 19980.750, 21698.615, 23550.275, 25546.155, 28113.205,"

"FZ(:) =    80.841,   248.821,   429.882,   625.045,   835.409,  1062.158,
          1306.565,  1570.008,  1853.969,  2160.047,  2489.963,  2845.575,
          3228.883,  3642.044,  4087.384,  4567.409,  5084.820,  5642.530,
          6243.676,  6891.642,  7590.074,  8342.904,  9154.367, 10029.028,
         10971.815, 11988.030, 13083.390, 14264.060, 15536.685, 16908.430,
         18387.010, 19980.750, 21698.615, 23550.275, 25546.155, 28113.205,"
) # required num. of parameters for each domain

BUFFER_DZ=( "5000.0"   "5000.0"  ) # required parameters for each domain
BUFFER_DX=( "400000.0" "80000.0" ) # required parameters for each domain
BUFFER_DY=( "400000.0" "80000.0" ) # required parameters for each domain

MPRJ_BASEPOINT_LON="135.220404"
MPRJ_BASEPOINT_LAT="34.653396"
MPRJ_TYPE="LC"
MPRJ_LC_LAT1="30.0"
MPRJ_LC_LAT2="40.0"

#################################################
#
# &PARAM_TRACER  (run config)
# &PARAM_ATMOS   (run config)
# &PARAM_OCEAN   (run config)
# &PARAM_LAND    (run config)
# &PARAM_URBAN   (run config)
#
#################################################

ATMOS_DYN_TYPE=(    "HEVI"     "HEVI"     ) # required parameters for each domain
ATMOS_PHY_MP_TYPE=( "TOMITA08" "TOMITA08" ) # required parameters for each domain
ATMOS_PHY_RD_TYPE=( "MSTRNX"   "MSTRNX"   ) # required parameters for each domain
ATMOS_PHY_SF_TYPE=( "COUPLE"   "COUPLE"   ) # required parameters for each domain
ATMOS_PHY_TB_TYPE=( "MYNN"     "MYNN"     ) # required parameters for each domain

OCEAN_TYPE=( "CONST" "CONST" ) # required parameters for each domain
LAND_TYPE=(  "SLAB"  "SLAB"  ) # required parameters for each domain
URBAN_TYPE=( "SLC"   "SLC"   ) # required parameters for each domain

#################################################
#
# &HISTITEM (run config)
#
#################################################

HIST_ITEMS_SNAPSHOT_2D=(
  "MSLP"
)
HIST_ITEMS_SNAPSHOT_3D=(
  "U" "V"
)
HIST_ITEMS_AVERAGE_2D=(
  "PREC"
)
HIST_ITEMS_AVERAGE_3D=(
)

#################################################
#
# &PARAM_MKINIT             (init config)
# &PARAM_MKINIT_REAL_ATMOS  (init config)
# &PARAM_MKINIT_REAL_OCEAN  (init config)
# &PARAM_MKINIT_REAL_LAND   (init config)
#
#################################################

INIT_BASENAME="init"

BASENAME_ORG="namelist.grads_boundary.FNL.grib1"
FILETYPE_ORG="GrADS"
PARENT_MP_TYPE=3
USE_FILE_DENSITY=".false."
USE_FILE_LANDWATER=".true."

#################################################
#
# &PARAM_CNVTOPO       (pp config)
# &PARAM_CNVTOPO_*     (pp config)
# &PARAM_COPYTOPO      (pp config)
# &PARAM_CNVLANDUSE    (pp config)
# &PARAM_CNVLANDUSE_*  (pp config)
#
#################################################

TOPODIR="${SCALE_DB}/topo"
LANDUSEDIR="${SCALE_DB}/landuse"

TOPOTYPE=(    "GTOPO30" "GTOPO30" ) # required parameters for each domain
LANDUSETYPE=( "GLCCv2"  "GLCCv2"  ) # required parameters for each domain
COPYTOPO=(    ".false." ".true."  ) # required parameters for each domain

LIMIT_URBAN_FRACTION="0.3"

#################################################
#
# &INFO  (net2g config)
#
#################################################

POPSCA_PLEV=( 850 500 200 )
