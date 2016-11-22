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

NUM_DOMAIN=1 # set number of domains

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

TIME_DT=(               "90.0" ) # required parameters for each domain - unit: SEC only
TIME_DT_ATMOS_DYN=(     "45.0" ) # required parameters for each domain - unit: SEC only
TIME_DT_ATMOS_PHY_MP=(  "90.0" ) # required parameters for each domain - unit: SEC only
TIME_DT_ATMOS_PHY_RD=( "900.0" ) # required parameters for each domain - unit: SEC only
TIME_DT_ATMOS_PHY_SF=(  "90.0" ) # required parameters for each domain - unit: SEC only
TIME_DT_ATMOS_PHY_TB=(  "90.0" ) # required parameters for each domain - unit: SEC only
TIME_DT_OCEAN=(        "450.0" ) # required parameters for each domain - unit: SEC only
TIME_DT_LAND=(         "450.0" ) # required parameters for each domain - unit: SEC only
TIME_DT_URBAN=(        "450.0" ) # required parameters for each domain - unit: SEC only

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

PRC_NUM_X=( 2 ) # required parameters for each domain
PRC_NUM_Y=( 2 ) # required parameters for each domain
 
KMAX=( 36 ) # required parameters for each domain
IMAX=( 45 ) # required parameters for each domain
JMAX=( 45 ) # required parameters for each domain

LKMAX=5
LDZ=( "0.05" "0.15" "0.30" "0.50" "1.00" ) # required parameters for LKMAX

UKMAX=5
UDZ=( "0.01" "0.01" "0.03" "0.05" "0.10" ) # required parameters for UKMAX

DX=( "20000.0" ) # required parameters for each domain
DY=( "20000.0" ) # required parameters for each domain

DEF_Z=( 
"FZ(:) =     80.8410,   248.8210,   429.8820,   625.0450,   835.4090,
           1062.1580,  1306.5650,  1570.0080,  1853.9690,  2160.0470,
           2489.9630,  2845.5750,  3228.8830,  3642.0440,  4087.3840,
           4567.4090,  5084.8200,  5642.5300,  6243.6760,  6891.6420,
           7590.0740,  8342.9040,  9154.3670, 10029.0280, 10971.8150,
          11988.0300, 13083.3900, 14264.0600, 15536.6850, 16908.4300,
          18387.0100, 19980.7500, 21698.6150, 23550.2750, 25546.1550,
          28113.2050,"
) # required num. of parameters for each domain

BUFFER_DZ=( "5000.0"   ) # required parameters for each domain
BUFFER_DX=( "400000.0" ) # required parameters for each domain
BUFFER_DY=( "400000.0" ) # required parameters for each domain

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

ATMOS_DYN_TYPE=(    "HEVI"     ) # required parameters for each domain
ATMOS_PHY_MP_TYPE=( "TOMITA08" ) # required parameters for each domain
ATMOS_PHY_RD_TYPE=( "MSTRNX"   ) # required parameters for each domain
ATMOS_PHY_SF_TYPE=( "COUPLE"   ) # required parameters for each domain
ATMOS_PHY_TB_TYPE=( "HYBRID"   ) # required parameters for each domain

OCEAN_TYPE=( "CONST" ) # required parameters for each domain
LAND_TYPE=(  "SLAB"  ) # required parameters for each domain
URBAN_TYPE=( "SLC"   ) # required parameters for each domain

#################################################
#
# &HISTITEM (run config)
#
#################################################

HIST_ITEMS_SNAPSHOT_2D=(
  "MSLP" "PREC" "OLR" "U10" "V10" "T2" "Q2" "SFC_PRES" "SFC_TEMP"
)
HIST_ITEMS_SNAPSHOT_3D=(
  "DENS" "QV" "QHYD" "PRES" "U" "V" "T" "W" "Uabs" "PT" "RH"
)
HIST_ITEMS_AVERAGE_2D=(
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

TOPOTYPE=(    "GTOPO30" ) # required parameters for each domain
LANDUSETYPE=( "GLCCv2"  ) # required parameters for each domain
COPYTOPO=(    ".false." ) # required parameters for each domain

MAXSLOPE_RATIO="1.0"
LIMIT_URBAN_FRACTION="0.3"

#################################################
#
# &INFO  (net2g config)
#
#################################################

POPSCA_PLEV=( 850 500 200 )
