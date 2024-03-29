#!/bin/bash

#################################################
#
# &PARAM_COMM_CARTESC_NEST(ALL configs)
# &PARAM_TIME             (ALL configs)
# &PARAM_FILE_HISTORY     (run config)
# &PARAM_ATMOS_BOUNDARY   (run config)
# &PARAM_ATMOS_REFSTATE   (run config)
#
#################################################

NUM_DOMAIN=1 # set number of domains
DOM_NUMBER=2

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
TIME_DT_BOUNDARY="900.0"  # unit: SEC only
TIME_DT_REFSTATE="10800.0"  # unit: SEC only
TIME_DT_HISTORY_2D="150.0"  # unit: SEC only
TIME_DT_HISTORY_3D="150.0"  # unit: SEC only

TIME_DT=(               "15.0" ) # required parameters for each domain - unit: SEC only
TIME_DT_ATMOS_DYN=(      "7.5" ) # required parameters for each domain - unit: SEC only
TIME_DT_ATMOS_PHY_CP=(  "15.0" ) # required parameters for each domain - unit: SEC only
TIME_DT_ATMOS_PHY_MP=(  "15.0" ) # required parameters for each domain - unit: SEC only
TIME_DT_ATMOS_PHY_RD=( "150.0" ) # required parameters for each domain - unit: SEC only
TIME_DT_ATMOS_PHY_SF=(  "15.0" ) # required parameters for each domain - unit: SEC only
TIME_DT_ATMOS_PHY_TB=(  "15.0" ) # required parameters for each domain - unit: SEC only
TIME_DT_ATMOS_PHY_BL=(  "15.0" ) # required parameters for each domain - unit: SEC only
TIME_DT_OCEAN=(         "15.0" ) # required parameters for each domain - unit: SEC only
TIME_DT_LAND=(          "15.0" ) # required parameters for each domain - unit: SEC only
TIME_DT_URBAN=(         "15.0" ) # required parameters for each domain - unit: SEC only

#################################################
#
# &PARAM_PRC                      (ALL configs)
# &PARAM_ATMOS_GRID_CARTESC_INDEX (ALL configs)
# &PARAM_ATMOS_GRID_CARTESC       (ALL configs)
# &PARAM_OCEAN_GRID_CARTESC_INDEX (ALL configs)
# &PARAM_OCEAN_GRID_CARTESC       (ALL configs)
# &PARAM_LAND_GRID_CARTESC_INDEX  (ALL configs)
# &PARAM_LAND_GRID_CARTESC        (ALL configs)
# &PARAM_URBAN_CARTESC_INDEX      (ALL configs)
# &PARAM_URBAN_GRID_CARTESC       (ALL configs)
# &PARAM_MAPPROJECTION            (ALL configs)
#
#################################################

PRC_NUM_X=( 4 ) # required parameters for each domain
PRC_NUM_Y=( 4 ) # required parameters for each domain

KMAX=(   60 ) # required parameters for each domain
IMAXG=( 128 ) # required parameters for each domain
JMAXG=( 128 ) # required parameters for each domain

OKMAX=1
ODZ=( "10.0" ) # required parameters for OKMAX

LKMAX=7
LDZ=( "0.05" "0.15" "0.30" "0.50" "1.00" "2.00" "4.00" ) # required parameters for LKMAX

UKMAX=5
UDZ=( "0.01" "0.01" "0.03" "0.05" "0.10" ) # required parameters for UKMAX

DX=( "7000.0" ) # required parameters for each domain
DY=( "7000.0" ) # required parameters for each domain

DEF_Z=(
"FZ(:) =     40.0000,    80.0000,   120.0000,   160.0000,   200.0000,
            243.2800,   290.1089,   340.7779,   395.6016,   454.9209,
            519.1044,   588.5510,   663.6921,   744.9949,   832.9644,
            928.1475,  1031.1355,  1142.5686,  1263.1392,  1393.5966,
           1534.7515,  1687.4811,  1852.7345,  2031.5387,  2225.0049,
           2434.3352,  2660.8306,  2905.8984,  3171.0618,  3457.9685,
           3768.4016,  4104.2900,  4467.7212,  4860.9536,  5286.4312,
           5746.7979,  6244.9146,  6783.8770,  7367.0342,  7998.0103,
           8680.7266,  9419.4258, 10218.6982, 11083.5107, 12019.2383,
          13019.2383, 14019.2383, 15019.2383, 16019.2383, 17019.2383,
          18019.2383, 19019.2383, 20019.2383, 21019.2383, 22019.2383,
          23019.2383, 24019.2383, 25019.2383, 26019.2383, 27019.2383,"
) # required num. of parameters for each domain

BUFFER_DZ=( "0.0"      ) # required parameters for each domain
BUFFER_DX=( "140000.0" ) # required parameters for each domain
BUFFER_DY=( "140000.0" ) # required parameters for each domain

MAPPROJECTION_BASEPOINT_LON="135.220404"
MAPPROJECTION_BASEPOINT_LAT="34.653396"
MAPPROJECTION_TYPE="LC"
MAPPROJECTION_LC_LAT1="30.0"
MAPPROJECTION_LC_LAT2="40.0"

#################################################
#
# &PARAM_ATMOS   (run config)
# &PARAM_OCEAN   (run config)
# &PARAM_LAND    (run config)
# &PARAM_URBAN   (run config)
#
#################################################

ATMOS_DYN_TYPE=(    "HEVI"        ) # required parameters for each domain
ATMOS_PHY_CP_TYPE=( "KF"          ) # required parameters for each domain
ATMOS_PHY_MP_TYPE=( "TOMITA08"    ) # required parameters for each domain
ATMOS_PHY_RD_TYPE=( "MSTRNX"      ) # required parameters for each domain
ATMOS_PHY_SF_TYPE=( "COUPLE"      ) # required parameters for each domain
ATMOS_PHY_TB_TYPE=( "SMAGORINSKY" ) # required parameters for each domain
ATMOS_PHY_BL_TYPE=( "MYNN"        ) # required parameters for each domain

OCEAN_DYN_TYPE=( "INIT"     ) # required parameters for each domain
LAND_DYN_TYPE=(  "BUCKET"   ) # required parameters for each domain
LAND_SFC_TYPE=(  "SKIN"     ) # required parameters for each domain
URBAN_DYN_TYPE=( "KUSAKA01" ) # required parameters for each domain

#################################################
#
# &HISTORY_ITEM (run config)
#
#################################################

HIST_ITEMS_SNAPSHOT_2D=(
  "MSLP" "U10m" "V10m" "U10" "V10" "T2" "Q2"
  "PREC" "OLR" "SFC_PRES" "SFC_TEMP"
  "OCEAN_SFC_TEMP" "LAND_SFC_TEMP" "URBAN_SFC_TEMP"
  "OCEAN_SFC_ALB_IR_dir"  "OCEAN_SFC_ALB_IR_dif"
  "OCEAN_SFC_ALB_NIR_dir" "OCEAN_SFC_ALB_NIR_dif"
  "OCEAN_SFC_ALB_VIS_dir" "OCEAN_SFC_ALB_VIS_dif"
  "LAND_SFC_ALB_IR_dir"   "LAND_SFC_ALB_IR_dif"
  "LAND_SFC_ALB_NIR_dir"  "LAND_SFC_ALB_NIR_dif"
  "LAND_SFC_ALB_VIS_dir"  "LAND_SFC_ALB_VIS_dif"
  "OCEAN_TEMP" "OCEAN_SFC_Z0M" "LAND_TEMP" "LAND_WATER"
)
HIST_ITEMS_SNAPSHOT_3D=(
  "DENS" "MOMZ" "MOMX" "MOMY" "RHOT"
  "QV" "QC" "QR" "QI" "QS" "QG" "QHYD"
  "PRES" "Umet" "Vmet" "U" "V" "T" "W" "Uabs" "PT" "RH"
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

BASENAME_ORG="history_d01"
FILETYPE_ORG="NetCDF"
PARENT_PRC_NUM_X=2
PARENT_PRC_NUM_Y=2
LATLON_CATALOGUE="latlon_domain_catalogue_d01.txt"
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

TOPOTYPE=(     "GTOPO30" ) # required parameters for each domain
LANDUSETYPE=(  "GLCCv2"  ) # required parameters for each domain
COPYTOPO=(     ".false." ) # required parameters for each domain
SMOOTH_LOCAL=( ".true."  ) # required parameters for each domain

SMOOTH_ITELIM="10000"
MAXSLOPE_RATIO="5.0"
LIMIT_URBAN_FRACTION="0.3"

#################################################
#
# &INFO  (net2g config)
#
#################################################

POPSCA_PLEV=( 850 500 200 )
