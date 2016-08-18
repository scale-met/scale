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

NUM_DOMAIN=3 # set number of domains

RUN_DATE_YEAR=2011
RUN_DATE_MON=9
RUN_DATE_DAY=18
RUN_DATE_HOUR=0
RUN_DATE_MIN=0
RUN_DATE_SEC=0
RUN_DATE_MSEC=0

TIME_DURATION="432000.0"
TIME_UNIT="SEC" # SEC only

TIME_DT_BOUNDARY="21600.0"
TIME_DT_REFSTATE="10800.0"
TIME_DT_HISTORY_2D="3600.0"
TIME_DT_HISTORY_3D="21600.0"

TIME_DT=(               "60.0"  "15.0"   "3.0" ) # required parameters for each domain
TIME_DT_ATMOS_DYN=(     "30.0"   "7.5"   "1.5" ) # required parameters for each domain
TIME_DT_ATMOS_PHY_MP=(  "60.0"  "15.0"   "3.0" ) # required parameters for each domain
TIME_DT_ATMOS_PHY_RD=( "300.0" "300.0" "300.0" ) # required parameters for each domain
TIME_DT_ATMOS_PHY_SF=(  "60.0"  "60.0"  "60.0" ) # required parameters for each domain
TIME_DT_ATMOS_PHY_TB=(  "60.0"  "15.0"   "3.0" ) # required parameters for each domain
TIME_DT_OCEAN=(         "60.0"  "60.0"  "60.0" ) # required parameters for each domain
TIME_DT_LAND=(          "60.0"  "60.0"  "60.0" ) # required parameters for each domain
TIME_DT_URBAN=(         "60.0"  "60.0"  "60.0" ) # required parameters for each domain

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

PRC_NUM_X=(  6  8 16 ) # required parameters for each domain
PRC_NUM_Y=(  6 12 18 ) # required parameters for each domain
 
KMAX=( 48 60 80 ) # required parameters for each domain
IMAX=( 40 36 27 ) # required parameters for each domain
JMAX=( 40 24 24 ) # required parameters for each domain

LKMAX=5
LDZ=( "0.05" "0.15" "0.30" "0.50" "1.00" ) # required parameters for LKMAX

UKMAX=5
UDZ=( "0.01" "0.01" "0.03" "0.05" "0.10" ) # required parameters for UKMAX

DX=( "15000.0" "5000.0" "1000.0" ) # required parameters for each domain
DY=( "15000.0" "5000.0" "1000.0" ) # required parameters for each domain

DEF_Z=( 
"FZ(:) =    160.0000,   320.0000,   480.0000,   640.0000,   800.0000,
            971.2000,  1154.3840,  1350.3910,  1560.1184,  1784.5267,
           2024.6437,  2281.5688,  2556.4788,  2850.6323,  3165.3767,
           3502.1533,  3862.5044,  4248.0801,  4660.6460,  5102.0913,
           5574.4380,  6079.8491,  6620.6392,  7199.2847,  7818.4355,
           8480.9268,  9189.7920,  9948.2773, 10759.8564, 11628.2461,
          12528.2461, 13428.2461, 14328.2461, 15228.2461, 16128.2461,
          17028.2461, 17928.2461, 18828.2461, 19728.2461, 20628.2461,
          21580.7617, 22588.8574, 23655.7754, 24784.9473, 25980.0059,
          27244.7969, 28583.3887, 30000.0879,"

"FZ(:) =    110.0000,   220.0000,   330.0000,   440.0000,   550.0000,
            666.1600,   788.8249,   918.3591,  1055.1472,  1199.5955,
           1352.1328,  1513.2123,  1683.3123,  1862.9379,  2052.6226,
           2252.9297,  2464.4541,  2687.8240,  2923.7026,  3172.7905,
           3435.8274,  3713.5942,  4006.9160,  4316.6636,  4643.7568,
           4989.1675,  5353.9209,  5739.1006,  6145.8501,  6575.3774,
           7028.9585,  7507.9399,  8013.7441,  8547.8730,  9111.9131,
           9707.5391, 10336.5205, 11000.7246, 11702.1240, 12442.8018,
          13242.8018, 14042.8018, 14842.8018, 15642.8018, 16442.8008,
          17242.8008, 18042.8008, 18842.8008, 19642.8008, 20442.8008,
          21242.8008, 22042.8008, 22842.8008, 23642.8008, 24442.8008,
          25242.8008, 26042.8008, 26842.8008, 27642.8008, 28442.8008,"

"FZ(:) =     60.0000,   120.0000,   180.0000,   240.0000,   300.0000,
            362.8800,   428.7783,   497.8396,   570.2159,   646.0663,
            725.5576,   808.8643,   896.1698,   987.6660,  1083.5540,
           1184.0446,  1289.3588,  1399.7280,  1515.3950,  1636.6140,
           1763.6515,  1896.7867,  2036.3125,  2182.5354,  2335.7771,
           2496.3745,  2664.6807,  2841.0654,  3025.9167,  3219.6409,
           3422.6638,  3635.4319,  3858.4128,  4092.0969,  4336.9980,
           4593.6543,  4862.6299,  5144.5161,  5439.9331,  5749.5303,
           6073.9883,  6414.0205,  6770.3745,  7143.8335,  7535.2188,
           7945.3906,  8375.2510,  8825.7451,  9297.8633,  9792.6436,
          10292.6436, 10792.6436, 11292.6436, 11792.6436, 12292.6436,
          12792.6436, 13292.6436, 13792.6436, 14292.6436, 14792.6436,
          15292.6436, 15792.6436, 16292.6436, 16792.6445, 17292.6445,
          17792.6445, 18292.6445, 18792.6445, 19292.6445, 19792.6445,
          20292.6445, 20792.6445, 21292.6445, 21792.6445, 22292.6445,
          22792.6445, 23292.6445, 23792.6445, 24292.6445, 24792.6445,"
) # required num. of parameters for each domain

BUFFER_DZ=( "5000.0"   "5000.0"   "5000.0"  ) # required parameters for each domain
BUFFER_DX=( "300000.0" "100000.0" "20000.0" ) # required parameters for each domain
BUFFER_DY=( "300000.0" "100000.0" "20000.0" ) # required parameters for each domain

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

ATMOS_DYN_TYPE=(    "HEVI"     "HEVI"     "HEVI"     ) # required parameters for each domain
ATMOS_PHY_MP_TYPE=( "TOMITA08" "TOMITA08" "TOMITA08" ) # required parameters for each domain
ATMOS_PHY_RD_TYPE=( "MSTRNX"   "MSTRNX"   "MSTRNX"   ) # required parameters for each domain
ATMOS_PHY_SF_TYPE=( "COUPLE"   "COUPLE"   "COUPLE"   ) # required parameters for each domain
ATMOS_PHY_TB_TYPE=( "MYNN"     "MYNN"     "MYNN"     ) # required parameters for each domain

OCEAN_TYPE=( "CONST" "CONST" "CONST" ) # required parameters for each domain
LAND_TYPE=(  "SLAB"  "SLAB"  "SLAB"  ) # required parameters for each domain
URBAN_TYPE=( "SLC"   "SLC"   "SLC"   ) # required parameters for each domain

#################################################
#
# &HISTITEM (run config)
#
#################################################

HIST_ITEMS_SNAPSHOT_2D=(
  "MSLP" "SFC_PRES" "U10" "V10" "T2" "Q2"
  "LWP" "IWP" "PW" "OLR" "LHFLX" "SHFLX" "GHFLX"
)
HIST_ITEMS_SNAPSHOT_3D=(
  "DENS" "QV" "QC" "QR" "QI" "QS" "QG"
  "T" "PRES" "U" "V" "W" "RH"
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
BASENAME_ORG="namelist.grads_boundary"
FILETYPE_ORG="GrADS"
PARENT_MP_TYPE=3
USE_FILE_DENSITY=".false."
USE_FILE_LANDWATER=".true."

#################################################
#
# &PARAM_CNVTOPO     (pp config)
# &PARAM_COPYTOPO    (pp config)
# &PARAM_CNVLANDUSE  (pp config)
#
#################################################

TOPOTYPE=(    "GTOPO30" "GTOPO30" "DEM50M" ) # required parameters for each domain
LANDUSETYPE=( "GLCCv2"  "GLCCv2"  "LU100M" ) # required parameters for each domain

MAXSLOPE=( "0.6"     "1.2"     "3.4"    ) # required parameters for each domain
COPYTOPO=( ".false." ".true."  ".true." ) # required parameters for each domain

LIMIT_URBAN_FRACTION="0.3"

#################################################
#
# &INFO  (net2g config)
#
#################################################

POPSCA_PLEV=( 850 500 200 )
