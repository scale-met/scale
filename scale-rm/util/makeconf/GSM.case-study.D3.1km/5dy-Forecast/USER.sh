#!/bin/bash

#################################################
#
# set user defined parameters
#
#################################################

# set number of domains
NUM_DOMAIN=3

# set date of run
YEAR=2011
MON=9
DAY=18
HOUR=0
MIN=0
SEC=0
MSEC=0

# time parameters
TIME_DURATION="432000.D0"
TIME_UNIT="SEC" # SEC only

TIME_DT_BOUNDARY="21600.D0"
TIME_DT_REFSTATE="10800.D0"
TIME_DT_HISTORY_2D="3600.D0"
TIME_DT_HISTORY_3D="21600.D0"

TIME_DT=(               "60.00D0"  "15.00D0"   "3.00D0" ) # required num. of parameters for each domain
TIME_DT_ATMOS_DYN=(     "30.00D0"   "7.50D0"   "1.50D0" ) # required num. of parameters for each domain
TIME_DT_ATMOS_PHY_MP=(  "60.00D0"  "15.00D0"   "3.00D0" ) # required num. of parameters for each domain
TIME_DT_ATMOS_PHY_RD=( "300.00D0" "300.00D0" "300.00D0" ) # required num. of parameters for each domain
TIME_DT_ATMOS_PHY_SF=(  "60.00D0"  "60.00D0"  "60.00D0" ) # required num. of parameters for each domain
TIME_DT_ATMOS_PHY_TB=(  "60.00D0"  "15.00D0"   "3.00D0" ) # required num. of parameters for each domain
TIME_DT_OCEAN=(         "60.00D0"  "60.00D0"  "60.00D0" ) # required num. of parameters for each domain
TIME_DT_LAND=(          "60.00D0"  "60.00D0"  "60.00D0" ) # required num. of parameters for each domain
TIME_DT_URBAN=(         "60.00D0"  "60.00D0"  "60.00D0" ) # required num. of parameters for each domain

# region parameters
PRC_NUM_X=(  6  8 16 ) # required num. of parameters for each domain
PRC_NUM_Y=(  6 12 18 ) # required num. of parameters for each domain
 
KMAX=( 48 60 80 ) # required num. of parameters for each domain
IMAX=( 40 36 27 ) # required num. of parameters for each domain
JMAX=( 40 24 24 ) # required num. of parameters for each domain

LKMAX=5
LDZ=( "0.05D0" "0.15D0" "0.30D0" "0.50D0" "1.00D0" ) # required num. of parameters for LKMAX

UKMAX=5
UDZ=( "0.01D0" "0.01D0" "0.03D0" "0.05D0" "0.10D0" ) # required num. of parameters for UKMAX

DX=( "15000.0D0" "5000.0D0" "1000.0D0" ) # required num. of parameters for each domain
DY=( "15000.0D0" "5000.0D0" "1000.0D0" ) # required num. of parameters for each domain

DEF_Z=( 
"FZ(:) =    160.0000D0,   320.0000D0,   480.0000D0,   640.0000D0,   800.0000D0,
            971.2000D0,  1154.3840D0,  1350.3910D0,  1560.1184D0,  1784.5267D0,
           2024.6437D0,  2281.5688D0,  2556.4788D0,  2850.6323D0,  3165.3767D0,
           3502.1533D0,  3862.5044D0,  4248.0801D0,  4660.6460D0,  5102.0913D0,
           5574.4380D0,  6079.8491D0,  6620.6392D0,  7199.2847D0,  7818.4355D0,
           8480.9268D0,  9189.7920D0,  9948.2773D0, 10759.8564D0, 11628.2461D0,
          12528.2461D0, 13428.2461D0, 14328.2461D0, 15228.2461D0, 16128.2461D0,
          17028.2461D0, 17928.2461D0, 18828.2461D0, 19728.2461D0, 20628.2461D0,
          21580.7617D0, 22588.8574D0, 23655.7754D0, 24784.9473D0, 25980.0059D0,
          27244.7969D0, 28583.3887D0, 30000.0879D0,"

"FZ(:) =    110.0000D0,   220.0000D0,   330.0000D0,   440.0000D0,   550.0000D0,
            666.1600D0,   788.8249D0,   918.3591D0,  1055.1472D0,  1199.5955D0,
           1352.1328D0,  1513.2123D0,  1683.3123D0,  1862.9379D0,  2052.6226D0,
           2252.9297D0,  2464.4541D0,  2687.8240D0,  2923.7026D0,  3172.7905D0,
           3435.8274D0,  3713.5942D0,  4006.9160D0,  4316.6636D0,  4643.7568D0,
           4989.1675D0,  5353.9209D0,  5739.1006D0,  6145.8501D0,  6575.3774D0,
           7028.9585D0,  7507.9399D0,  8013.7441D0,  8547.8730D0,  9111.9131D0,
           9707.5391D0, 10336.5205D0, 11000.7246D0, 11702.1240D0, 12442.8018D0,
          13242.8018D0, 14042.8018D0, 14842.8018D0, 15642.8018D0, 16442.8008D0,
          17242.8008D0, 18042.8008D0, 18842.8008D0, 19642.8008D0, 20442.8008D0,
          21242.8008D0, 22042.8008D0, 22842.8008D0, 23642.8008D0, 24442.8008D0,
          25242.8008D0, 26042.8008D0, 26842.8008D0, 27642.8008D0, 28442.8008D0,"

"FZ(:) =     60.0000D0,   120.0000D0,   180.0000D0,   240.0000D0,   300.0000D0,
            362.8800D0,   428.7783D0,   497.8396D0,   570.2159D0,   646.0663D0,
            725.5576D0,   808.8643D0,   896.1698D0,   987.6660D0,  1083.5540D0,
           1184.0446D0,  1289.3588D0,  1399.7280D0,  1515.3950D0,  1636.6140D0,
           1763.6515D0,  1896.7867D0,  2036.3125D0,  2182.5354D0,  2335.7771D0,
           2496.3745D0,  2664.6807D0,  2841.0654D0,  3025.9167D0,  3219.6409D0,
           3422.6638D0,  3635.4319D0,  3858.4128D0,  4092.0969D0,  4336.9980D0,
           4593.6543D0,  4862.6299D0,  5144.5161D0,  5439.9331D0,  5749.5303D0,
           6073.9883D0,  6414.0205D0,  6770.3745D0,  7143.8335D0,  7535.2188D0,
           7945.3906D0,  8375.2510D0,  8825.7451D0,  9297.8633D0,  9792.6436D0,
          10292.6436D0, 10792.6436D0, 11292.6436D0, 11792.6436D0, 12292.6436D0,
          12792.6436D0, 13292.6436D0, 13792.6436D0, 14292.6436D0, 14792.6436D0,
          15292.6436D0, 15792.6436D0, 16292.6436D0, 16792.6445D0, 17292.6445D0,
          17792.6445D0, 18292.6445D0, 18792.6445D0, 19292.6445D0, 19792.6445D0,
          20292.6445D0, 20792.6445D0, 21292.6445D0, 21792.6445D0, 22292.6445D0,
          22792.6445D0, 23292.6445D0, 23792.6445D0, 24292.6445D0, 24792.6445D0,"
) # required num. of parameters for each domain

BUFFER_DZ=( "5000.0D0"   "5000.0D0"   "5000.0D0"  ) # required num. of parameters for each domain
BUFFER_DX=( "300000.0D0" "100000.0D0" "20000.0D0" ) # required num. of parameters for each domain
BUFFER_DY=( "300000.0D0" "100000.0D0" "20000.0D0" ) # required num. of parameters for each domain

MPRJ_BASEPOINT_LON="135.220404D0"
MPRJ_BASEPOINT_LAT="34.653396D0"
MPRJ_TYPE="LC"
MPRJ_LC_LAT1="30.0D0"
MPRJ_LC_LAT2="40.0D0"

# type parameters
CONST_THERMODYN_TYPE="SIMPLE"

ATMOS_DYN_TYPE=(    "HEVI"     "HEVI"     "HEVI"     ) # required num. of parameters for each domain
ATMOS_PHY_MP_TYPE=( "TOMITA08" "TOMITA08" "TOMITA08" ) # required num. of parameters for each domain
ATMOS_PHY_RD_TYPE=( "MSTRNX"   "MSTRNX"   "MSTRNX"   ) # required num. of parameters for each domain
ATMOS_PHY_SF_TYPE=( "COUPLE"   "COUPLE"   "COUPLE"   ) # required num. of parameters for each domain
ATMOS_PHY_TB_TYPE=( "MYNN"     "MYNN"     "MYNN"     ) # required num. of parameters for each domain

ATMOS_DYN_TINTEG_LARGE_TYPE="EULER"               # EULER / RK3
ATMOS_DYN_TINTEG_SHORT_TYPE="RK4"                 # RK3 / RK4
ATMOS_DYN_TINTEG_TRACER_TYPE="RK3WS2002"          # EULER / RK3WS2002
ATMOS_DYN_FVM_FLUX_TYPE="CD4"                     # UD3KOREN1993 / UD3 / CD4
ATMOS_DYN_FVM_FLUX_TRACER_TYPE="UD3KOREN1993"     # UD3KOREN1993 / UD3 / CD4
ATMOS_DYN_NUMERICAL_DIFF_COEF="1.0D-2"
ATMOS_DYN_NUMERICAL_DIFF_COEF_TRACER="0.0D0"
ATMOS_DYN_FLAG_FCT_TRACER=".false."

OCEAN_TYPE=( "CONST" "CONST" "CONST" ) # required num. of parameters for each domain
LAND_TYPE=(  "SLAB"  "SLAB"  "SLAB"  ) # required num. of parameters for each domain
URBAN_TYPE=( "SLC"   "SLC"   "SLC"   ) # required num. of parameters for each domain

ATMOS_BOUNDARY_TAUX=( "600.0D0" "150.0D0" "30.0D0" ) # required num. of parameters for each domain
ATMOS_BOUNDARY_TAUY=( "600.0D0" "150.0D0" "30.0D0" ) # required num. of parameters for each domain

ATMOS_BOUNDARY_ALPHAFACT_DENS=( "1.0D0" "1.0D0" "1.0D0" ) # required num. of parameters for each domain

ATMOS_PHY_TB_SMG_consistent_tke=( ".false." ".false." ".false." ) # required num. of parameters for each domain
ATMOS_PHY_TB_SMG_implicit=(       ".true."  ".true."  ".true."  ) # required num. of parameters for each domain

# history parameters
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

##### for init
INIT_BASENAME="init"
BASENAME_ORG="namelist.grads_boundary"
FILETYPE_ORG="GrADS"
PARENT_MP_TYPE=3
USE_FILE_DENSITY=".false."
USE_FILE_LANDWATER=".true."

##### for pp
TOPOTYPE=(    "GTOPO30" "GTOPO30" "DEM50M" ) # required num. of parameters for each domain: GTOPO30, DEM50M
LANDUSETYPE=( "GLCCv2"  "GLCCv2"  "LU100M" ) # required num. of parameters for each domain: GLCCv2, LU100M

MAXSLOPE=( "0.6D0"   "1.2D0"   "3.4D0"  ) # required num. of parameters for each domain
COPYTOPO=( ".false." ".true."  ".true." ) # required num. of parameters for each domain

LIMIT_URBAN_FRACTION="0.3D0"

##### for net2g
POPSCA_PLEV=( 850 500 200 )
