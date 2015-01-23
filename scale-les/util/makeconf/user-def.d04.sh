#!/bin/bash
#
#################################################
#
# set user defined parameters
#
#################################################

NUM_DOMAIN=4

# use nesting or not
RUN_USE_NESTING=".true."

# time parameters
DEF_YEAR=${1:-0000}   # Year
DEF_MON=${2:-1}       # Mon
DEF_DAY=${3:-1}       # Day
DEF_HOUR=${4:-0}      # Hour
DEF_MIN=${5:-0}       # Min
DEF_SEC=${6:-0}       # Sec
DEF_MSEC=${7:-0}      # milli-Sec

STARTDATE=( ${DEF_YEAR} ${DEF_MON} ${DEF_DAY} ${DEF_HOUR} ${DEF_MIN} ${DEF_SEC} ${DEF_MSEC} )

TIME_DURATION="410400.D0"
TIME_UNIT="SEC"

TIME_DT_RESTART="86400.D0"
TIME_DT_HISTORY="3600.D0"
TIME_DT_HISTORY_HI="600.D0"
TIME_DT_BOUNDARY="21600.D0"
TIME_DT_REFSTATE="10800.D0"

TIME_DT=(               "30.0D0"  "10.0D0"   "2.0D0"   "0.4D0" ) # required num. of parameters for each domain
TIME_DT_ATMOS_DYN=(      "7.5D0"   "2.5D0"   "0.5D0"   "0.1D0" ) # required num. of parameters for each domain
TIME_DT_ATMOS_PHY_MP=(  "30.0D0"  "10.0D0"   "2.0D0"   "0.4D0" ) # required num. of parameters for each domain
TIME_DT_ATMOS_PHY_RD=( "300.0D0" "300.0D0" "300.0D0" "300.0D0" ) # required num. of parameters for each domain
TIME_DT_ATMOS_PHY_SF=(  "30.0D0"  "10.0D0"   "2.0D0"   "0.4D0" ) # required num. of parameters for each domain
TIME_DT_ATMOS_PHY_TB=(  "30.0D0"  "10.0D0"   "2.0D0"   "0.4D0" ) # required num. of parameters for each domain
TIME_DT_OCEAN=(        "300.0D0" "300.0D0" "300.0D0" "300.0D0" ) # required num. of parameters for each domain
TIME_DT_LAND=(         "300.0D0" "300.0D0" "300.0D0" "300.0D0" ) # required num. of parameters for each domain
TIME_DT_URBAN=(        "300.0D0" "300.0D0" "300.0D0" "300.0D0" ) # required num. of parameters for each domain

# region parameters
PRC_NUM_X=(  6  9 16 40 ) # required num. of parameters for each domain
PRC_NUM_Y=(  6 12 18 54 ) # required num. of parameters for each domain
 
KMAX=( 36 60 80 80 ) # required num. of parameters for each domain
IMAX=( 56 48 32 24 ) # required num. of parameters for each domain
JMAX=( 56 32 24 16 ) # required num. of parameters for each domain

LKMAX=5
LDZ=( "0.05D0" "0.15D0" "0.30D0" "0.50D0" "1.00D0" ) # required num. of parameters for LKMAX

UKMAX=5
UDZ=( "0.01D0" "0.01D0" "0.03D0" "0.05D0" "0.10D0" ) # required num. of parameters for UKMAX

DX=( "7500.0D0" "2500.0D0" "500.0D0" "100.0D0" ) # required num. of parameters for each domain
DY=( "7500.0D0" "2500.0D0" "500.0D0" "100.0D0" ) # required num. of parameters for each domain

DEF_Z=( 
"FZ(:) =    161.6830D0,   335.9580D0,   523.8060D0,   726.2850D0,   944.5340D0,
           1179.7810D0,  1433.3490D0,  1706.6670D0,  2001.2720D0,  2318.8220D0,
           2661.1040D0,  3030.0450D0,  3427.7200D0,  3856.3680D0,  4318.4000D0,
           4816.4180D0,  5353.2230D0,  5931.8370D0,  6555.5160D0,  7227.7690D0,
           7952.3800D0,  8733.4280D0,  9575.3060D0, 10482.7500D0, 11460.8800D0,
          12515.1800D0, 13651.6000D0, 14876.5200D0, 16196.8500D0, 17620.0100D0,
          19154.0100D0, 20807.4900D0, 22589.7400D0, 24510.8100D0, 26581.5000D0,
          29644.9100D0,"

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

BUFFER_DZ=(  "5000.0D0"  "5000.0D0" "5000.0D0" "1000.0D0" ) # required num. of parameters for each domain
BUFFER_DX=( "75000.0D0" "25000.0D0" "5000.0D0" "1000.0D0" ) # required num. of parameters for each domain
BUFFER_DY=( "75000.0D0" "25000.0D0" "5000.0D0" "1000.0D0" ) # required num. of parameters for each domain

MPRJ_BASEPOINT_LON="135.220404D0"
MPRJ_BASEPOINT_LAT="34.653396D0"
MPRJ_TYPE="LC"
MPRJ_LC_LAT1="30.0D0"
MPRJ_LC_LAT2="40.0D0"

# type parameters
CONST_THERMODYN_TYPE="SIMPLE"

ATMOS_DYN_TYPE=(    "HEVI"     "HEVI"     "HEVI"     "HEVI"        ) # required num. of parameters for each domain
ATMOS_PHY_MP_TYPE=( "TOMITA08" "TOMITA08" "TOMITA08" "TOMITA08"    ) # required num. of parameters for each domain
ATMOS_PHY_RD_TYPE=( "MSTRNX"   "MSTRNX"   "MSTRNX"   "MSTRNX"      ) # required num. of parameters for each domain
ATMOS_PHY_SF_TYPE=( "COUPLE"   "COUPLE"   "COUPLE"   "COUPLE"      ) # required num. of parameters for each domain
ATMOS_PHY_TB_TYPE=( "MYNN"     "MYNN"     "MYNN"     "SMAGORINSKY" ) # required num. of parameters for each domain

OCEAN_TYPE=( "CONST" "CONST" "CONST" "CONST" ) # required num. of parameters for each domain
LAND_TYPE=(  "SLAB"  "SLAB"  "SLAB"  "SLAB"  ) # required num. of parameters for each domain
URBAN_TYPE=( "SLC"   "SLC"   "SLC"   "SLC"   ) # required num. of parameters for each domain

# history parameters
HIST_ITEMS_SNAPSHOT=(
# ATMOS history
  "DENS" "MOMZ" "MOMX" "MOMY" "RHOT"
  "QV" "QC" "QR" "QI" "QS" "QG" "QHYD" "QLIQ" "QICE"
  "T" "PRES" "U" "V" "W" "Uabs" "PT" "RH" "RHL" "RHI"
  "SFC_PRES" "SFC_TEMP" "SFC_Z0M"
  "U10" "V10" "Uabs10" "T2" "Q2" "MSLP" "LHFLX" "SHFLX" "GHFLX"
  "SFLX_LW_up" "SFLX_LW_dn" "SFLX_SW_up" "SFLX_SW_dn"
  "TOAFLX_LW_up" "TOAFLX_LW_dn" "TOAFLX_SW_up" "TOAFLX_SW_dn"
  "OSR" "OLR" "SLR" "SSR" "RADFLUX_SWUP" "RADFLUX_SWDN"
# OCEAN history
  "OCEAN_TEMP" "OCEAN_SFC_TEMP" "OCEAN_ALB_SW" "OCEAN_ALB_LW" "OCEAN_SFC_Z0M" "OCEAN_SFC_Z0H" "OCEAN_SFC_Z0E"
# LAND history
  "LAND_TEMP" "LAND_WATER" "LAND_SFC_TEMP" "LAND_ALB_SW" "LAND_ALB_LW"
# URBAN history
  "URBAN_TC" "URBAN_SFC_TEMP"
)
HIST_ITEMS_SNAPSHOT_HI=(
)
HIST_ITEMS_AVERAGE=(
  "PREC"
)
HIST_ITEMS_AVERAGE_HI=(
  "RAIN" "SNOW"
)

##### for init
BASENAME_ORG=""
FILETYPE_ORG="NICAM-NETCDF"
NUMBER_OF_FILES=6
WRF_FILE_TYPE=""

INIT_USE_NESTING=".false."
PARENT_PRC_NUM_X=""        # used if INIT_USE_NESTING = .true.
PARENT_PRC_NUM_Y=""        # used if INIT_USE_NESTING = .true.
PARENT_KMAX=""             # used if INIT_USE_NESTING = .true.
PARENT_IMAX=""             # used if INIT_USE_NESTING = .true.
PARENT_JMAX=""             # used if INIT_USE_NESTING = .true.
PARENT_LKMAX=""            # used if INIT_USE_NESTING = .true.

# linker
LINK_BASEDIR="/data/share005/data/nicam/amip/trimmed_data"
LINK_OUTPUTDIR="input"
LINK_BOUNDARY_DT=86400
LINK_FTYPE='peall'

##### for pp
TOPOTYPE=(    "GTOPO30" "GTOPO30" "DEM50M" "DEM50M" ) # required num. of parameters for each domain: GTOPO30, DEM50M
LANDUSETYPE=( "GLCCv2"  "GLCCv2"  "LU100M" "LU100M" ) # required num. of parameters for each domain: GLCCv2, LU100M

MAXSLOPE=(     "0.5D0" "1.5D0" "7.0D0" "14.0D0" ) # required num. of parameters for each domain
MAXSLOPE_BND=( "0.5D0" "1.5D0" "7.0D0" "14.0D0" ) # required num. of parameters for each domain

LIMIT_URBAN_FRACTION="0.3D0"

##### for staging
TRIMDIR=true

INDIR="input"   # used if TRIMDIR = false
OUTDIR="output" # used if TRIMDIR = false
