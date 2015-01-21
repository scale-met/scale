#!/bin/bash
#
# SCALE-LES configuration file generator
#
#################################################
#
# set user defined parameters
#
#################################################

NUM_DOMAIN=3

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

TIME_DT=(               "30.0D0"  "10.0D0"   "2.0D0" ) # required num. of parameters for each domain
TIME_DT_ATMOS_DYN=(      "7.5D0"   "2.5D0"   "0.5D0" ) # required num. of parameters for each domain
TIME_DT_ATMOS_PHY_MP=(  "30.0D0"  "10.0D0"   "2.0D0" ) # required num. of parameters for each domain
TIME_DT_ATMOS_PHY_RD=( "300.0D0" "300.0D0" "300.0D0" ) # required num. of parameters for each domain
TIME_DT_ATMOS_PHY_SF=(  "30.0D0"  "10.0D0"   "2.0D0" ) # required num. of parameters for each domain
TIME_DT_ATMOS_PHY_TB=(  "30.0D0"  "10.0D0"   "2.0D0" ) # required num. of parameters for each domain
TIME_DT_OCEAN=(        "300.0D0" "300.0D0" "300.0D0" ) # required num. of parameters for each domain
TIME_DT_LAND=(         "300.0D0" "300.0D0" "300.0D0" ) # required num. of parameters for each domain
TIME_DT_URBAN=(        "300.0D0" "300.0D0" "300.0D0" ) # required num. of parameters for each domain

# region parameters
PRC_NUM_X=(  6  9 16 ) # required num. of parameters for each domain
PRC_NUM_Y=(  6 12 18 ) # required num. of parameters for each domain
 
KMAX=( 36 60 80 ) # required num. of parameters for each domain
IMAX=( 56 48 32 ) # required num. of parameters for each domain
JMAX=( 56 32 24 ) # required num. of parameters for each domain

LKMAX=5
LDZ=( "0.05D0" "0.15D0" "0.30D0" "0.50D0" "1.00D0" ) # required num. of parameters for LKMAX

UKMAX=5
UDZ=( "0.01D0" "0.01D0" "0.03D0" "0.05D0" "0.10D0" ) # required num. of parameters for UKMAX

DX=( "7500.0D0" "2500.0D0" "500.0D0" ) # required num. of parameters for each domain
DY=( "7500.0D0" "2500.0D0" "500.0D0" ) # required num. of parameters for each domain

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
) # required num. of parameters for each domain

BUFFER_DZ=(  "5000.0D0"  "5000.0D0" "5000.0D0" ) # required num. of parameters for each domain
BUFFER_DX=( "75000.0D0" "25000.0D0" "5000.0D0" ) # required num. of parameters for each domain
BUFFER_DY=( "75000.0D0" "25000.0D0" "5000.0D0" ) # required num. of parameters for each domain

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

OCEAN_TYPE=( "CONST" "CONST" "CONST" ) # required num. of parameters for each domain
LAND_TYPE=(  "SLAB"  "SLAB"  "SLAB"  ) # required num. of parameters for each domain
URBAN_TYPE=( "SLC"   "SLC"   "SLC"   ) # required num. of parameters for each domain

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
)
HIST_ITEMS_AVERAGE_HI=(
  "RAIN" "SNOW" "PREC"
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
TOPOTYPE=(    "GTOPO30" "GTOPO30" "DEM50M" ) # required num. of parameters for each domain: GTOPO30, DEM50M
LANDUSETYPE=( "GLCCv2"  "GLCCv2"  "LU100M" ) # required num. of parameters for each domain: GLCCv2, LU100M

MAXSLOPE=(     "0.5D0" "1.5D0" "7.0D0" ) # required num. of parameters for each domain
MAXSLOPE_BND=( "0.5D0" "1.5D0" "7.0D0" ) # required num. of parameters for each domain

##### for staging
TRIMDIR=true
#TRIMDIR=false

INDIR="input"   # used if TRIMDIR = false
OUTDIR="output" # used if TRIMDIR = false

#################################################
#
# check parameters
#
#################################################

if [ 7 -ne ${#STARTDATE[*]} ]; then echo "Error: Wrong array size (STARTDATE)."; exit 1; fi

if [ ${NUM_DOMAIN} -ne ${#TIME_DT[*]} ];              then echo "Error: Wrong array size (TIME_DT).";              exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#TIME_DT_ATMOS_DYN[*]} ];    then echo "Error: Wrong array size (TIME_DT_ATMOS_DYN).";    exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#TIME_DT_ATMOS_PHY_MP[*]} ]; then echo "Error: Wrong array size (TIME_DT_ATMOS_PHY_MP)."; exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#TIME_DT_ATMOS_PHY_RD[*]} ]; then echo "Error: Wrong array size (TIME_DT_ATMOS_PHY_RD)."; exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#TIME_DT_ATMOS_PHY_SF[*]} ]; then echo "Error: Wrong array size (TIME_DT_ATMOS_PHY_SF)."; exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#TIME_DT_ATMOS_PHY_TB[*]} ]; then echo "Error: Wrong array size (TIME_DT_ATMOS_PHY_TB)."; exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#TIME_DT_OCEAN[*]} ];        then echo "Error: Wrong array size (TIME_DT_OCEAN).";        exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#TIME_DT_LAND[*]} ];         then echo "Error: Wrong array size (TIME_DT_LAND).";         exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#TIME_DT_URBAN[*]} ];        then echo "Error: Wrong array size (TIME_DT_URBAN).";        exit 1; fi

if [ ${NUM_DOMAIN} -ne ${#PRC_NUM_X[*]} ]; then echo "Error: Wrong array size (PRC_NUM_X)."; exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#PRC_NUM_Y[*]} ]; then echo "Error: Wrong array size (PRC_NUM_Y)."; exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#KMAX[*]} ];      then echo "Error: Wrong array size (KMAX).";      exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#IMAX[*]} ];      then echo "Error: Wrong array size (IMAX).";      exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#JMAX[*]} ];      then echo "Error: Wrong array size (JMAX).";      exit 1; fi
 
if [ ${LKMAX} -ne ${#LDZ[*]} ]; then echo "Error: Wrong array size (LDZ)."; exit 1; fi
if [ ${UKMAX} -ne ${#UDZ[*]} ]; then echo "Error: Wrong array size (UDZ)."; exit 1; fi

if [ ${NUM_DOMAIN} -ne ${#DX[*]} ];    then echo "Error: Wrong array size (DX).";    exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#DY[*]} ];    then echo "Error: Wrong array size (DY).";    exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#DEF_Z[*]} ]; then echo "Error: Wrong array size (DEF_Z)."; exit 1; fi

if [ ${NUM_DOMAIN} -ne ${#BUFFER_DZ[*]} ]; then echo "Error: Wrong array size (BUFFER_DZ)."; exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#BUFFER_DX[*]} ]; then echo "Error: Wrong array size (BUFFER_DX)."; exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#BUFFER_DY[*]} ]; then echo "Error: Wrong array size (BUFFER_DY)."; exit 1; fi

if [ ${NUM_DOMAIN} -ne ${#ATMOS_DYN_TYPE[*]} ];    then echo "Error: Wrong array size (ATMOS_DYN_TYPE).";    exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#ATMOS_PHY_MP_TYPE[*]} ]; then echo "Error: Wrong array size (ATMOS_PHY_MP_TYPE)."; exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#ATMOS_PHY_RD_TYPE[*]} ]; then echo "Error: Wrong array size (ATMOS_PHY_RD_TYPE)."; exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#ATMOS_PHY_SF_TYPE[*]} ]; then echo "Error: Wrong array size (ATMOS_PHY_SF_TYPE)."; exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#ATMOS_PHY_TB_TYPE[*]} ]; then echo "Error: Wrong array size (ATMOS_PHY_TB_TYPE)."; exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#OCEAN_TYPE[*]} ];        then echo "Error: Wrong array size (OCEAN_TYPE).";        exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#LAND_TYPE[*]} ];         then echo "Error: Wrong array size (LAND_TYPE).";         exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#URBAN_TYPE[*]} ];        then echo "Error: Wrong array size (URBAN_TYPE).";        exit 1; fi

if [ ${NUM_DOMAIN} -ne ${#TOPOTYPE[*]} ];     then echo "Error: Wrong array size (TOPOTYPE).";     exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#LANDUSETYPE[*]} ];  then echo "Error: Wrong array size (LANDUSETYPE).";  exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#MAXSLOPE[*]} ];     then echo "Error: Wrong array size (MAXSLOPE).";     exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#MAXSLOPE_BND[*]} ]; then echo "Error: Wrong array size (MAXSLOPE_BND)."; exit 1; fi

#################################################
#
# set common parameters
#
#################################################

PP_CONF="conf/base.pp.conf.sh"
INIT_CONF="conf/base.init.conf.sh"
RUN_CONF="conf/base.run.conf.sh"
PARAM_CONF="conf/param.conf.sh"
LAUNCH_CONF="conf/base.launch.conf.sh"
JOB_SHELL="conf/base.jobshell.sh"
LINK_SHELL="conf/base.linkinput.nicam.sh"

LANDPARAMDIR="data"
LAND_PROPERTY_IN_FILENAME="${LANDPARAMDIR}/param.bucket.conf"

LATLON_CATALOGUEDIR="data"
LATLON_CATALOGUE_FNAME="${LATLON_CATALOGUEDIR}/latlon_domain_catalogue.txt"

RADDIR="data"
RD_MSTRN_GASPARA_IN_FILENAME="${RADDIR}/PARAG.29"
RD_MSTRN_AEROPARA_IN_FILENAME="${RADDIR}/PARAPC.29"
RD_MSTRN_HYGROPARA_IN_FILENAME="${RADDIR}/VARDATA.RM29"
RD_PROFILE_CIRA86_IN_FILENAME="${RADDIR}/cira.nc"
RD_PROFILE_MIPAS2001_IN_BASENAME="${RADDIR}/MIPAS"

# set time parameters
TIME_STARTDATE="${STARTDATE[0]}, ${STARTDATE[1]}, ${STARTDATE[2]}, ${STARTDATE[3]}, ${STARTDATE[4]}, ${STARTDATE[5]}"
TIME_STARTMS="${STARTDATE[6]}.D0"

TIME1=`date --utc --date="${STARTDATE[0]}-${STARTDATE[1]}-${STARTDATE[2]} ${STARTDATE[3]}:${STARTDATE[4]:${STARTDATE[5]}}" +%s`
TIME0=`date --utc --date="${STARTDATE[0]}-01-01 00:00:00" +%s`
DAYSEC=`expr ${TIME1} - ${TIME0}`
INITSEC=`printf "%011d" ${DAYSEC}`.`printf "%03d" ${STARTDATE[6]}`

IFS="," eval 'LIST_LDZ="${LDZ[*]}"'
IFS="," eval 'LIST_UDZ="${UDZ[*]}"'

DNUM=1
while [ $DNUM -le $NUM_DOMAIN ]
do ### domain loop
D=`expr $DNUM - 1`

FNUM=`printf "%02d" $DNUM`

# set numbers of domain process
eval 'PRC_DOMAINS[$D]=`expr ${PRC_NUM_X[$D]} \* ${PRC_NUM_Y[$D]}`'

# set names of run.conf
eval 'CONF_FILES[$D]="run.d${FNUM}.conf"'

# set vertical axis
LINE_Z="${DEF_Z[$D]}"

# set filenames for each domain
PP_IO_LOG_BASENAME="${INDIR}/pp_LOG_d${FNUM}"
TOPO_IN_BASENAME="${INDIR}/topo_d${FNUM}"
LANDUSE_IN_BASENAME="${INDIR}/landuse_d${FNUM}"

TOPODIR="${SCALE_DB}/topo"
TOPO_IN_CATALOGUE="${TOPOTYPE[$D]}_catalogue.txt"
TOPO_IN_DIR="${TOPODIR}/${TOPOTYPE[$D]}/Products"

LANDUSEDIR="${SCALE_DB}/landuse"
LANDUSE_IN_CATALOGUE="${LANDUSETYPE[$D]}_catalogue.txt"
LANDUSE_IN_DIR="${LANDUSEDIR}/${LANDUSETYPE[$D]}/Products"

INIT_IO_LOG_BASENAME="${INDIR}/init_LOG_d${FNUM}"
INIT_RESTART_OUT_BASENAME="${INDIR}/init_d${FNUM}"
BASENAME_BOUNDARY="${INDIR}/boundary_d${FNUM}"

RUN_IO_LOG_BASENAME="${OUTDIR}/LOG_d${FNUM}"
RESTART_OUT_BASENAME="${OUTDIR}/restart_d${FNUM}"
HISTORY_DEFAULT_BASENAME="${OUTDIR}/history_d${FNUM}"

if [ ${TRIMDIR} = true ]; then
  LAND_PROPERTY_IN_FILENAME="${LAND_PROPERTY_IN_FILENAME##*/}"
  LATLON_CATALOGUE_FNAME="${LATLON_CATALOGUE_FNAME##*/}"

  PP_IO_LOG_BASENAME="${PP_IO_LOG_BASENAME##*/}"
  TOPO_IN_BASENAME="${TOPO_IN_BASENAME##*/}"
  LANDUSE_IN_BASENAME="${LANDUSE_IN_BASENAME##*/}"

  INIT_IO_LOG_BASENAME="${INIT_IO_LOG_BASENAME##*/}"
  INIT_RESTART_OUT_BASENAME="${INIT_RESTART_OUT_BASENAME##*/}"
  BASENAME_BOUNDARY="${BASENAME_BOUNDARY##*/}"

  RUN_IO_LOG_BASENAME="${RUN_IO_LOG_BASENAME##*/}"
  HISTORY_DEFAULT_BASENAME="${HISTORY_DEFAULT_BASENAME##*/}"
  RESTART_OUT_BASENAME="${RESTART_OUT_BASENAME##*/}"

  RD_MSTRN_GASPARA_IN_FILENAME="${RD_MSTRN_GASPARA_IN_FILENAME##*/}"
  RD_MSTRN_AEROPARA_IN_FILENAME="${RD_MSTRN_AEROPARA_IN_FILENAME##*/}"
  RD_MSTRN_HYGROPARA_IN_FILENAME="${RD_MSTRN_HYGROPARA_IN_FILENAME##*/}"
  RD_PROFILE_CIRA86_IN_FILENAME="${RD_PROFILE_CIRA86_IN_FILENAME##*/}"
  RD_PROFILE_MIPAS2001_IN_BASENAME="${RD_PROFILE_MIPAS2001_IN_BASENAME##*/}"
fi

# copy parameters
TOPO_OUT_BASENAME="${TOPO_IN_BASENAME}"
LANDUSE_OUT_BASENAME="${LANDUSE_IN_BASENAME}"
RUN_RESTART_IN_BASENAME="${INIT_RESTART_OUT_BASENAME}_${INITSEC}"
ATMOS_BOUNDARY_IN_BASENAME="${BASENAME_BOUNDARY}"
ATMOS_BOUNDARY_START_DATE="${TIME_STARTDATE}"

# set nesting parameters
if [ $DNUM -lt $NUM_DOMAIN ]; then
  IAM_PARENT=".true."
else
  IAM_PARENT=".false."
fi
if [ $DNUM -gt 1 ]; then
  IAM_DAUGHTER=".true."
else
  IAM_DAUGHTER=".false."
fi

# set boundary parameters
if [ $DNUM -gt 1 ]; then
  ATMOS_BOUNDARY_IN_BASENAME=""
  ATMOS_BOUNDARY_START_DATE=""
  ATMOS_BOUNDARY_UPDATE_DT="0.D0"
  ATMOS_BOUNDARY_USE_DENS=".false."
  ATMOS_BOUNDARY_USE_QHYD=".true."
else
  ATMOS_BOUNDARY_UPDATE_DT="${TIME_DT_BOUNDARY}"
  ATMOS_BOUNDARY_USE_DENS=".true."
  ATMOS_BOUNDARY_USE_QHYD=".false."
fi

#################################################
#
# make config files
#
#################################################

source ${PP_CONF}
source ${INIT_CONF}
source ${RUN_CONF}
source ${PARAM_CONF}

cat base.pp.conf \
    param.region.conf \
    param.admin.conf \
> pp.d${FNUM}.conf

cat base.init.conf \
    param.region.conf \
    param.admin.conf \
> init.d${FNUM}.conf

cat base.run.conf \
    param.region.conf \
    param.admin.conf \
    param.physics.conf \
    param.history.conf \
> run.d${FNUM}.conf

rm -f base.*.conf param.*.conf

DNUM=`expr $DNUM + 1`
done ### end domain loop

#################################################
#
# make launcher
#
#################################################

IFS="," eval 'LIST_PRC_DOMAINS="${PRC_DOMAINS[*]}"'
IFS="," eval 'LIST_CONF_FILES="${CONF_FILES[*]}"'

source ${LAUNCH_CONF}
mv -f base.launch.conf launch.conf

#################################################
#
# make jobshell
#
#################################################

TOTALNODE=0
for NODE in ${PRC_DOMAINS[*]}
do
  TOTALNODE=`expr ${TOTALNODE} + ${NODE}`
done

if [ ${TOTALNODE} -ge 36865 ]; then
  RSCGRP="huge"
elif [ ${TOTALNODE} -ge 385 ]; then
  RSCGRP="large"
else
  RSCGRP="small"
fi

source ${JOB_SHELL}
mv -f base.jobshell run.nest.sh

#################################################
#
# make initial link shell
#
#################################################

source ${LINK_SHELL}
mv -f base.linkinput linkinput.sh

#################################################
