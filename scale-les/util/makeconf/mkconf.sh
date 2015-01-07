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
STARTDATE=( 1999  5 16  6  0  0  0 ) # Year Mon Day Hour Min Sec milli-Sec
#STARTDATE=( 1999  6 10  6  0  0  0 ) # Year Mon Day Hour Min Sec milli-Sec
#STARTDATE=( 1999  7 29  6  0  0  0 ) # Year Mon Day Hour Min Sec milli-Sec
#STARTDATE=( 1999  8 21  6  0  0  0 ) # Year Mon Day Hour Min Sec milli-Sec

TIME_DURATION="604800.D0"
TIME_UNIT="SEC"

TIME_DT_RESTART="86400.D0"
TIME_DT_HISTORY="900.D0"
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
"FZ(:) =     80.841D0,    248.821D0,    429.882D0,    625.045D0,    835.409D0,
           1062.158D0,   1306.565D0,   1570.008D0,   1853.969D0,   2160.047D0,
           2489.963D0,   2845.575D0,   3228.883D0,   3642.044D0,   4087.384D0,
           4567.409D0,   5084.820D0,   5642.530D0,   6243.676D0,   6891.642D0,
           7590.074D0,   8342.904D0,   9154.367D0,  10029.028D0,  10971.815D0,
          11988.030D0,  13083.390D0,  14264.060D0,  15536.685D0,  16908.430D0,
          18387.010D0,  19980.750D0,  21698.615D0,  23550.275D0,  25546.155D0,
          28113.205D0,"

"FZ(:) =    80.0000D0,   160.0000D0,   240.0000D0,   320.0000D0,   400.0000D0,
           485.4400D0,   576.6899D0,   674.1449D0,   778.2268D0,   889.3863D0,
          1008.1046D0,  1134.8958D0,  1270.3087D0,  1414.9298D0,  1569.3851D0,
          1734.3434D0,  1910.5188D0,  2098.6741D0,  2299.6240D0,  2514.2385D0,
          2743.4468D0,  2988.2412D0,  3249.6816D0,  3528.9001D0,  3827.1055D0,
          4145.5889D0,  4485.7290D0,  4848.9985D0,  5236.9702D0,  5651.3242D0,
          6093.8545D0,  6566.4771D0,  7071.2378D0,  7610.3223D0,  8186.0645D0,
          8800.9570D0,  9457.6621D0, 10159.0234D0, 10908.0771D0, 11708.0664D0,
         12508.0664D0, 13308.0664D0, 14108.0664D0, 14908.0664D0, 15708.0664D0,
         16508.0664D0, 17308.0664D0, 18108.0664D0, 18908.0664D0, 19708.0664D0,
         20508.0664D0, 21308.0664D0, 22108.0664D0, 22908.0664D0, 23708.0664D0,
         24508.0664D0, 25308.0664D0, 26108.0664D0, 26908.0664D0, 27708.0664D0,"

"FZ(:) =    60.0000D0,   120.0000D0,   180.0000D0,   240.0000D0,   300.0000D0,
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

OCEAN_TYPE=( "SLAB" "SLAB" "SLAB" ) # required num. of parameters for each domain
LAND_TYPE=(  "SLAB" "SLAB" "SLAB" ) # required num. of parameters for each domain
URBAN_TYPE=( "SLC"  "SLC"  "SLC"  ) # required num. of parameters for each domain

# history parameters
HIST_ITEMS_SNAPSHOT=(
# ATMOS history
  "DENS" "MOMZ" "MOMX" "MOMY"
  "QV" "QC" "QR" "QI" "QS" "QG" "QHYD" "QLIQ" "QICE"
  "T" "PRES" "U" "V" "W" "PT" "RH" "RHL" "RHI"
  "SFC_PRES" "SFC_TEMP" "SFC_Z0M"
  "U10" "V10" "T2" "Q2" "LHFLX" "SHFLX" "GHFLX"
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
HIST_ITEMS_AVERAGE=(
  "PREC"
)

##### for init
BASENAME_ORG=""
FILETYPE_ORG="NICAM-NETCDF"
NUMBER_OF_FILES="8"
WRF_FILE_TYPE=""

INIT_USE_NESTING=".false."
PARENT_PRC_NUM_X=""        # used if INIT_USE_NESTING = .true.
PARENT_PRC_NUM_Y=""        # used if INIT_USE_NESTING = .true.
PARENT_KMAX=""             # used if INIT_USE_NESTING = .true.
PARENT_IMAX=""             # used if INIT_USE_NESTING = .true.
PARENT_JMAX=""             # used if INIT_USE_NESTING = .true.
PARENT_LKMAX=""            # used if INIT_USE_NESTING = .true.

##### for pp
TOPOTYPE=(    "GTOPO30" "GTOPO30" "DEM50M" ) # required num. of parameters for each domain: GTOPO30, DEM50M
LANDUSETYPE=( "GLCCv2"  "GLCCv2"  "LU100M" ) # required num. of parameters for each domain: GLCCv2, LU100M

MAXSLOPE=(     "0.5D0" "1.5D0" "7.0D0" ) # required num. of parameters for each domain
MAXSLOPE_BND=( "0.5D0" "0.5D0" "1.5D0" ) # required num. of parameters for each domain

##### for staging
TRIMDIR=true

INDIR="input"   # used if TRIMDIR = false
OUTDIR="output" # used if TRIMDIR = false

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
  ATMOS_BOUNDARY_USE_QHYD=".true."
  ATMOS_BOUNDARY_IN_BASENAME=""
  ATMOS_BOUNDARY_UPDATE_DT="0.D0"
else
  ATMOS_BOUNDARY_USE_QHYD=".false."
  ATMOS_BOUNDARY_UPDATE_DT="${TIME_DT_BOUNDARY}"
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
