#!/bin/bash
#
#################################################
#
# SCALE-LES configuration file generator
#
#################################################

USERDEF_FILE="./user-def.sh"

if [ ! -f "${USERDEF_FILE}" ]; then
  echo "Error: User defined file was not found: "${USERDEF_FILE}
  exit 1
else
  # include user defined variables
  source ${USERDEF_FILE}
fi

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

INPUT_CONFIGDIR="input"
OUTPUT_CONFIGDIR="output"

PP_CONF="${INPUT_CONFIGDIR}/base.pp.conf.sh"
INIT_CONF="${INPUT_CONFIGDIR}/base.init.conf.sh"
RUN_CONF="${INPUT_CONFIGDIR}/base.run.conf.sh"
PARAM_CONF="${INPUT_CONFIGDIR}/param.conf.sh"
LAUNCH_CONF="${INPUT_CONFIGDIR}/base.launch.conf.sh"
JOB_SHELL="${INPUT_CONFIGDIR}/base.jobshell.sh"
LINK_SHELL="${INPUT_CONFIGDIR}/base.linkinput.nicam.sh"

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

eval 'STARTDATE[0]=`printf "%04d" ${STARTDATE[0]}`'
eval 'STARTDATE[1]=`printf "%02d" ${STARTDATE[1]}`'
eval 'STARTDATE[2]=`printf "%02d" ${STARTDATE[2]}`'
eval 'STARTDATE[3]=`printf "%02d" ${STARTDATE[3]}`'
eval 'STARTDATE[4]=`printf "%02d" ${STARTDATE[4]}`'
eval 'STARTDATE[5]=`printf "%02d" ${STARTDATE[5]}`'
eval 'STARTDATE[6]=`printf "%03d" ${STARTDATE[6]}`'

echo "STARTDATE: ${STARTDATE[0]}/${STARTDATE[1]}/${STARTDATE[2]} ${STARTDATE[3]}:${STARTDATE[4]}:${STARTDATE[5]}.${STARTDATE[6]}"

TIME1=`date --utc --date="${STARTDATE[0]}-${STARTDATE[1]}-${STARTDATE[2]} ${STARTDATE[3]}:${STARTDATE[4]:${STARTDATE[5]}}" +%s`
TIME0=`date --utc --date="${STARTDATE[0]}-01-01 00:00:00" +%s`
DAYSEC=`expr ${TIME1} - ${TIME0}`
INITSEC=`printf "%011d" ${DAYSEC}`.${STARTDATE[6]}

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
  ATMOS_BOUNDARY_USE_QHYD=".true."
else
  ATMOS_BOUNDARY_UPDATE_DT="${TIME_DT_BOUNDARY}"
  ATMOS_BOUNDARY_USE_QHYD=".false."
fi

#################################################
#
# make config files
#
#################################################

mkdir -p ${OUTPUT_CONFIGDIR}

source ${PP_CONF}
source ${INIT_CONF}
source ${RUN_CONF}
source ${PARAM_CONF}

cat base.pp.conf \
    param.region.conf \
    param.admin.conf \
> ${OUTPUT_CONFIGDIR}/pp.d${FNUM}.conf

cat base.init.conf \
    param.region.conf \
    param.admin.conf \
> ${OUTPUT_CONFIGDIR}/init.d${FNUM}.conf

cat base.run.conf \
    param.region.conf \
    param.admin.conf \
    param.physics.conf \
    param.history.conf \
> ${OUTPUT_CONFIGDIR}/run.d${FNUM}.conf

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
mv -f base.launch.conf ${OUTPUT_CONFIGDIR}/launch.conf

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
mv -f base.jobshell ${OUTPUT_CONFIGDIR}/run.nest.sh

#################################################
#
# make initial link shell
#
#################################################

source ${LINK_SHELL}
mv -f base.linkinput ${OUTPUT_CONFIGDIR}/linkinput.sh

#################################################
