#!/bin/bash

#################################################
#
# SCALE-LES configuration file generator
#
#################################################

if [ "${USERDEF_FILE}x" = "x" ]; then
  USERDEF_FILE="./USER.sh"
fi
echo ${OUTPUT_CONFIGDIR}
if [ "${OUTPUT_CONFIGDIR}x" = "x" ]; then
  OUTPUT_CONFIGDIR="experiment"
fi

if [ ! -f "${USERDEF_FILE}" ]; then
  echo "Error: User defined file was not found: "${USERDEF_FILE}
  exit 1
else
  # include user defined variables
  source ${USERDEF_FILE}
fi

# debug information
FMT_YEAR=`printf '%04d' ${RUN_DATE_YEAR}`
FMT_MON=`printf '%02d' ${RUN_DATE_MON}`
FMT_DAY=`printf '%02d' ${RUN_DATE_DAY}`
FMT_HOUR=`printf '%02d' ${RUN_DATE_HOUR}`
FMT_MIN=`printf '%02d' ${RUN_DATE_MIN}`
FMT_SEC=`printf '%02d' ${RUN_DATE_SEC}`
FMT_MSEC=`printf '%03d' ${RUN_DATE_MSEC}`

echo ""
echo "#################################################"
echo "#     SCALE-RM Configuration File Generator     #"
echo "#################################################"
echo ""
echo "START DATE: ${FMT_YEAR}/${FMT_MON}/${FMT_DAY} - ${FMT_HOUR}:${FMT_MIN}:${FMT_SEC}.${FMT_MSEC}"
echo ""

#################################################
#
# Support functions
#
#################################################

multiply() {
  factor=$1
  shift
  for x in "$@"; do
    bc -l <<< "$factor * $x"
  done
}

replicate() {
  n=$1
  value=$2
  for (( i = 0; i < n; ++i )); do
    echo "$value"
  done
}

#################################################
#
# check parameters
#
#################################################

if [ ! -n "${SCALE_DB-}" ]; then echo "Error: SCALE_DB is not defined. Check!"; exit 1; fi

if [ ${NUM_DOMAIN} -ne ${#TIME_DT[*]} ];              then echo "Error: Wrong array size (TIME_DT).";              exit 1; fi

if [[ -z ${TIME_DT_ATMOS_DYN+x} ]];    then TIME_DT_ATMOS_DYN=( $(multiply 0.5 "${TIME_DT[@]}") );   fi
if [[ -z ${TIME_DT_ATMOS_PHY_CP+x} ]]; then TIME_DT_ATMOS_PHY_CP=( "${TIME_DT[@]}" );                fi
if [[ -z ${TIME_DT_ATMOS_PHY_MP+x} ]]; then TIME_DT_ATMOS_PHY_MP=( "${TIME_DT[@]}" );                fi
if [[ -z ${TIME_DT_ATMOS_PHY_RD+x} ]]; then TIME_DT_ATMOS_PHY_RD=( $(multiply 10 "${TIME_DT[@]}") ); fi
if [[ -z ${TIME_DT_ATMOS_PHY_SF+x} ]]; then TIME_DT_ATMOS_PHY_SF=( "${TIME_DT[@]}" );                fi
if [[ -z ${TIME_DT_ATMOS_PHY_TB+x} ]]; then TIME_DT_ATMOS_PHY_TB=( "${TIME_DT[@]}" );                fi
if [[ -z ${TIME_DT_ATMOS_PHY_BL+x} ]]; then TIME_DT_ATMOS_PHY_BL=( "${TIME_DT[@]}" );                fi
if [[ -z ${TIME_DT_OCEAN+x} ]];        then TIME_DT_OCEAN=( "${TIME_DT[@]}" );                       fi
if [[ -z ${TIME_DT_LAND+x} ]];         then TIME_DT_LAND=(  "${TIME_DT[@]}" );                       fi
if [[ -z ${TIME_DT_URBAN+x} ]];        then TIME_DT_URBAN=( "${TIME_DT[@]}" );                       fi

if [ ${NUM_DOMAIN} -ne ${#TIME_DT_ATMOS_DYN[*]} ];    then echo "Error: Wrong array size (TIME_DT_ATMOS_DYN).";    exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#TIME_DT_ATMOS_PHY_CP[*]} ]; then echo "Error: Wrong array size (TIME_DT_ATMOS_PHY_CP)."; exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#TIME_DT_ATMOS_PHY_MP[*]} ]; then echo "Error: Wrong array size (TIME_DT_ATMOS_PHY_MP)."; exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#TIME_DT_ATMOS_PHY_RD[*]} ]; then echo "Error: Wrong array size (TIME_DT_ATMOS_PHY_RD)."; exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#TIME_DT_ATMOS_PHY_SF[*]} ]; then echo "Error: Wrong array size (TIME_DT_ATMOS_PHY_SF)."; exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#TIME_DT_ATMOS_PHY_TB[*]} ]; then echo "Error: Wrong array size (TIME_DT_ATMOS_PHY_TB)."; exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#TIME_DT_ATMOS_PHY_BL[*]} ]; then echo "Error: Wrong array size (TIME_DT_ATMOS_PHY_BL)."; exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#TIME_DT_OCEAN[*]} ];        then echo "Error: Wrong array size (TIME_DT_OCEAN).";        exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#TIME_DT_LAND[*]} ];         then echo "Error: Wrong array size (TIME_DT_LAND).";         exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#TIME_DT_URBAN[*]} ];        then echo "Error: Wrong array size (TIME_DT_URBAN).";        exit 1; fi

if [ ${NUM_DOMAIN} -ne ${#PRC_NUM_X[*]} ]; then echo "Error: Wrong array size (PRC_NUM_X)."; exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#PRC_NUM_Y[*]} ]; then echo "Error: Wrong array size (PRC_NUM_Y)."; exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#KMAX[*]} ];      then echo "Error: Wrong array size (KMAX).";      exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#IMAXG[*]} ];     then echo "Error: Wrong array size (IMAXG).";     exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#JMAXG[*]} ];     then echo "Error: Wrong array size (JMAXG).";     exit 1; fi

if [ ${OKMAX} -ne ${#ODZ[*]} ]; then echo "Error: Wrong array size (ODZ)."; exit 1; fi
if [ ${LKMAX} -ne ${#LDZ[*]} ]; then echo "Error: Wrong array size (LDZ)."; exit 1; fi
if [ ${UKMAX} -ne ${#UDZ[*]} ]; then echo "Error: Wrong array size (UDZ)."; exit 1; fi

if [ ${NUM_DOMAIN} -ne ${#DX[*]} ];    then echo "Error: Wrong array size (DX).";    exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#DY[*]} ];    then echo "Error: Wrong array size (DY).";    exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#DEF_Z[*]} ]; then echo "Error: Wrong array size (DEF_Z)."; exit 1; fi

if [[ $NUM_DOMAIN -gt 1 ]]; then
  if [[ ${#BUFFER_DZ[@]} -eq 1 ]]; then BUFFER_DZ=( $(replicate "$NUM_DOMAIN" "$BUFFER_DZ") ); fi
fi
if [[ -z ${BUFFER_DX+x} ]]; then BUFFER_DX=( $(multiply 20 "${DX[@]}") ); fi
if [[ -z ${BUFFER_DY+x} ]]; then BUFFER_DY=( $(multiply 20 "${DY[@]}") ); fi
if [ ${NUM_DOMAIN} -ne ${#BUFFER_DZ[*]} ]; then echo "Error: Wrong array size (BUFFER_DZ)."; exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#BUFFER_DX[*]} ]; then echo "Error: Wrong array size (BUFFER_DX)."; exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#BUFFER_DY[*]} ]; then echo "Error: Wrong array size (BUFFER_DY)."; exit 1; fi

if [[ $NUM_DOMAIN -gt 1 ]]; then
  if [[ ${#MAPPROJECTION_BASEPOINT_LON[@]} -eq 1 ]]; then MAPPROJECTION_BASEPOINT_LON=( $(replicate "$NUM_DOMAIN" "$MAPPROJECTION_BASEPOINT_LON") ); fi
  if [[ ${#MAPPROJECTION_BASEPOINT_LAT[@]} -eq 1 ]]; then MAPPROJECTION_BASEPOINT_LAT=( $(replicate "$NUM_DOMAIN" "$MAPPROJECTION_BASEPOINT_LAT") ); fi
fi
if [ ${NUM_DOMAIN} -ne ${#MAPPROJECTION_BASEPOINT_LON[*]} ]; then echo "Error: Wrong array size (MAPPROJECTION_BASEPOINT_LON)."; exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#MAPPROJECTION_BASEPOINT_LAT[*]} ]; then echo "Error: Wrong array size (MAPPROJECTION_BASEPOINT_LAT)."; exit 1; fi

if [[ $NUM_DOMAIN -gt 1 ]]; then
  if [[ ${#ATMOS_DYN_TYPE[@]} -eq 1 ]];    then ATMOS_DYN_TYPE=( $(replicate "$NUM_DOMAIN" "$ATMOS_DYN_TYPE") );       fi
  if [[ ${#ATMOS_PHY_CP_TYPE[@]} -eq 1 ]]; then ATMOS_PHY_CP_TYPE=( $(replicate "$NUM_DOMAIN" "$ATMOS_PHY_CP_TYPE") ); fi
  if [[ ${#ATMOS_PHY_MP_TYPE[@]} -eq 1 ]]; then ATMOS_PHY_MP_TYPE=( $(replicate "$NUM_DOMAIN" "$ATMOS_PHY_MP_TYPE") ); fi
  if [[ ${#ATMOS_PHY_RD_TYPE[@]} -eq 1 ]]; then ATMOS_PHY_RD_TYPE=( $(replicate "$NUM_DOMAIN" "$ATMOS_PHY_RD_TYPE") ); fi
  if [[ ${#ATMOS_PHY_SF_TYPE[@]} -eq 1 ]]; then ATMOS_PHY_SF_TYPE=( $(replicate "$NUM_DOMAIN" "$ATMOS_PHY_SF_TYPE") ); fi
  if [[ ${#ATMOS_PHY_TB_TYPE[@]} -eq 1 ]]; then ATMOS_PHY_TB_TYPE=( $(replicate "$NUM_DOMAIN" "$ATMOS_PHY_TB_TYPE") ); fi
  if [[ ${#ATMOS_PHY_BL_TYPE[@]} -eq 1 ]]; then ATMOS_PHY_BL_TYPE=( $(replicate "$NUM_DOMAIN" "$ATMOS_PHY_BL_TYPE") ); fi
  if [[ ${#OCEAN_DYN_TYPE[@]} -eq 1 ]];    then OCEAN_DYN_TYPE=( $(replicate "$NUM_DOMAIN" "$OCEAN_DYN_TYPE") );       fi
  if [[ ${#LAND_DYN_TYPE[@]} -eq 1 ]];     then LAND_DYN_TYPE=( $(replicate "$NUM_DOMAIN" "$LAND_DYN_TYPE") );         fi
  if [[ ${#LAND_SFC_TYPE[@]} -eq 1 ]];     then LAND_SFC_TYPE=( $(replicate "$NUM_DOMAIN" "$LAND_SFC_TYPE") );         fi
  if [[ ${#URBAN_DYN_TYPE[@]} -eq 1 ]];    then URBAN_DYN_TYPE=( $(replicate "$NUM_DOMAIN" "$URBAN_DYN_TYPE") );       fi
fi
if [ ${NUM_DOMAIN} -ne ${#ATMOS_DYN_TYPE[*]} ];    then echo "Error: Wrong array size (ATMOS_DYN_TYPE).";    exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#ATMOS_PHY_CP_TYPE[*]} ]; then echo "Error: Wrong array size (ATMOS_PHY_CP_TYPE)."; exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#ATMOS_PHY_MP_TYPE[*]} ]; then echo "Error: Wrong array size (ATMOS_PHY_MP_TYPE)."; exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#ATMOS_PHY_RD_TYPE[*]} ]; then echo "Error: Wrong array size (ATMOS_PHY_RD_TYPE)."; exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#ATMOS_PHY_SF_TYPE[*]} ]; then echo "Error: Wrong array size (ATMOS_PHY_SF_TYPE)."; exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#ATMOS_PHY_TB_TYPE[*]} ]; then echo "Error: Wrong array size (ATMOS_PHY_TB_TYPE)."; exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#ATMOS_PHY_BL_TYPE[*]} ]; then echo "Error: Wrong array size (ATMOS_PHY_BL_TYPE)."; exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#OCEAN_DYN_TYPE[*]} ];    then echo "Error: Wrong array size (OCEAN_DYN_TYPE).";    exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#LAND_DYN_TYPE[*]} ];     then echo "Error: Wrong array size (LAND_DYN_TYPE).";     exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#LAND_SFC_TYPE[*]} ];     then echo "Error: Wrong array size (LAND_SFC_TYPE).";     exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#URBAN_DYN_TYPE[*]} ];    then echo "Error: Wrong array size (URBAN_DYN_TYPE).";    exit 1; fi

if [[ $NUM_DOMAIN -gt 1 ]]; then
  if [[ ${#TOPOTYPE[@]} -eq 1 ]];    then TOPOTYPE=( $(replicate "$NUM_DOMAIN" "$TOPOTYPE") );       fi
  if [[ ${#LANDUSETYPE[@]} -eq 1 ]]; then LANDUSETYPE=( $(replicate "$NUM_DOMAIN" "$LANDUSETYPE") ); fi
fi
if [ ${NUM_DOMAIN} -ne ${#TOPOTYPE[*]} ];    then echo "Error: Wrong array size (TOPOTYPE).";     exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#LANDUSETYPE[*]} ]; then echo "Error: Wrong array size (LANDUSETYPE).";  exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#COPYTOPO[*]} ];    then echo "Error: Wrong array size (COPYTOPO).";     exit 1; fi
if [ ${NUM_DOMAIN} -ne ${#SMOOTH_LOCAL[*]} ]; then echo "Error: Wrong array size (COPYTOPO).";    exit 1; fi

if [[ -z ${ONLINE_AGGRESSIVE_COMM+x} ]]; then ONLINE_AGGRESSIVE_COMM=".false."; fi

#################################################
#
# set common parameters
#
#################################################

# set formatted date
YEAR=`printf '%1.f' ${RUN_DATE_YEAR}`
MON=`printf '%1.f' ${RUN_DATE_MON}`
DAY=`printf '%1.f' ${RUN_DATE_DAY}`
HOUR=`printf '%1.f' ${RUN_DATE_HOUR}`
MIN=`printf '%1.f' ${RUN_DATE_MIN}`
SEC=`printf '%1.f' ${RUN_DATE_SEC}`
MSEC=`printf '%1.f' ${RUN_DATE_MSEC}`

STARTDATE=( ${YEAR} ${MON} ${DAY} ${HOUR} ${MIN} ${SEC} ${MSEC} )

# set number of files for boundary
case ${TIME_DURATION_UNIT} in
  "DAY"  ) TIME_DURATION_UNIT_SEC=86400 ;;
  "HOUR" ) TIME_DURATION_UNIT_SEC=3600  ;;
  "MIN"  ) TIME_DURATION_UNIT_SEC=60    ;;
  "SEC"  ) TIME_DURATION_UNIT_SEC=1     ;;
  *      ) echo 'Error: Undefined type of TIME_DURATION_UNIT.' ;;
esac
INT_DURATION=`expr ${TIME_DURATION%%.*} \* ${TIME_DURATION_UNIT_SEC}`
INT_BOUNDARY_DT=`echo ${TIME_DT_BOUNDARY%%.*}`
if [ ${FILETYPE_ORG} = "NetCDF" ]; then
  NUMBER_OF_FILES=1
  NUMBER_OF_TSTEPS=`expr ${INT_DURATION} / ${INT_BOUNDARY_DT} + 1`
  if [ ${NUMBER_OF_TSTEPS} -le 1 ]; then
    NUMBER_OF_TSTEPS=2
  fi
else
  NUMBER_OF_FILES=`expr ${INT_DURATION} / ${INT_BOUNDARY_DT} + 1`
  if [ ${NUMBER_OF_FILES} -le 1 ]; then
    NUMBER_OF_FILES=2
  fi
  NUMBER_OF_TSTEPS=1
fi

# set maximum number of steps
HISTORY_2D_INTERVAL=`echo ${TIME_DT_HISTORY_2D%%.*}`
HISTORY_3D_INTERVAL=`echo ${TIME_DT_HISTORY_3D%%.*}`
MAXSTEP_2D=`expr ${INT_DURATION} / ${HISTORY_2D_INTERVAL} + 1`
MAXSTEP_3D=`expr ${INT_DURATION} / ${HISTORY_3D_INTERVAL} + 1`

INPUT_CONFIGDIR="config"

PPDIR="../pp"
INITDIR="../init"
RUNDIR="../run"

PP_CONF="${INPUT_CONFIGDIR}/base.pp.conf.sh"
INIT_CONF="${INPUT_CONFIGDIR}/base.init.conf.sh"
RUN_CONF="${INPUT_CONFIGDIR}/base.run.conf.sh"
PARAM_CONF="${INPUT_CONFIGDIR}/param.conf.sh"
SNO_CONF="${INPUT_CONFIGDIR}/base.sno.conf.sh"
LAUNCH_CONF="${INPUT_CONFIGDIR}/base.launch.conf.sh"

# set time parameters
TIME_STARTDATE="${STARTDATE[0]}, ${STARTDATE[1]}, ${STARTDATE[2]}, ${STARTDATE[3]}, ${STARTDATE[4]}, ${STARTDATE[5]}"
TIME_STARTMS="${STARTDATE[6]}.0"

TIME_DT_UNIT="SEC" # SEC only

eval 'STARTDATE[0]=`printf "%04d" ${STARTDATE[0]}`'
eval 'STARTDATE[1]=`printf "%02d" ${STARTDATE[1]}`'
eval 'STARTDATE[2]=`printf "%02d" ${STARTDATE[2]}`'
eval 'STARTDATE[3]=`printf "%02d" ${STARTDATE[3]}`'
eval 'STARTDATE[4]=`printf "%02d" ${STARTDATE[4]}`'
eval 'STARTDATE[5]=`printf "%02d" ${STARTDATE[5]}`'
eval 'STARTDATE[6]=`printf "%03d" ${STARTDATE[6]}`'

INITTIME="${STARTDATE[0]}${STARTDATE[1]}${STARTDATE[2]}-${STARTDATE[3]}${STARTDATE[4]}${STARTDATE[5]}.${STARTDATE[6]}"

IFS="," eval 'LIST_ODZ="${ODZ[*]}"'
IFS="," eval 'LIST_LDZ="${LDZ[*]}"'
IFS="," eval 'LIST_UDZ="${UDZ[*]}"'

if [[ -z ${DOM_NUMBER+x} ]]; then DOM_NUMBER=1; fi

DNUM=1
TPROC=0
while [ $DNUM -le $NUM_DOMAIN ]
do
  D=`expr $DNUM - 1`
  PD=`expr $DNUM - 2`

  DN=`expr $DNUM + $DOM_NUMBER - 1`
  PDNUM=`expr $DN - 1`

  FNUM=`printf "%02d" $DN`
  PFNUM=`printf "%02d" $PDNUM`

  # set numbers of domain process
  eval 'PRC_DOMAINS[$D]=`expr ${PRC_NUM_X[$D]} \* ${PRC_NUM_Y[$D]}`'
  eval 'TPROC=`expr ${TPROC} + ${PRC_DOMAINS[$D]}`'

  # set names of config files
  eval 'PP_CONF_FILES[$D]="pp.d${FNUM}.conf"'
  eval 'INIT_CONF_FILES[$D]="init.d${FNUM}.conf"'
  eval 'RUN_CONF_FILES[$D]="run.d${FNUM}.conf"'
  eval 'SNO_CONF_FILES[$D]="sno.vgridope.d${FNUM}.conf,sno.hgridope.d${FNUM}.conf"'
  eval 'PRC_SNO[$D]="${PRC_DOMAINS[$D]},1"'

  # set vertical axis
  LINE_Z="${DEF_Z[$D]}"

  # set output directory
  if [ -z "${OUT_DIR_PP}" ]; then OUT_DIR_PP="."; fi
  if [ -z "${OUT_DIR_INIT}" ]; then OUT_DIR_INIT="."; fi
  if [ -z "${OUT_DIR_RUN}" ]; then OUT_DIR_RUN="."; fi
  if [ -z "${OUT_DIR_SNO}" ]; then OUT_DIR_SNO="."; fi

  # set filenames for each domain
  PP_IO_LOG_BASENAME="pp_LOG_d${FNUM}"
  TOPOGRAPHY_OUT_BASENAME="${OUT_DIR_PP}/topo_d${FNUM}"
  COPYTOPO_IN_BASENAME="${OUT_DIR_PP}/topo_d${PFNUM}"
  LANDUSE_OUT_BASENAME="${OUT_DIR_PP}/landuse_d${FNUM}"

  TOPO_IN_DIR="${TOPODIR}/${TOPOTYPE[$D]}/Products"
  TOPO_IN_CATALOGUE="${TOPOTYPE[$D]}_catalogue.txt"

  LANDUSE_IN_DIR="${LANDUSEDIR}/${LANDUSETYPE[$D]}/Products"
  LANDUSE_IN_CATALOGUE="${LANDUSETYPE[$D]}_catalogue.txt"

  INIT_TOPOGRAPHY_IN_BASENAME="${PPDIR}/${OUT_DIR_PP}/topo_d${FNUM}"
  INIT_LANDUSE_IN_BASENAME="${PPDIR}/${OUT_DIR_PP}/landuse_d${FNUM}"
  INIT_IO_LOG_BASENAME="init_LOG_d${FNUM}"
  INIT_RESTART_OUT_BASENAME="${OUT_DIR_INIT}/init_d${FNUM}"
  BASENAME_BOUNDARY="${OUT_DIR_INIT}/boundary_d${FNUM}"

  RUN_TOPOGRAPHY_IN_BASENAME="${PPDIR}/${OUT_DIR_PP}/topo_d${FNUM}"
  RUN_LANDUSE_IN_BASENAME="${PPDIR}/${OUT_DIR_PP}/landuse_d${FNUM}"
  RUN_IO_LOG_BASENAME="LOG_d${FNUM}"
  RUN_RESTART_IN_BASENAME="${INITDIR}/${OUT_DIR_INIT}/${INIT_BASENAME}_d${FNUM}_${INITTIME}"
  ATMOS_BOUNDARY_IN_BASENAME="${INITDIR}/${BASENAME_BOUNDARY}"
  RESTART_OUT_BASENAME="${OUT_DIR_RUN}/restart_d${FNUM}"
  FILE_HISTORY_DEFAULT_BASENAME="${OUT_DIR_RUN}/history_d${FNUM}"

  # set nesting parameters
  if [ $DNUM -lt $NUM_DOMAIN ]; then
    IAM_PARENT=".true."
  else
    IAM_PARENT=".false."
  fi
  if [ $DNUM -gt 1 ]; then
    IAM_DAUGHTER=".true."
    PARENT_BASENAME=${COPYTOPO_IN_BASENAME}
    PARENT_PRC_NUM_X=${PRC_NUM_X[$PD]}
    PARENT_PRC_NUM_Y=${PRC_NUM_Y[$PD]}
  else
    IAM_DAUGHTER=".false."
  fi

  # set boundary parameters
  if [ $DNUM -gt 1 ]; then
    ATMOS_BOUNDARY_IN_BASENAME=""
    ATMOS_BOUNDARY_USE_QHYD=".true."
  else
    ATMOS_BOUNDARY_USE_QHYD=".false."
  fi

  #################################################
  #
  # make config files
  #
  #################################################

  mkdir -p ${OUTPUT_CONFIGDIR}/pp
  mkdir -p ${OUTPUT_CONFIGDIR}/init
  mkdir -p ${OUTPUT_CONFIGDIR}/run
  mkdir -p ${OUTPUT_CONFIGDIR}/sno

  source ${PP_CONF}

  if [ ${FILETYPE_ORG} = "NetCDF" ]; then
    PARENT_BASENAME="${BASENAME_ORG}"
  else
    PARENT_BASENAME=""
  fi

  source ${INIT_CONF}
  source ${RUN_CONF}
  source ${PARAM_CONF}
  source ${SNO_CONF}

  cat base.pp.conf \
      param.region.conf \
      param.admin.conf \
  > ${OUTPUT_CONFIGDIR}/pp/pp.d${FNUM}.conf

  cat base.init.conf \
      param.region.conf \
      param.admin.conf \
  > ${OUTPUT_CONFIGDIR}/init/init.d${FNUM}.conf

  cat base.run.conf \
      param.region.conf \
      param.admin.conf \
      param.physics.conf \
      param.history.conf \
  > ${OUTPUT_CONFIGDIR}/run/run.d${FNUM}.conf

  cat base.sno.vgridope.conf > ${OUTPUT_CONFIGDIR}/sno/sno.vgridope.d${FNUM}.conf
  cat base.sno.hgridope.conf > ${OUTPUT_CONFIGDIR}/sno/sno.hgridope.d${FNUM}.conf

  rm -f base.*.conf param.*.conf

  DNUM=`expr $DNUM + 1`
done


#################################################
#
# make launcher
#
#################################################

IFS="," eval 'LIST_PRC_DOMAINS="${PRC_DOMAINS[*]}"'
IFS="," eval 'LIST_PP_CONF_FILES="${PP_CONF_FILES[*]}"'
IFS="," eval 'LIST_INIT_CONF_FILES="${INIT_CONF_FILES[*]}"'
IFS="," eval 'LIST_RUN_CONF_FILES="${RUN_CONF_FILES[*]}"'
IFS="," eval 'LIST_SNO_CONF_FILES="${SNO_CONF_FILES[*]}"'

LIST_INIT_CONF_FILES=`echo ${LIST_INIT_CONF_FILES} | sed -e "s/,/\",\"/g"`
LIST_RUN_CONF_FILES=`echo ${LIST_RUN_CONF_FILES} | sed -e "s/,/\",\"/g"`

source ${LAUNCH_CONF}

mv -f base.init.launch.conf ${OUTPUT_CONFIGDIR}/init/init.launch.conf
mv -f base.run.launch.conf ${OUTPUT_CONFIGDIR}/run/run.launch.conf

#################################################
#
# set-up experimental environment
#
#################################################

if [ $DNUM -eq 1 ]; then
  INIT_CONF_FILE=init.d01.conf
  RUN_CONF_FILE=run.d01.conf
else
  INIT_CONF_FILE=init.launch.conf
  RUN_CONF_FILE=run.launch.conf
fi

if [ ${FILETYPE_ORG} = "NetCDF" ]; then
  eval 'NP=`expr ${PARENT_PRC_NUM_X} "*" ${PARENT_PRC_NUM_Y}`'
  DATPARAM="\"\""
  DATDISTS="\" [${NP} ../../data/${BASENAME_ORG} ${BASENAME_ORG}] \""
else
  DATPARAM="\" [../../data/${BASENAME_ORG} ${BASENAME_ORG}] \""
fi

source ${INPUT_CONFIGDIR}/mklink.sh

source ${INPUT_CONFIGDIR}/mkMakefile.sh

