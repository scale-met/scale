#!/bin/bash

RPATH_RUN="./../run"
RPATH_TOOLS="./../../tools/FNL_output"
RPATH_BIN="./../../../bin"
RPATH_DATA="./../../../../data"
RPATH_NET2G="./../../../../../util/netcdf2grads_h"

# symbolic link for initial data
LINK_BOUNDARY_DT=`echo ${TIME_DT_BOUNDARY%%.*}`

fn=0
while [ ${fn} -lt ${NUMBER_OF_FILES} ]
do
  init_sec=`date --utc --date="${STARTDATE[0]}/${STARTDATE[1]}/${STARTDATE[2]} ${STARTDATE[3]}:${STARTDATE[4]}:${STARTDATE[5]}" +%s`
  unix_sec=`expr ${init_sec} + ${LINK_BOUNDARY_DT} \* ${fn}`

  year=`date -u -d "@${unix_sec}" +%Y`
  month=`date -u -d "@${unix_sec}" +%m`
  day=`date -u -d "@${unix_sec}" +%d`
  hour=`date -u -d "@${unix_sec}" +%H`
  min=`date -u -d "@${unix_sec}" +%M`
  sec=`date -u -d "@${unix_sec}" +%S`

  fmtd_fn=`printf "%05d" $fn`

  ln -sf ${RPATH_FNL_OUTPUT}/${year}${month}/FNLatm_${year}${month}${day}${hour}.grd  ${OUTPUT_CONFIGDIR}/init/FNLatm_${fmtd_fn}.grd
  ln -sf ${RPATH_FNL_OUTPUT}/${year}${month}/FNLsfc_${year}${month}${day}${hour}.grd  ${OUTPUT_CONFIGDIR}/init/FNLsfc_${fmtd_fn}.grd
  ln -sf ${RPATH_FNL_OUTPUT}/${year}${month}/FNLland_${year}${month}${day}${hour}.grd ${OUTPUT_CONFIGDIR}/init/FNLland_${fmtd_fn}.grd

  fn=`expr ${fn} + 1`
done

# symbolic link for database
ln -sf ${RPATH_DATA}/land/param.bucket.conf ${OUTPUT_CONFIGDIR}/init/
ln -sf ${RPATH_DATA}/land/param.bucket.conf ${OUTPUT_CONFIGDIR}/run/
ln -sf ${RPATH_DATA}/rad/MIPAS              ${OUTPUT_CONFIGDIR}/run/
ln -sf ${RPATH_DATA}/rad/PARAG.29           ${OUTPUT_CONFIGDIR}/run/
ln -sf ${RPATH_DATA}/rad/PARAPC.29          ${OUTPUT_CONFIGDIR}/run/
ln -sf ${RPATH_DATA}/rad/VARDATA.RM29       ${OUTPUT_CONFIGDIR}/run/
ln -sf ${RPATH_DATA}/rad/cira.nc            ${OUTPUT_CONFIGDIR}/run/

# symbolic link for net2g
DNUM=1
while [ ${DNUM} -le ${NUM_DOMAIN} ]
do
  ln -sf ${RPATH_RUN}/run.d`printf "%02d" $DNUM`.conf ${OUTPUT_CONFIGDIR}/net2g/
  DNUM=`expr ${DNUM} + 1`
done

# symbolic link for executable binary
ln -sf ${RPATH_BIN}/scale-rm_pp   ${OUTPUT_CONFIGDIR}/pp/
ln -sf ${RPATH_BIN}/scale-rm_init ${OUTPUT_CONFIGDIR}/init/
ln -sf ${RPATH_BIN}/scale-rm      ${OUTPUT_CONFIGDIR}/run/
ln -sf ${RPATH_NET2G}/net2g       ${OUTPUT_CONFIGDIR}/net2g/
