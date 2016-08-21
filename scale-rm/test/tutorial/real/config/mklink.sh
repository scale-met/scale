#!/bin/bash

RPATH_BIN="./../../../bin"
RPATH_DATA="./../../../../data"
RPATH_NET2G="./../../../../../util/netcdf2grads_h"

# symbolic link for database
ln -sf ${RPATH_DATA}/land/param.bucket.conf ${OUTPUT_CONFIGDIR}/init/
ln -sf ${RPATH_DATA}/land/param.bucket.conf ${OUTPUT_CONFIGDIR}/run/
ln -sf ${RPATH_DATA}/rad/MIPAS              ${OUTPUT_CONFIGDIR}/run/
ln -sf ${RPATH_DATA}/rad/PARAG.29           ${OUTPUT_CONFIGDIR}/run/
ln -sf ${RPATH_DATA}/rad/PARAPC.29          ${OUTPUT_CONFIGDIR}/run/
ln -sf ${RPATH_DATA}/rad/VARDATA.RM29       ${OUTPUT_CONFIGDIR}/run/
ln -sf ${RPATH_DATA}/rad/cira.nc            ${OUTPUT_CONFIGDIR}/run/

# symbolic link for executable binary
ln -sf ${RPATH_BIN}/scale-rm_pp   ${OUTPUT_CONFIGDIR}/pp/
ln -sf ${RPATH_BIN}/scale-rm_init ${OUTPUT_CONFIGDIR}/init/
ln -sf ${RPATH_BIN}/scale-rm      ${OUTPUT_CONFIGDIR}/run/
ln -sf ${RPATH_NET2G}/net2g       ${OUTPUT_CONFIGDIR}/net2g/
