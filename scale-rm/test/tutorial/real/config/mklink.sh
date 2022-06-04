#!/bin/bash

TOPDIR="../../../../../.."

RPATH_BIN="${TOPDIR}/bin"
RPATH_DAT="${TOPDIR}/data"

# symbolic link for database
ln -sf ${RPATH_DAT}/land/param.bucket.conf    ${OUTPUT_CONFIGDIR}/init/
ln -sf ${RPATH_DAT}/land/param.bucket.conf    ${OUTPUT_CONFIGDIR}/run/
ln -sf ${RPATH_DAT}/urban/param.kusaka01.dat  ${OUTPUT_CONFIGDIR}/run/
ln -sf ${RPATH_DAT}/rad/MIPAS                 ${OUTPUT_CONFIGDIR}/run/
ln -sf ${RPATH_DAT}/rad/PARAG.29              ${OUTPUT_CONFIGDIR}/run/
ln -sf ${RPATH_DAT}/rad/PARAPC.29             ${OUTPUT_CONFIGDIR}/run/
ln -sf ${RPATH_DAT}/rad/VARDATA.RM29          ${OUTPUT_CONFIGDIR}/run/
ln -sf ${RPATH_DAT}/rad/cira.nc               ${OUTPUT_CONFIGDIR}/run/

# symbolic link for executable binary
ln -sf ${RPATH_BIN}/scale-rm_pp   ${OUTPUT_CONFIGDIR}/pp/
ln -sf ${RPATH_BIN}/scale-rm_init ${OUTPUT_CONFIGDIR}/init/
ln -sf ${RPATH_BIN}/scale-rm      ${OUTPUT_CONFIGDIR}/run/
ln -sf ${RPATH_BIN}/sno           ${OUTPUT_CONFIGDIR}/sno/
