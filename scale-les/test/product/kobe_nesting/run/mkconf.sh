#!/bin/bash

#################################################
#
# SCALE-LES configuration file generator
#
#################################################

NUM_DOMAIN=1

INDIR="../input"
OUTDIR="../output"

RADDIR="../../../data/rad"
TOPODIR="/data2/scale_database/topo"
TOPOTYPE="GTOPO30"
LANDUSEDIR="/data2/scale_database/landuse"
LANDUSETYPE="GLCCv2"

BASENAME_ORG="${INDIR}/wrfout_d01_2011-09"

#################################################

mkdir -p ${INDIR}
mkdir -p ${OUTDIR}

NUM=1
while [ $NUM -le $NUM_DOMAIN ]
do
FNUM=`printf "%02d" $NUM`

#################################################
# 
# domain-separated config files
#
#################################################

if [ -f "conf/base.pp.d${FNUM}.conf.sh" ]; then
  BASE_PP="conf/base.pp.d${FNUM}.conf.sh"
else
  BASE_PP="conf/base.pp.conf.sh"
fi

if [ -f "conf/base.init.d${FNUM}.conf.sh" ]; then
  BASE_INIT="conf/base.init.d${FNUM}.conf.sh"
else
  BASE_INIT="conf/base.init.conf.sh"
fi

if [ -f "conf/base.run.d${FNUM}.conf.sh" ]; then
  BASE_RUN="conf/base.run.d${FNUM}.conf.sh"
else
  BASE_RUN="conf/base.run.conf.sh"
fi

if [ -f "conf/param.region.d${FNUM}.conf.sh" ]; then
  CONF_REGION="conf/param.region.d${FNUM}.conf.sh"
else
  CONF_REGION="conf/param.region.conf.sh"
fi

if [ -f "conf/param.admin.d${FNUM}.conf.sh" ]; then
  CONF_ADMIN="conf/param.admin.d${FNUM}.conf.sh"
else
  CONF_ADMIN="conf/param.admin.conf.sh"
fi

if [ -f "conf/param.physics.d${FNUM}.conf.sh" ]; then
  CONF_PHYSICS="conf/param.physics.d${FNUM}.conf.sh"
else
  CONF_PHYSICS="conf/param.physics.conf.sh"
fi

if [ -f "conf/param.bucket.d${FNUM}.conf.sh" ]; then
  CONF_BUCKET="conf/param.bucket.d${FNUM}.conf.sh"
else
  CONF_BUCKET="conf/param.bucket.conf.sh"
fi

if [ -f "conf/param.history.d${FNUM}.conf.sh" ]; then
  CONF_HISTORY="conf/param.history.d${FNUM}.conf.sh"
else
  CONF_HISTORY="conf/param.history.conf.sh"
fi

if [ -f "conf/param.monitor.d${FNUM}.conf.sh" ]; then
  CONF_MONITOR="conf/param.monitor.d${FNUM}.conf.sh"
else
  CONF_MONITOR="conf/param.monitor.conf.sh"
fi

#################################################
# 
# make pp.conf
#
#################################################

IO_LOG_BASENAME="${INDIR}/pp_LOG_d${FNUM}"
TOPO_OUT_BASENAME="${INDIR}/topo_d${FNUM}"
LANDUSE_OUT_BASENAME="${INDIR}/landuse_d${FNUM}"

source ${BASE_PP}
source ${CONF_REGION}
source ${CONF_ADMIN}

cat conf/base.pp.conf \
    conf/param.region.conf \
    conf/param.admin.conf \
> pp.d${FNUM}.conf

rm -f conf/*.conf

#################################################
# 
# make init.conf
#
#################################################

IO_LOG_BASENAME="${INDIR}/init_LOG_d${FNUM}"
RESTART_OUT_BASENAME="${INDIR}/init_d${FNUM}"
TOPO_IN_BASENAME="${INDIR}/topo_d${FNUM}"
LANDUSE_IN_BASENAME="${INDIR}/landuse_d${FNUM}"
BASENAME_BOUNDARY="${INDIR}/boundary_d${FNUM}"

source ${BASE_INIT}
source ${CONF_REGION}
source ${CONF_ADMIN}

cat conf/base.init.conf \
    conf/param.region.conf \
    conf/param.admin.conf \
> init.d${FNUM}.conf

rm -f conf/*.conf

#################################################
# 
# make run.conf
#
#################################################

IO_LOG_BASENAME="${OUTDIR}/LOG_d${FNUM}"
RESTART_IN_BASENAME="${INDIR}/init_d${FNUM}_00000000000.000"
TOPO_IN_BASENAME="${INDIR}/topo_d${FNUM}"
LANDUSE_IN_BASENAME="${INDIR}/landuse_d${FNUM}"
ATMOS_BOUNDARY_IN_BASENAME="${INDIR}/boundary_d${FNUM}"
HISTORY_DEFAULT_BASENAME="${OUTDIR}/history_d${FNUM}"
MONITOR_OUT_BASENAME="${OUTDIR}/monitor_d${FNUM}"

source ${BASE_RUN}
source ${CONF_REGION}
source ${CONF_ADMIN}
source ${CONF_PHYSICS}
source ${CONF_BUCKET}
source ${CONF_HISTORY}
source ${CONF_MONITOR}

cat conf/base.run.conf \
    conf/param.region.conf \
    conf/param.admin.conf \
    conf/param.physics.conf \
    conf/param.bucket.conf \
    conf/param.history.conf \
    conf/param.monitor.conf \
> run.d${FNUM}.conf

rm -f conf/*.conf

#################################################

NUM=`expr $NUM + 1`
done
