#!/bin/bash

cat << EOF > conf/base.pp.conf

#################################################
#
# model configuration: pp.conf only
#
#################################################

&PARAM_IO
 IO_LOG_BASENAME = "${IO_LOG_BASENAME}",
 IO_LOG_ALLNODE  = .true.,
/

&PARAM_TIME
 TIME_STARTDATE = 0000, 1, 1, 0, 0, 0,
 TIME_STARTMS   = 0.D0,
/

&PARAM_STATISTICS
 STATISTICS_checktotal     = .true.,
 STATISTICS_use_globalcomm = .true.,
/

&PARAM_TOPO
 TOPO_OUT_BASENAME = "${TOPO_OUT_BASENAME}",
/

&PARAM_LANDUSE
 LANDUSE_OUT_BASENAME = "${LANDUSE_OUT_BASENAME}",
/

&PARAM_CONVERT
 CONVERT_TOPO    = .true.,
 CONVERT_LANDUSE = .true.,
/

&PARAM_CNVTOPO
 CNVTOPO_name = "${TOPOTYPE}",
/

&PARAM_CNVTOPO_${TOPOTYPE}
 ${TOPOTYPE}_IN_CATALOGUE = "${TOPO_IN_CATALOGUE}",
 ${TOPOTYPE}_IN_DIR       = "${TOPO_IN_DIR}",
/

&PARAM_CNVLANDUSE
 CNVLANDUSE_name = "${LANDUSETYPE}",
/

&PARAM_CNVLANDUSE_${LANDUSETYPE}
 ${LANDUSETYPE}_IN_CATALOGUE = "${LANDUSE_IN_CATALOGUE}",
 ${LANDUSETYPE}_IN_DIR       = "${LANDUSE_IN_DIR}",
/
EOF
