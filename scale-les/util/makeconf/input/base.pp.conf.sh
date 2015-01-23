#!/bin/bash

cat << EOF > base.pp.conf

#################################################
#
# model configuration: pp.conf only
#
#################################################

&PARAM_IO
 IO_LOG_BASENAME = "${PP_IO_LOG_BASENAME}",
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
 CNVTOPO_name                = "${TOPOTYPE[$D]}",
 CNVTOPO_SMOOTH_MAXSLOPE     = ${MAXSLOPE[$D]},
 CNVTOPO_SMOOTH_MAXSLOPE_BND = ${MAXSLOPE_BND[$D]},
 CNVTOPO_SMOOTH_itelim       = 10000,
/

&PARAM_CNVTOPO_${TOPOTYPE[$D]}
 ${TOPOTYPE[$D]}_IN_CATALOGUE = "${TOPO_IN_CATALOGUE}",
 ${TOPOTYPE[$D]}_IN_DIR       = "${TOPO_IN_DIR}",
/

&PARAM_CNVLANDUSE
 CNVLANDUSE_name = "${LANDUSETYPE[$D]}",
/

&PARAM_CNVLANDUSE_${LANDUSETYPE[$D]}
 ${LANDUSETYPE[$D]}_IN_CATALOGUE = "${LANDUSE_IN_CATALOGUE}",
 ${LANDUSETYPE[$D]}_IN_DIR       = "${LANDUSE_IN_DIR}",
 limit_urban_fraction            = ${LIMIT_URBAN_FRACTION},
/
EOF
