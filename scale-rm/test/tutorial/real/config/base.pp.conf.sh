#!/bin/bash

cat << EOF > base.pp.conf

#################################################
#
# model configuration: pp.conf only
#
#################################################


&PARAM_TOPOGRAPHY
 TOPOGRAPHY_OUT_BASENAME = "${TOPOGRAPHY_OUT_BASENAME}",
/

&PARAM_LANDUSE
 LANDUSE_OUT_BASENAME = "${LANDUSE_OUT_BASENAME}",
/

&PARAM_CONVERT
 CONVERT_TOPO    = .true.,
 CONVERT_LANDUSE = .true.,
/

&PARAM_CNVTOPO
 CNVTOPO_name                  = "${TOPOTYPE[$D]}",
 CNVTOPO_smooth_local          = ${SMOOTH_LOCAL[$D]},
 CNVTOPO_smooth_itelim         = ${SMOOTH_ITELIM},
 CNVTOPO_smooth_maxslope_ratio = ${MAXSLOPE_RATIO},
 CNVTOPO_copy_parent           = ${COPYTOPO[$D]},
/

&PARAM_COPYTOPO
 COPYTOPO_IN_BASENAME   = "${COPYTOPO_IN_BASENAME}",
 COPYTOPO_IN_FILETYPE   = "SCALE",
 COPYTOPO_ENTIRE_REGION = .false.,
 COPYTOPO_LINEAR_H      = .true.,
/

&PARAM_CNVTOPO_${TOPOTYPE[$D]}
 ${TOPOTYPE[$D]}_IN_DIR       = "${TOPO_IN_DIR}",
 ${TOPOTYPE[$D]}_IN_CATALOGUE = "${TOPO_IN_CATALOGUE}",
/

&PARAM_CNVLANDUSE
 CNVLANDUSE_name = "${LANDUSETYPE[$D]}",
 CNVLANDUSE_limit_urban_fraction = ${LIMIT_URBAN_FRACTION},
/

&PARAM_CNVLANDUSE_${LANDUSETYPE[$D]}
 ${LANDUSETYPE[$D]}_IN_DIR        = "${LANDUSE_IN_DIR}",
 ${LANDUSETYPE[$D]}_IN_CATALOGUE  = "${LANDUSE_IN_CATALOGUE}",
/

&PARAM_IO
 IO_LOG_BASENAME = "${PP_IO_LOG_BASENAME}",
/

EOF
