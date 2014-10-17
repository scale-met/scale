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
 CNVTOPO_name = "DEM50M",
/

&PARAM_CNVTOPO_DEM50M
 DEM50M_IN_CATALOGUE = "DEM50m_catalogue.txt",
 DEM50M_IN_DIR       = "${TOPODIR}",
/

&PARAM_CNVLANDUSE
 CNVLANDUSE_name = "LU100M",
/

&PARAM_CNVLANDUSE_LU100M
 LU100M_IN_CATALOGUE = "LU100m_catalogue.txt",
 LU100M_IN_DIR       = "${LANDUSEDIR}",
/
EOF
