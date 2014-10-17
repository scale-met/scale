#!/bin/bash

cat << EOF > conf/base.init.conf

#################################################
#
# model configuration: init.conf only
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

&PARAM_RESTART
 RESTART_OUTPUT         = .true.,
 RESTART_OUT_BASENAME   = "${RESTART_OUT_BASENAME}",
/

&PARAM_TOPO
 TOPO_IN_BASENAME     = "${TOPO_IN_BASENAME}",
/

&PARAM_LANDUSE
 LANDUSE_IN_BASENAME  = "${LANDUSE_IN_BASENAME}",
/

&PARAM_MKINIT
 MKINIT_initname = "REAL",
/

&PARAM_MKINIT_REAL
 BASENAME_BOUNDARY   = "${BASENAME_BOUNDARY}",
 BASENAME_ORG        = "${BASENAME_ORG}",
 FILETYPE_ORG        = "WRF-ARW",
 NUMBER_OF_FILES     = 37,
 BOUNDARY_UPDATE_DT  = 600.D0,
 INTERP_SERC_DIV_NUM = 20,
 WRF_FILE_TYPE       = .true.,
/
EOF
