#!/bin/bash

cat << EOF > conf/base.run.conf

#################################################
#
# model configuration: run.conf only
#
#################################################

&PARAM_IO
 IO_LOG_BASENAME = "${IO_LOG_BASENAME}",
 IO_LOG_ALLNODE = .true.,
/

&PARAM_TIME
 TIME_STARTDATE             = 2011, 9,19, 6, 0, 0,
 TIME_STARTMS               = 0.D0,
 TIME_DURATION              = 6.0D0,
 TIME_DURATION_UNIT         = "HOUR",
 TIME_DT                    = 3.D0,
 TIME_DT_UNIT               = "SEC",
 TIME_DT_ATMOS_DYN          = 0.6D0,
 TIME_DT_ATMOS_DYN_UNIT     = "SEC",
 TIME_DT_ATMOS_PHY_RD       = 300.D0,
 TIME_DT_ATMOS_PHY_RD_UNIT  = "SEC",
 TIME_DT_OCEAN              = 30.D0,
 TIME_DT_OCEAN_UNIT         = "SEC",
 TIME_DT_LAND               = 30.D0,
 TIME_DT_LAND_UNIT          = "SEC",
 TIME_DT_URBAN              = 30.D0,
 TIME_DT_URBAN_UNIT         = "SEC",
/

&PARAM_STATISTICS
 STATISTICS_checktotal     = .false.,
 STATISTICS_use_globalcomm = .true.,
/

&PARAM_RESTART
 RESTART_OUTPUT      = .false.,
 RESTART_IN_BASENAME = "${RESTART_IN_BASENAME}",
/

&PARAM_TOPO
 TOPO_IN_BASENAME = "${TOPO_IN_BASENAME}",
/

&PARAM_LANDUSE
 LANDUSE_IN_BASENAME = "${LANDUSE_IN_BASENAME}",
/
EOF
