#!/bin/bash

cat << EOF > param.admin.conf

#################################################
#
# model configuration: administrator
#
#################################################

&PARAM_CONST
 CONST_THERMODYN_TYPE = "${CONST_THERMODYN_TYPE}",
/

&PARAM_TRACER
 TRACER_TYPE = "${ATMOS_PHY_MP_TYPE[$D]}",
/

&PARAM_ATMOS
 ATMOS_DYN_TYPE    = "${ATMOS_DYN_TYPE[$D]}",
 ATMOS_PHY_MP_TYPE = "${ATMOS_PHY_MP_TYPE[$D]}",
 ATMOS_PHY_RD_TYPE = "${ATMOS_PHY_RD_TYPE[$D]}",
 ATMOS_PHY_SF_TYPE = "${ATMOS_PHY_SF_TYPE[$D]}",
 ATMOS_PHY_TB_TYPE = "${ATMOS_PHY_TB_TYPE[$D]}",
/

&PARAM_OCEAN
 OCEAN_TYPE = "${OCEAN_TYPE[$D]}",
/

&PARAM_LAND
 LAND_TYPE = "${LAND_TYPE[$D]}",
/

&PARAM_URBAN
 URBAN_TYPE = "${URBAN_TYPE[$D]}",
/
EOF

cat << EOF > param.region.conf

#################################################
#
# model configuration: process
#
#################################################

&PARAM_PRC
 PRC_NUM_X      = ${PRC_NUM_X[$D]},
 PRC_NUM_Y      = ${PRC_NUM_Y[$D]},
 PRC_PERIODIC_X = .false.,
 PRC_PERIODIC_Y = .false.,
/

#################################################
#
# model configuration: region
#
#################################################

&PARAM_INDEX
 KMAX = ${KMAX[$D]},
 IMAX = ${IMAX[$D]},
 JMAX = ${JMAX[$D]},
/

&PARAM_LAND_INDEX
 LKMAX = ${LKMAX},
/

&PARAM_URBAN_INDEX
 UKMAX = ${UKMAX},
/

&PARAM_LAND_GRID
 LDZ = ${LIST_LDZ},
/

&PARAM_URBAN_GRID
 UDZ = ${LIST_UDZ},
/

&PARAM_GRID
 DX = ${DX[$D]},
 DY = ${DY[$D]},
 ${LINE_Z}
 BUFFER_DZ = ${BUFFER_DZ[$D]},
 BUFFER_DX = ${BUFFER_DX[$D]},
 BUFFER_DY = ${BUFFER_DY[$D]},
/

&PARAM_MAPPROJ
 MPRJ_basepoint_lon = ${MPRJ_BASEPOINT_LON},
 MPRJ_basepoint_lat = ${MPRJ_BASEPOINT_LAT},
 MPRJ_type          = "${MPRJ_TYPE}",
 MPRJ_LC_lat1       = ${MPRJ_LC_LAT1},
 MPRJ_LC_lat2       = ${MPRJ_LC_LAT2},
/
EOF

cat << EOF > param.physics.conf

#################################################
#
# model configuration: atmosphere
#
#################################################

&PARAM_ATMOS_VARS
 ATMOS_VARS_CHECKRANGE = .false.,
/

&PARAM_ATMOS_REFSTATE
 ATMOS_REFSTATE_TYPE        = "INIT",
 ATMOS_REFSTATE_UPDATE_FLAG = .true.,
 ATMOS_REFSTATE_UPDATE_DT   = ${TIME_DT_REFSTATE},
/

&PARAM_ATMOS_BOUNDARY
 ATMOS_BOUNDARY_TYPE        = "REAL",
 ATMOS_BOUNDARY_IN_BASENAME = "${ATMOS_BOUNDARY_IN_BASENAME}",
 ATMOS_BOUNDARY_UPDATE_DT   = ${ATMOS_BOUNDARY_UPDATE_DT},
 ATMOS_BOUNDARY_USE_VELZ    = .true.,
 ATMOS_BOUNDARY_USE_QHYD    = ${ATMOS_BOUNDARY_USE_QHYD},
 ATMOS_BOUNDARY_VALUE_VELZ  = 0.0D0,
 ATMOS_BOUNDARY_LINEAR_H    = .false.,
 ATMOS_BOUNDARY_EXP_H       = 2.d0,
/

&PARAM_ATMOS_DYN
 ATMOS_DYN_NUMERICAL_DIFF_COEF   = 1.D-2,
 ATMOS_DYN_NUMERICAL_DIFF_COEF_Q = 1.D-2,
 ATMOS_DYN_enable_coriolis       = .true.,
/

&PARAM_ATMOS_PHY_RD_MSTRN
 ATMOS_PHY_RD_MSTRN_KADD                  = 30,
 ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME   = "${RD_MSTRN_GASPARA_IN_FILENAME}",
 ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME  = "${RD_MSTRN_AEROPARA_IN_FILENAME}",
 ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME = "${RD_MSTRN_HYGROPARA_IN_FILENAME}",
 ATMOS_PHY_RD_MSTRN_NBAND                 = 29,
/

&PARAM_ATMOS_PHY_RD_PROFILE
 ATMOS_PHY_RD_PROFILE_TOA                   = 100.D0,
 ATMOS_PHY_RD_PROFILE_CIRA86_IN_FILENAME    = "${RD_PROFILE_CIRA86_IN_FILENAME}",
 ATMOS_PHY_RD_PROFILE_MIPAS2001_IN_BASENAME = "${RD_PROFILE_MIPAS2001_IN_BASENAME}",
/

#################################################
#
# model configuration: ocean
#
#################################################

&PARAM_OCEAN_VARS
 OCEAN_VARS_CHECKRANGE = .false.,
/

&PARAM_OCEAN_PHY_SLAB
 OCEAN_PHY_SLAB_DEPTH = 10.D0,
/

#################################################
#
# model configuration: land
#
#################################################

&PARAM_LAND_VARS
 LAND_VARS_CHECKRANGE = .false.,
/

&PARAM_LAND_PHY_SLAB
 LAND_PHY_UPDATE_BOTTOM_TEMP  = .false.,
 LAND_PHY_UPDATE_BOTTOM_WATER = .true.,
/

#################################################
#
# model configuration: urban
#
#################################################

&PARAM_URBAN_VARS
 URBAN_VARS_CHECKRANGE = .false.,
/

&PARAM_URBAN_PHY_SLC
 ZR         = 15.0D0,
 roof_width = 7.5D0,
 road_width = 22.5D0,
 AH         = 0.0D0,
 ALH        = 0.0D0,
 STRGR      = 0.24D0,
 STRGB      = 0.009D0,
 STRGG      = 0.24D0,
 AKSR       = 2.28D0,
 AKSB       = 2.28D0,
 AKSG       = 2.28D0,
 ALBR       = 0.20D0,
 ALBB       = 0.20D0,
 ALBG       = 0.20D0,
 EPSR       = 0.97D0,
 EPSB       = 0.97D0,
 EPSG       = 0.97D0,
 Z0R        = 0.005D0,
 Z0B        = 0.005D0,
 Z0G        = 0.005D0,
 CAPR       = 2.01D6,
 CAPB       = 2.01D6,
 CAPG       = 2.01D6,
/
EOF

cat << EOF > param.history.conf

#################################################
#
# model configuration: history
#
#################################################

&PARAM_HISTORY
 HISTORY_DEFAULT_BASENAME  = "${HISTORY_DEFAULT_BASENAME}",
 HISTORY_DEFAULT_TINTERVAL = ${TIME_DT_HISTORY},
 HISTORY_DEFAULT_TUNIT     = "${TIME_UNIT}",
 HISTORY_DEFAULT_TAVERAGE  = .false.,
 HISTORY_DEFAULT_DATATYPE  = "REAL4",
 HISTORY_DEFAULT_ZINTERP   = .false.,
/

&PARAM_HIST
 HIST_BND = .false.,
/

EOF

if [ ${#HIST_ITEMS_SNAPSHOT[*]} -ge 1 ]; then
  for VAR in ${HIST_ITEMS_SNAPSHOT[*]}
  do
    echo "&HISTITEM item=\"${VAR}\" /" >> param.history.conf
  done
fi
if [ ${#HIST_ITEMS_AVERAGE[*]} -ge 1 ]; then
  for VAR in ${HIST_ITEMS_AVERAGE[*]}
  do
    echo "&HISTITEM item=\"${VAR}\", taverage=.true. /" >> param.history.conf
  done
fi
