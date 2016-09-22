#!/bin/bash

cat << EOF > param.admin.conf

#################################################
#
# model configuration: administrator
#
#################################################

&PARAM_CONST
 CONST_THERMODYN_TYPE = "SIMPLE",
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
 ATMOS_VARS_CHECKRANGE = .true.,
/

&PARAM_ATMOS_REFSTATE
 ATMOS_REFSTATE_TYPE        = "INIT",
 ATMOS_REFSTATE_UPDATE_FLAG = .true.,
 ATMOS_REFSTATE_UPDATE_DT   = ${TIME_DT_REFSTATE},
/

&PARAM_ATMOS_BOUNDARY
 ATMOS_BOUNDARY_TYPE           = "REAL",
 ATMOS_BOUNDARY_IN_BASENAME    = "${ATMOS_BOUNDARY_IN_BASENAME}",
 ATMOS_BOUNDARY_START_DATE     = ${ATMOS_BOUNDARY_START_DATE},
 ATMOS_BOUNDARY_UPDATE_DT      = ${ATMOS_BOUNDARY_UPDATE_DT},
 ATMOS_BOUNDARY_USE_DENS       = .true.,
 ATMOS_BOUNDARY_USE_VELZ       = .true.,
 ATMOS_BOUNDARY_USE_QHYD       = ${ATMOS_BOUNDARY_USE_QHYD},
 ATMOS_BOUNDARY_VALUE_VELZ     = 0.0,
 ATMOS_BOUNDARY_ALPHAFACT_DENS = 1.0,
 ATMOS_BOUNDARY_LINEAR_H       = .false.,
 ATMOS_BOUNDARY_EXP_H          = 2.0,
/

&PARAM_ATMOS_DYN
 ATMOS_DYN_TINTEG_LARGE_TYPE          = "EULER",
 ATMOS_DYN_TINTEG_SHORT_TYPE          = "RK4",
 ATMOS_DYN_TINTEG_TRACER_TYPE         = "RK3WS2002",
 ATMOS_DYN_FVM_FLUX_TYPE              = "CD4",
 ATMOS_DYN_FVM_FLUX_TRACER_TYPE       = "UD3KOREN1993",
 ATMOS_DYN_NUMERICAL_DIFF_COEF        = 0.01,
 ATMOS_DYN_NUMERICAL_DIFF_COEF_TRACER = 0.0,
 ATMOS_DYN_enable_coriolis            = .true.,
 ATMOS_DYN_FLAG_FCT_TRACER            = .false.,
/

&PARAM_ATMOS_PHY_RD_MSTRN
 ATMOS_PHY_RD_MSTRN_KADD                  = 30,
 ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME   = "PARAG.29",
 ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME  = "PARAPC.29",
 ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME = "VARDATA.RM29",
/

&PARAM_ATMOS_PHY_RD_PROFILE
 ATMOS_PHY_RD_PROFILE_CIRA86_IN_FILENAME    = "cira.nc",
 ATMOS_PHY_RD_PROFILE_MIPAS2001_IN_BASENAME = "MIPAS",
/

#################################################
#
# model configuration: ocean
#
#################################################

&PARAM_OCEAN_VARS
 OCEAN_VARS_CHECKRANGE = .true.,
/

&PARAM_OCEAN_PHY_SLAB
 OCEAN_PHY_SLAB_DEPTH = 10.0,
/

#################################################
#
# model configuration: land
#
#################################################

&PARAM_LAND_VARS
 LAND_VARS_CHECKRANGE = .true.,
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
! URBAN_VARS_CHECKRANGE = .true.,
/

&PARAM_URBAN_PHY_SLC
 STRGR = 0.0,
 STRGB = 0.0,
 STRGG = 0.0,
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
 HISTORY_DEFAULT_TINTERVAL = ${TIME_DT_HISTORY_2D},
 HISTORY_DEFAULT_TUNIT     = "${TIME_DT_UNIT}",
 HISTORY_DEFAULT_TAVERAGE  = .false.,
 HISTORY_DEFAULT_DATATYPE  = "REAL4",
 HISTORY_DEFAULT_ZDIM      = "native",
 HISTORY_OUTPUT_STEP0      = .true.,
/

&PARAM_HIST
 HIST_BND = .false.,
/

EOF

if [ ${#HIST_ITEMS_SNAPSHOT_2D[*]} -ge 1 ]; then
  for VAR in ${HIST_ITEMS_SNAPSHOT_2D[*]}
  do
    echo "&HISTITEM item=\"${VAR}\" /" >> param.history.conf
  done
fi
if [ ${#HIST_ITEMS_SNAPSHOT_3D[*]} -ge 1 ]; then
  for VAR in ${HIST_ITEMS_SNAPSHOT_3D[*]}
  do
    echo "&HISTITEM item=\"${VAR}\", tinterval=${TIME_DT_HISTORY_3D} /" >> param.history.conf
  done
fi
if [ ${#HIST_ITEMS_AVERAGE_2D[*]} -ge 1 ]; then
  for VAR in ${HIST_ITEMS_AVERAGE_2D[*]}
  do
    echo "&HISTITEM item=\"${VAR}\", taverage=.true. /" >> param.history.conf
  done
fi
if [ ${#HIST_ITEMS_AVERAGE_3D[*]} -ge 1 ]; then
  for VAR in ${HIST_ITEMS_AVERAGE_3D[*]}
  do
    echo "&HISTITEM item=\"${VAR}\", taverage=.true., tinterval=${TIME_DT_HISTORY_3D} /" >> param.history.conf
  done
fi
