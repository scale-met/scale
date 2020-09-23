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

&PARAM_ATMOS
 ATMOS_DYN_TYPE    = "${ATMOS_DYN_TYPE[$D]}",
 ATMOS_PHY_CP_TYPE = "${ATMOS_PHY_CP_TYPE[$D]}",
 ATMOS_PHY_MP_TYPE = "${ATMOS_PHY_MP_TYPE[$D]}",
 ATMOS_PHY_RD_TYPE = "${ATMOS_PHY_RD_TYPE[$D]}",
 ATMOS_PHY_SF_TYPE = "${ATMOS_PHY_SF_TYPE[$D]}",
 ATMOS_PHY_TB_TYPE = "${ATMOS_PHY_TB_TYPE[$D]}",
 ATMOS_PHY_BL_TYPE = "${ATMOS_PHY_BL_TYPE[$D]}",
/

&PARAM_OCEAN
 OCEAN_DYN_TYPE = "${OCEAN_DYN_TYPE[$D]}",
/

&PARAM_LAND
 LAND_DYN_TYPE = "${LAND_DYN_TYPE[$D]}",
 LAND_SFC_TYPE = "${LAND_SFC_TYPE[$D]}",
/

&PARAM_URBAN
 URBAN_DYN_TYPE = "${URBAN_DYN_TYPE[$D]}",
/
EOF

cat << EOF > param.region.conf

#################################################
#
# model configuration: process
#
#################################################

&PARAM_PRC_CARTESC
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

&PARAM_ATMOS_GRID_CARTESC_INDEX
 KMAX  = ${KMAX[$D]},
 IMAXG = ${IMAXG[$D]},
 JMAXG = ${JMAXG[$D]},
/

&PARAM_OCEAN_GRID_CARTESC_INDEX
 OKMAX = ${OKMAX},
/

&PARAM_LAND_GRID_CARTESC_INDEX
 LKMAX = ${LKMAX},
/

&PARAM_URBAN_GRID_CARTESC_INDEX
 UKMAX = ${UKMAX},
/

&PARAM_OCEAN_GRID_CARTESC
 ODZ = ${LIST_ODZ},
/

&PARAM_LAND_GRID_CARTESC
 LDZ = ${LIST_LDZ},
/

&PARAM_URBAN_GRID_CARTESC
 UDZ = ${LIST_UDZ},
/

&PARAM_ATMOS_GRID_CARTESC
 DX = ${DX[$D]},
 DY = ${DY[$D]},
 ${LINE_Z}
 BUFFER_DZ = ${BUFFER_DZ[$D]},
 BUFFER_DX = ${BUFFER_DX[$D]},
 BUFFER_DY = ${BUFFER_DY[$D]},
/

&PARAM_MAPPROJECTION
 MAPPROJECTION_basepoint_lon = ${MAPPROJECTION_BASEPOINT_LON},
 MAPPROJECTION_basepoint_lat = ${MAPPROJECTION_BASEPOINT_LAT},
 MAPPROJECTION_type          = "${MAPPROJECTION_TYPE}",
 MAPPROJECTION_LC_lat1       = ${MAPPROJECTION_LC_LAT1},
 MAPPROJECTION_LC_lat2       = ${MAPPROJECTION_LC_LAT2},
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

&PARAM_ATMOS_HYDROSTATIC
 HYDROSTATIC_barometric_law_mslp_kref = 2,
/

&PARAM_ATMOS_REFSTATE
 ATMOS_REFSTATE_TYPE      = "INIT",
 ATMOS_REFSTATE_UPDATE_DT = ${TIME_DT_REFSTATE},
/

&PARAM_ATMOS_BOUNDARY
 ATMOS_BOUNDARY_TYPE           = "REAL",
 ATMOS_BOUNDARY_IN_BASENAME    = "${ATMOS_BOUNDARY_IN_BASENAME}",
 ATMOS_BOUNDARY_START_DATE     = ${ATMOS_BOUNDARY_START_DATE},
 ATMOS_BOUNDARY_UPDATE_DT      = ${ATMOS_BOUNDARY_UPDATE_DT},
 ATMOS_BOUNDARY_USE_QHYD       = ${ATMOS_BOUNDARY_USE_QHYD},
 ATMOS_BOUNDARY_LINEAR_H       = .false.,
 ATMOS_BOUNDARY_EXP_H          = 2.0,
 ATMOS_GRID_NUDGING_uv         = .false.,
 ATMOS_GRID_NUDGING_pt         = .false.,
 ATMOS_GRID_NUDGING_qv         = .false.,
 ATMOS_GRID_NUDGING_tau        = 864000.,
/

&PARAM_ATMOS_DYN
 ATMOS_DYN_TINTEG_LARGE_TYPE          = "EULER",
 ATMOS_DYN_TINTEG_SHORT_TYPE          = "RK4",
 ATMOS_DYN_TINTEG_TRACER_TYPE         = "RK3WS2002",
 ATMOS_DYN_FVM_FLUX_TYPE              = "UD3",
 ATMOS_DYN_FVM_FLUX_TRACER_TYPE       = "UD3KOREN1993",
 ATMOS_DYN_NUMERICAL_DIFF_COEF        = 0.0,
 ATMOS_DYN_NUMERICAL_DIFF_COEF_TRACER = 0.0,
 ATMOS_DYN_FLAG_FCT_TRACER            = .false.,
 ATMOS_DYN_WDAMP_HEIGHT               = 15.D3,
/

&PARAM_CORIOLIS
 CORIOLIS_type = "SPHERE",
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

EOF

if [ ${ATMOS_PHY_BL_TYPE[$D]} = "MYNN" ]; then
if [ ${ATMOS_PHY_TB_TYPE[$D]} = "SMAGORINSKY" ]; then
  cat <<EOF >> param.physics.conf
&PARAM_ATMOS_PHY_TB_SMG
 ATMOS_PHY_TB_SMG_horizontal = .true.,
/

EOF
fi
fi

cat <<EOF >> param.physics.conf
#################################################
#
# model configuration: ocean
#
#################################################

&PARAM_OCEAN_VARS
 OCEAN_VARS_CHECKRANGE = .true.,
/

#################################################
#
# model configuration: land
#
#################################################

&PARAM_LAND_VARS
 LAND_VARS_CHECKRANGE = .true.,
/

&PARAM_LAND_DYN_BUCKET
 LAND_DYN_BUCKET_UPDATE_BOTTOM_TEMP  = .false.,
 LAND_DYN_BUCKET_UPDATE_BOTTOM_WATER = .true.,
/

#################################################
#
# model configuration: urban
#
#################################################

&PARAM_URBAN_VARS
 URBAN_VARS_CHECKRANGE = .true.,
/

&PARAM_URBAN_DYN_KUSAKA01
 URBAN_DYN_KUSAKA01_PARAM_IN_FILENAME = 'param.kusaka01.dat',
/
EOF

cat << EOF > param.history.conf

#################################################
#
# model configuration: history
#
#################################################

&PARAM_FILE_HISTORY
 FILE_HISTORY_DEFAULT_BASENAME  = "${FILE_HISTORY_DEFAULT_BASENAME}",
 FILE_HISTORY_DEFAULT_TINTERVAL = ${TIME_DT_HISTORY_2D},
 FILE_HISTORY_DEFAULT_TUNIT     = "${TIME_DT_UNIT}",
 FILE_HISTORY_DEFAULT_TAVERAGE  = .false.,
 FILE_HISTORY_DEFAULT_DATATYPE  = "REAL4",
 FILE_HISTORY_DEFAULT_ZCOORD    = "model",
 FILE_HISTORY_OUTPUT_STEP0      = .true.,
/

&PARAM_FILE_HISTORY_CARTESC
 FILE_HISTORY_CARTESC_BOUNDARY = .false.,
/

EOF

if [ ${#HIST_ITEMS_SNAPSHOT_2D[*]} -ge 1 ]; then
  for VAR in ${HIST_ITEMS_SNAPSHOT_2D[*]}
  do
    echo "&HISTORY_ITEM name=\"${VAR}\" /" >> param.history.conf
  done
fi
if [ ${#HIST_ITEMS_SNAPSHOT_3D[*]} -ge 1 ]; then
  for VAR in ${HIST_ITEMS_SNAPSHOT_3D[*]}
  do
    echo "&HISTORY_ITEM name=\"${VAR}\", tinterval=${TIME_DT_HISTORY_3D} /" >> param.history.conf
  done
fi
if [ ${#HIST_ITEMS_AVERAGE_2D[*]} -ge 1 ]; then
  for VAR in ${HIST_ITEMS_AVERAGE_2D[*]}
  do
    echo "&HISTORY_ITEM name=\"${VAR}\", taverage=.true. /" >> param.history.conf
  done
fi
if [ ${#HIST_ITEMS_AVERAGE_3D[*]} -ge 1 ]; then
  for VAR in ${HIST_ITEMS_AVERAGE_3D[*]}
  do
    echo "&HISTORY_ITEM name=\"${VAR}\", taverage=.true., tinterval=${TIME_DT_HISTORY_3D} /" >> param.history.conf
  done
fi
