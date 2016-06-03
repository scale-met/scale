#!/bin/bash

cat << EOF > conf/param.physics.conf

#################################################
#
# model configuration: atmosphere
#
#################################################

&PARAM_ATMOS_VARS
 ATMOS_VARS_CHECKRANGE = .true.,
/

&PARAM_ATMOS_REFSTATE
 ATMOS_REFSTATE_TYPE = "INIT",
/

&PARAM_ATMOS_BOUNDARY
 ATMOS_BOUNDARY_TYPE        = "REAL",
 ATMOS_BOUNDARY_IN_BASENAME = "${ATMOS_BOUNDARY_IN_BASENAME}",
 ATMOS_BOUNDARY_USE_VELZ    = .true.,
 ATMOS_BOUNDARY_USE_QHYD    = .true.,
 ATMOS_BOUNDARY_VALUE_VELZ  = 0.0D0,
 ATMOS_BOUNDARY_UPDATE_DT   = 600.0D0,
/

&PARAM_ATMOS_DYN
 ATMOS_DYN_NUMERICAL_DIFF_COEF        = 1.D-2,
 ATMOS_DYN_NUMERICAL_DIFF_COEF_TRACER = 1.D-4,
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
 DEBUG                                      = .true.,
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
 OCEAN_PHY_SLAB_DEPTH = 10.D0,
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

#################################################
#
# model configuration: coupler
#
#################################################

&PARAM_URBAN_PHY_SLC
 STRGR = 0.24D0,
 STRGB = 0.009D0,
 STRGG = 0.24D0,
/
EOF
