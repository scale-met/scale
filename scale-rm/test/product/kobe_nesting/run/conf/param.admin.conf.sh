#!/bin/bash

cat << EOF > conf/param.admin.conf

#################################################
#
# model configuration: administrator
#
#################################################

&PARAM_TRACER
 TRACER_TYPE = "TOMITA08",
/

&PARAM_ATMOS
 ATMOS_DYN_TYPE    = "HEVI",
 ATMOS_PHY_MP_TYPE = "TOMITA08",
 ATMOS_PHY_RD_TYPE = "MSTRNX",
 ATMOS_PHY_SF_TYPE = "COUPLE",
 ATMOS_PHY_TB_TYPE = "SMAGORINSKY",
/

&PARAM_OCEAN
 OCEAN_TYPE = "SLAB",
/

&PARAM_LAND
 LAND_TYPE = "SLAB",
/

&PARAM_URBAN
 URBAN_TYPE = "SLC",
/
EOF
