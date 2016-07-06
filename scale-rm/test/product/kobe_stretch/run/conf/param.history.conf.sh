#!/bin/bash

cat << EOF > conf/param.history.conf

#################################################
#
# model configuration: history
#
#################################################

&PARAM_HISTORY
 HISTORY_DEFAULT_BASENAME  = "${HISTORY_DEFAULT_BASENAME}",
 HISTORY_DEFAULT_TINTERVAL = 300.D0,
 HISTORY_DEFAULT_TUNIT     = "SEC",
 HISTORY_DEFAULT_TAVERAGE  = .false.,
 HISTORY_DEFAULT_DATATYPE  = "REAL4",
 HISTORY_DEFAULT_ZINTERP   = .true.,
/

&HISTITEM item="DENS" /
&HISTITEM item="MOMZ" /
&HISTITEM item="MOMX" /
&HISTITEM item="MOMY" /
&HISTITEM item="RHOT" /

&HISTITEM item="QV"   /
&HISTITEM item="QC"   /
&HISTITEM item="QR"   /
&HISTITEM item="QI"   /
&HISTITEM item="QS"   /
&HISTITEM item="QG"   /
&HISTITEM item="QHYD" /
&HISTITEM item="QLIQ" /
&HISTITEM item="QICE" /

&HISTITEM item="T"    /
&HISTITEM item="PRES" /
&HISTITEM item="U"    /
&HISTITEM item="V"    /
&HISTITEM item="W"    /
&HISTITEM item="PT"   /
&HISTITEM item="RH"   /
&HISTITEM item="RHL"  /
&HISTITEM item="RHI"  /

&HISTITEM item="PREC" /

&HISTITEM item="SFC_PRES"   /
&HISTITEM item="SFC_TEMP"   /
&HISTITEM item="SFC_ALB_SW" /
&HISTITEM item="SFC_ALB_LW" /
&HISTITEM item="SFC_Z0M"    /

&HISTITEM item="U10" /
&HISTITEM item="V10" /
&HISTITEM item="T2"  /
&HISTITEM item="Q2"  /

&HISTITEM item="LHFLX" /
&HISTITEM item="SHFLX" /
&HISTITEM item="GHFLX" /

&HISTITEM item="RADFLUX_SWDN" /
&HISTITEM item="SFLX_LW_up"   /
&HISTITEM item="SFLX_LW_dn"   /
&HISTITEM item="SFLX_SW_up"   /
&HISTITEM item="SFLX_SW_dn"   /

&HISTITEM item="OCEAN_TEMP"     /
&HISTITEM item="OCEAN_SFC_TEMP" /
&HISTITEM item="OCEAN_ALB_SW"   /
&HISTITEM item="OCEAN_ALB_LW"   /
&HISTITEM item="OCEAN_SFC_Z0M"  /
&HISTITEM item="OCEAN_SFC_Z0H"  /
&HISTITEM item="OCEAN_SFC_Z0E"  /

&HISTITEM item="LAND_TEMP"     /
&HISTITEM item="LAND_WATER"    /
&HISTITEM item="LAND_SFC_TEMP" /
&HISTITEM item="LAND_ALB_SW"   /
&HISTITEM item="LAND_ALB_LW"   /

&HISTITEM item="URBAN_TC"       /
&HISTITEM item="URBAN_SFC_TEMP" /
EOF
