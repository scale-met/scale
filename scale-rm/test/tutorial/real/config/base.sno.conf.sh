#!/bin/bash

cat << EOF > base.sno.vgridope.conf

#################################################
#
# post process configuration: sno
#
#################################################

&PARAM_SNO
 basename_in     = '${RUNDIR}/history_d${FNUM}',
 dirpath_out     = '.',
 basename_out    = 'merged-p_history_d${FNUM}',
 nprocs_x_out    = ${PRC_NUM_X[$D]},
 nprocs_y_out    = ${PRC_NUM_Y[$D]},
 vars            = '',
/

&PARAM_SNOPLGIN_VGRIDOPE
 SNOPLGIN_vgridope_type     = 'PLEV',
 SNOPLGIN_vgridope_lev_num  = 3,
 SNOPLGIN_vgridope_lev_data = 850.e+2, 500.e+2, 200.e+2,
/
EOF

pi=3.14159265359
buffer=1.1
earth_radius=6371000 # meter

dlon=$(bc -l <<< "${DX[$D]} / (c($MAPPROJECTION_BASEPOINT_LAT / (180 / $pi)) * $earth_radius) * (180 / $pi)")
dlat=$(bc -l <<< "${DY[$D]} / $earth_radius * (180 / $pi)")

half_width=$(bc -l <<< "$buffer * (${IMAXG[$D]} / 2) * $dlon")
half_height=$(bc -l <<< "$buffer * (${JMAXG[$D]} / 2) * $dlat")

cat << EOF > base.sno.hgridope.conf

#################################################
#
# post process configuration: sno
#
#################################################

&PARAM_SNO
 basename_in     = 'merged-p_history_d${FNUM}',
 dirpath_out     = '.',
 basename_out    = 'merged-h_history_d${FNUM}',
 nprocs_x_out    = 1,
 nprocs_y_out    = 1,
 vars            = '',
 output_single   = .true.,
/

&PARAM_SNOPLGIN_HGRIDOPE
 SNOPLGIN_hgridope_type      = 'LATLON',
 SNOPLGIN_hgridope_lat_start = $(bc -l <<< "$MAPPROJECTION_BASEPOINT_LAT - $half_height")
 SNOPLGIN_hgridope_lat_end   = $(bc -l <<< "$MAPPROJECTION_BASEPOINT_LAT + $half_height")
 SNOPLGIN_hgridope_dlat      = $dlat
 SNOPLGIN_hgridope_lon_start = $(bc -l <<< "$MAPPROJECTION_BASEPOINT_LON - $half_width")
 SNOPLGIN_hgridope_lon_end   = $(bc -l <<< "$MAPPROJECTION_BASEPOINT_LON + $half_width")
 SNOPLGIN_hgridope_dlon      = $dlon
/
EOF
