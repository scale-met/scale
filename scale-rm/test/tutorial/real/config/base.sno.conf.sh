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
 nprocs_x_out    = ${PRC_NUM_X[$D]},
 nprocs_y_out    = ${PRC_NUM_Y[$D]},
 vars            = '',
 output_single   = .true.,
/

&PARAM_SNOPLGIN_HGRIDOPE
 SNOPLGIN_hgridope_type      = 'LATLON',
 SNOPLGIN_hgridope_lat_start = 20.0,
 SNOPLGIN_hgridope_lat_end   = 50.0,
 SNOPLGIN_hgridope_dlat      = 0.2,
 SNOPLGIN_hgridope_lon_start = 120.0,
 SNOPLGIN_hgridope_lon_end   = 150.0,
 SNOPLGIN_hgridope_dlon      = 0.2,
/
EOF
