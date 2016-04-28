#!/bin/bash

# path to history file from SCALE
history_dir=../../../../case_real/check_awajishima/run_WRF/output_new/

# path to output directory
output_dir=./

# time step and initial time of history file (These are used for ctl file)
time_step=2      # the number of history data
time_interval="5mn"  # (mn, hr, mon, day) 
start_time="00:05Z07Sep2011"

cat > namelist.in <<EOF
&info
timestep=${time_step}
idir="${history_dir}"
odir="${output_dir}"
vcount=13
&end
&vari
vname="DENS","RHOT","T","PRES","QV","QC","QR","QI","QS","QG","SFC_TEMP","SFC_ALB_LW","SFC_ALB_SW",
&end
&grads
delt="${time_interval}"
stime="${start_time}"
&end
EOF

rm -f convine
ifort -O3 -xHOST -convert big_endian -assume byterecl -I${NETCDF4}/include -L${NETCDF4}/lib -lnetcdff -lnetcdf make_grads_file.f90 -o convine
./convine

