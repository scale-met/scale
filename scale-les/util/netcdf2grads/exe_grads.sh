#!/bin/bash

output_grads_dir="./grads"
if [ ! -d ${output_grads_dir} ] ; then   mkdir -p ${output_grads_dir} ; fi

cat > namelist.in <<EOF
&info
timestep=6
idir="../../test/case_real/check_awajishima/run_WRF"
odir="${output_grads_dir}"
vcount=1
&end
&vari
vname="SHFLX",
&end
&grads
delt="5mn"
stime="00:05Z07sep2011"
&end
EOF

rm -f convine
ifort -O3 -xHOST -convert big_endian -assume byterecl -I${NETCDF4}/include -L${NETCDF4}/lib -lnetcdff -lnetcdf make_grads_file.f90 -o convine
./convine -s

#---------------------------------------------------------------------
# timestep : the number of data
# idir     : path to the directory of history.pe*.nc
# odir     : path to the directory of grads files (output)
# vcount   : the number of variable you want to convert to grads format
# vname    : the name of variable you want to convert
# delt     : time interval of data, used for grads ctl file. The format is 1 or 2 digits and unit(yr mo dy hr mn sc).
# stime    : initial time for grads ctl file. The format is "hh:mmZddmmmyyyy".
#---------------------------------------------------------------------

#vcount=16
#vname="DENS","QHYD","RH","QC","QR","QI","QS","QG","T","PRES","SFC_ALB_LW","SFC_ALB_SW","PRES","RADFLUX_SWDN","SFLX_SW_dn","TOAFLX_SW_dn",
#vcount=26
#vname="DENS","MOMZ","MOMX","MOMY","RHOT","QV","SFC_TEMP","DENS_t_advcv","MOMZ_t_advcv","MOMX_t_advcv","MOMY_t_advcv","RHOT_t_advcv","DENS_t_advch","MOMZ_t_advch","MOMX_t_advch","MOMY_t_advch","RHOT_t_advch","MOMZ_t_pg","MOMX_t_pg","MOMY_t_pg","MOMZ_t_ddiv","MOMX_t_ddiv","MOMY_t_ddiv","MOMZ_t_cf","MOMX_t_cf","MOMY_t_cf",
#vcount=24
#vname="TC_URB","QC_URB","SHFLX_urb","LHFLX_urb","GHFLX_urb","Rngrd_URB","TR_URB","TB_URB","TG_URB","SHR_URB","SHB_URB","SHG_URB","LHR_URB","LHB_URB","LHG_URB","RnR_URB","RnB_URB","RnG_URB","TRL_URB","TBL_URB","TGL_URB","GHR_URB","GHB_URB","GHG_URB",

