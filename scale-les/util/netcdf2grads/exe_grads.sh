#!/bin/bash
#
case1="19990516-19990522"
ver="20150120.fd705f8"
itime=115
stime="06:00z16may1999"
#----------------------------
for icase in $case1 ; do  #$case1 $case2 $case3 $case4
#
for ielem in DENS  ; do
#
for inest in 2 3 ; do
inn=`printf "%2.2d" $inest`
#
output_conf_dir="../${ver}/${icase}/run.d${inn}.conf"
output_history_dir="../${ver}/${icase}/"
output_grads_dir="../grads/${ver}/${icase}/"
#----------------------------

if [ ! -d ${output_grads_dir} ] ; then   mkdir -p ${output_grads_dir} ; fi

cat > namelist.in <<EOF
&info
timestep=${itime}
inest=${inest}
conffile="${output_conf_dir}"
idir="${output_history_dir}"
odir="${output_grads_dir}"
vcount=1
&end
&vari
vname="${ielem}",
&end
&grads
delt="10mn"
stime="${stime}"
&end
EOF

rm -f convine
ifort -O3 -xHOST -convert big_endian -assume byterecl -I${NETCDF4}/include -L${NETCDF4}/lib -lnetcdff -lnetcdf make_grads_file.f90 -o convine
./convine -s

done # inest

done # ielem

done # icase

#---------------------------------------------------------------------
# timestep : the number of data
# inest    : domain number (used for the name of ctl file)
# conffile : path to the configure file, run.d0X.conf
# idir     : path to the directory of history.peXXXXXX.nc
# odir     : path to the directory of grads files (output)
# vcount   : the number of variables you want to convert to grads format
# vname    : the name of variables you want to convert
# delt     : time interval of data, used for grads ctl file. The format is 1 or 2 digits and unit(yr mo dy hr mn sc).
# stime    : initial time for grads ctl file. The format is "hh:mmZddmmmyyyy".
#---------------------------------------------------------------------

#"PRES","SFC_PRES","DENS","T","SFC_TEMP","PT","U","V","PREC","MSLP","RAIN","SNOW","LHFLX",
#vcount=16
#vname="DENS","QHYD","RH","QC","QR","QI","QS","QG","T","PRES","SFC_ALB_LW","SFC_ALB_SW","PRES","RADFLUX_SWDN","SFLX_SW_dn","TOAFLX_SW_dn",
#vcount=26
#vname="DENS","MOMZ","MOMX","MOMY","RHOT","QV","SFC_TEMP","DENS_t_advcv","MOMZ_t_advcv","MOMX_t_advcv","MOMY_t_advcv","RHOT_t_advcv","DENS_t_advch","MOMZ_t_advch","MOMX_t_advch","MOMY_t_advch","RHOT_t_advch","MOMZ_t_pg","MOMX_t_pg","MOMY_t_pg","MOMZ_t_ddiv","MOMX_t_ddiv","MOMY_t_ddiv","MOMZ_t_cf","MOMX_t_cf","MOMY_t_cf",
#vcount=24
#vname="TC_URB","QC_URB","SHFLX_urb","LHFLX_urb","GHFLX_urb","Rngrd_URB","TR_URB","TB_URB","TG_URB","SHR_URB","SHB_URB","SHG_URB","LHR_URB","LHB_URB","LHG_URB","RnR_URB","RnB_URB","RnG_URB","TRL_URB","TBL_URB","TGL_URB","GHR_URB","GHB_URB","GHG_URB",
#"PRES","T","U","V","MSLP","DENS","SFC_TEMP","LHFLX","height",

