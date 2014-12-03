#!/bin/bash

FFLAGS="-O3 -assume byterecl -convert big_endian"
LIBS="-L${NETCDF}/lib -L${HDF5}/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lm -lz"
INCLUDES="-I${NETCDF}/include"

rm -f grads2netcdf *.mod
ifort ${FFLAGS}  src/mod_netcdf.f90 src/nicam_grads2netcdf.f90 -o grads2netcdf ${LIBS} ${INCLUDES}
rm -f *.mod

##
## topography
##

# grads format
gfile=topog/02560x01280.zorg.torg/topog.grd
# netcdf format
nfile=topog/02560x01280.zorg.torg/topog.nc
# reference
# copied from or path to ../19990501/02560x01280.zorg.torg/ss_tem_sfc.nc
rfile=reffile.nc

cat > conv.in <<EOF
&fileinfo
gfile='${gfile}'
nfile='${nfile}'
rfile='${rfile}'
/
&gridinfo
nx=2560
ny=1280
varn='topog'
vard='Topography'
varu='m'
/
EOF

./grads2netcdf


##
## landuse
##

# grads format
gfile=lsmask/02560x01280.zorg.torg/lsmask.grd
# netcdf format
nfile=lsmask/02560x01280.zorg.torg/lsmask.nc
# reference
# copied from or path to ../19990501/02560x01280.zorg.torg/ss_tem_sfc.nc
rfile=reffile.nc

cat > conv.in <<EOF
&fileinfo
gfile='${gfile}'
nfile='${nfile}'
rfile='${rfile}'
/
&gridinfo
nx=2560
ny=1280
varn='lsmask'
vard='LANDMASK'
varu='fraction'
/
EOF

./grads2netcdf

