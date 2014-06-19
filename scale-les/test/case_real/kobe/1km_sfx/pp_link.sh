#! /bin/bash -x

if [ ! -f ../../../framework/cnvtopo/DEM50m/boundary_topo.pe000000.nc ]; then
   cd ../../../framework/cnvtopo/DEM50m
   make     || exit 1
   make run || exit 1
   cd -
fi

for prc in `seq 0 3`
do
   PE=`printf %06d ${prc}`
   ln -svf ../../../framework/cnvtopo/DEM50m/boundary_topo.pe${PE}.nc .
done

if [ ! -f ../../../framework/cnvlanduse/LU100M/boundary_landuse.pe000000.nc ]; then
   cd ../../../framework/cnvlanduse/LU100M
   make     || exit 1
   make run || exit 1
   cd -
fi

for prc in `seq 0 3`
do
   PE=`printf %06d ${prc}`
   ln -svf ../../../framework/cnvlanduse/LU100M/boundary_landuse.pe${PE}.nc .
done

