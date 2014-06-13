#! /bin/bash -x

if [ ! -f ../../../framework/cnvtopo/DEM50M/boundary.pe000000.nc ]; then
   cd ../../../framework/cnvtopo/DEM50M
   make     || exit 1
   make run || exit 1
   cd -
fi

for prc in `seq 0 3`
do
   PE=`printf %06d ${prc}`
   ln -svf ../../../framework/cnvtopo/DEM50M/boundary.pe${PE}.nc .
done
