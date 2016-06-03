#! /bin/bash -x
################################################################################
#
# ------ FOR Linux64 & gnu C&fortran & openmpi -----
#
################################################################################
export FORT_FMT_RECL=400
export GFORTRAN_UNBUFFERED_ALL=Y

ln -svf ../../../data/rad/PARAG.29 .
ln -svf ../../../data/rad/PARAPC.29 .
ln -svf ../../../data/rad/VARDATA.RM29 .
ln -svf ../../../data/rad/cira.nc .
ln -svf ../../../data/rad/MIPAS/day.atm .
ln -svf ../../../data/rad/MIPAS/equ.atm .
ln -svf ../../../data/rad/MIPAS/sum.atm .
ln -svf ../../../data/rad/MIPAS/win.atm .
ln -svf ../../../data/land/param.bucket.conf .
ln -svf ../topo/topo.d01.pe000000.nc .
ln -svf ../landuse/landuse.d01.pe000000.nc .
ln -svf ../init/boundary.d01.pe000000.nc .
ln -svf ../init/init.d01_00011685600.000.pe000000.nc .
ln -svf ../topo/topo.d02.pe000000.nc .
ln -svf ../landuse/landuse.d02.pe000000.nc .
ln -svf ../init/init.d02_00011685600.000.pe000000.nc .
ln -svf ../topo/topo.d01.pe000001.nc .
ln -svf ../landuse/landuse.d01.pe000001.nc .
ln -svf ../init/boundary.d01.pe000001.nc .
ln -svf ../init/init.d01_00011685600.000.pe000001.nc .
ln -svf ../topo/topo.d02.pe000001.nc .
ln -svf ../landuse/landuse.d02.pe000001.nc .
ln -svf ../init/init.d02_00011685600.000.pe000001.nc .
ln -svf ../topo/topo.d01.pe000002.nc .
ln -svf ../landuse/landuse.d01.pe000002.nc .
ln -svf ../init/boundary.d01.pe000002.nc .
ln -svf ../init/init.d01_00011685600.000.pe000002.nc .
ln -svf ../topo/topo.d02.pe000002.nc .
ln -svf ../landuse/landuse.d02.pe000002.nc .
ln -svf ../init/init.d02_00011685600.000.pe000002.nc .
ln -svf ../topo/topo.d01.pe000003.nc .
ln -svf ../landuse/landuse.d01.pe000003.nc .
ln -svf ../init/boundary.d01.pe000003.nc .
ln -svf ../init/init.d01_00011685600.000.pe000003.nc .
ln -svf ../topo/topo.d02.pe000003.nc .
ln -svf ../landuse/landuse.d02.pe000003.nc .
ln -svf ../init/init.d02_00011685600.000.pe000003.nc .

# run

mpirun -np 8 /home/tyamaura/scale/scale-rm/test/case_real/check_mass/run/scale-rm launch.conf || exit

################################################################################
