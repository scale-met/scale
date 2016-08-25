#! /bin/bash -x
################################################################################
#
# ------ For Linux64 & GNU fortran&C & openmpi -----
#
################################################################################
export FORT_FMT_RECL=400


rm -f             ./input
ln -svf /data3/kenshi/scale_database/landuse/GLCCv2/Products ./input

# run
mpirun -np 4 /home/tyamaura/scale/scale-rm/test/case_real/check_mass/landuse/scale-rm_pp   pp.d01.conf   || exit 1
mpirun -np 4 /home/tyamaura/scale/scale-rm/test/case_real/check_mass/landuse/scale-rm_pp   pp.d02.conf   || exit 1
