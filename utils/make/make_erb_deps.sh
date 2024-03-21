#!/bin/bash -x

cd ../../scalelib/src

echo
echo "Update .erb -> .F90 (scalelib)"

erbs=(`find . -name "*.erb"`)

for ferb in ${erbs[@]}
do
   ferb_base=`basename ${ferb}`
   if [ "${ferb_base}" == "scale_atmos_dyn_fvm_flux_udcd.F90.erb" ];then
      udcds=(`find . -name "scale_atmos_dyn_fvm_flux_*.F90"`)
      for fudcd in ${udcds[@]}
      do
         echo
         echo "source : ${ferb}"
         echo "target : ${fudcd}"
         fname=${fudcd} erb -T - ${ferb#./} > ${fudcd}
      done
   else
      echo
      echo "source : ${ferb}"
      echo "target : ${ferb%.erb}"
      erb ${ferb#./} > ${ferb%.erb}
   fi
done

echo
echo "Update dependency (scalelib)"
./makedepend .

cd -

cd ../../scale-rm/src

echo
echo "Update dependency (scale-rm)"
./makedepend .

cd -

cd ../../scale-gm/src

echo
echo "Update dependency (scale-gm)"
./makedepend .

cd -
