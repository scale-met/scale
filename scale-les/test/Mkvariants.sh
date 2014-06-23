#! /bin/bash -x

for dir in $(ls ./case)
do
   for subdir in $(ls ./case/${dir})
   do
      if [ ${subdir} = "500m" -o ${subdir} = "1000m" -o ${subdir} = "2000m" ]; then
         echo ./case/${dir}/${subdir}

         rm -rf ./case/${dir}/${subdir}_hevi
         cp -r  ./case/${dir}/${subdir} ./case/${dir}/${subdir}_hevi
         sed -e "s/HEVE/HEVI/g" ./case/${dir}/${subdir}/run.conf > ./case/${dir}/${subdir}_hevi/run.conf

         rm -rf ./case/${dir}/${subdir}_hivi
         cp -r  ./case/${dir}/${subdir} ./case/${dir}/${subdir}_hivi
         sed -e "s/HEVE/HIVI/g" ./case/${dir}/${subdir}/run.conf > ./case/${dir}/${subdir}_hivi/run.conf

      fi
   done
done
