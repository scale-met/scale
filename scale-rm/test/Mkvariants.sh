#! /bin/bash -x

for dir in $(ls ./case)
do
   for subdir in $(ls ./case/${dir})
   do
      if [ ${subdir} = "500m" -o ${subdir} = "1000m" -o ${subdir} = "2000m" ]; then
         echo ./case/${dir}/${subdir}

         if [ ${dir} = "squallline" -o ${dir} = "supercell" -o ${dir} = "warmbubble" ]; then
            echo "--->kessler,tomita08"

            rm -rf ./case/${dir}/${subdir}_kessler
            cp -r  ./case/${dir}/${subdir} ./case/${dir}/${subdir}_kessler
            sed -e "s/SN14/KESSLER/g" ./case/${dir}/${subdir}/init.conf > ./case/${dir}/${subdir}_kessler/init.conf
            sed -e "s/SN14/KESSLER/g" ./case/${dir}/${subdir}/run.conf  > ./case/${dir}/${subdir}_kessler/run.conf

            rm -rf ./case/${dir}/${subdir}_tomita08
            cp -r  ./case/${dir}/${subdir} ./case/${dir}/${subdir}_tomita08
            sed -e "s/SN14/TOMITA08/g" ./case/${dir}/${subdir}/init.conf > ./case/${dir}/${subdir}_tomita08/init.conf
            sed -e "s/SN14/TOMITA08/g" ./case/${dir}/${subdir}/run.conf  > ./case/${dir}/${subdir}_tomita08/run.conf
         else
            rm -rf ./case/${dir}/${subdir}_kessler
            rm -rf ./case/${dir}/${subdir}_tomita08
         fi
      fi
   done
done
