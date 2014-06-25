#! /bin/bash -x

for dir in $(ls ./case)
do
   for subdir in $(ls ./case/${dir})
   do
      if [ ${subdir} = "500m" -o ${subdir} = "1000m" -o ${subdir} = "2000m" ]; then
         echo ./case/${dir}/${subdir}

         if [ ${dir} = "radiation" -o ${dir} = "rad-conv" -o ${dir} = "urban" ]; then
            echo "--->skip"
            rm -rf ./case/${dir}/${subdir}_hevi
            rm -rf ./case/${dir}/${subdir}_hivi
         else
            echo "--->HEVI,HIVI"

            rm -rf ./case/${dir}/${subdir}_hevi
            cp -r  ./case/${dir}/${subdir} ./case/${dir}/${subdir}_hevi
            sed -e "s/HEVE/HEVI/g" ./case/${dir}/${subdir}/run.conf > ./case/${dir}/${subdir}_hevi/run.conf

            rm -rf ./case/${dir}/${subdir}_hivi
            cp -r  ./case/${dir}/${subdir} ./case/${dir}/${subdir}_hivi
            sed -e "s/HEVE/HIVI/g" ./case/${dir}/${subdir}/run.conf > ./case/${dir}/${subdir}_hivi/run.conf
         fi

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
