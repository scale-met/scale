#!/bin/bash
#
#  This is converter from NCEP-FNL with grib2 format to grads format for SCALE.
#
#  The way to execute the shell is as follows.
#     $ sh convert_grib2grads.sh
#
#  The program reads original data in ${dir_path_to_grib} and then output grads format data to ${output_dir}.
#
#-------------------------------------------------------------------------------
# set date you want to convert
syear=2007
smonth=7
sday=14
eyear=2007
emonth=7
eday=15

# path to directory with GSM grib2 data
dir_path_to_grib="./"

output_dir="./FNL_output"
#---------------------------------------------------
function uruu() {
   endday=(0 31 28 31 30 31 30 31 31 30 31 30 31)
   if [ $(($1%4)) -eq 0 ] ;then
       endday[2]=29
       if [ $(($1%100)) -eq 0 ] ;then
           endday[2]=28
           if [ $(($1%400)) -eq 0 ] ;then
           endday[2]=29
           fi
       fi
   fi
   echo $1 ${endday[2]}
   return
}
#---------------------------------------------------
iyear=${syear}
while [ ${iyear} -le ${eyear} ] ;do

    uruu $iyear

    imonth=1
    end_month=12
    if [ ${iyear} -eq ${syear} ] ;then
        imonth=${smonth}
    fi
    if [ ${iyear} -eq ${eyear} ] ;then
        end_month=${emonth}
    fi

    while [ ${imonth} -le ${end_month} ] ;do

        yyyy=${iyear}
        mm=`printf "%2.2d" $imonth`
        if [ ! -d ${output_dir}/${yyyy}${mm} ] ;then
            mkdir -p ${output_dir}/${yyyy}${mm}
        fi

        iday=1
        end_day=${endday[$imonth]}
        if [ ${iyear} -eq ${syear} ]&&[ ${imonth} -eq ${smonth} ] ;then
            iday=${sday}
        fi
        if [ ${iyear} -eq ${eyear} ]&&[ ${imonth} -eq ${emonth} ] ;then
            end_day=${eday}
        fi

        while [ ${iday} -le ${end_day} ] ;do

            dd=`printf "%2.2d" $iday`

            for hh in 00 06 12 18 ;do

                gfile=${dir_path_to_grib}/fnl_${yyyy}${mm}${dd}_${hh}_00.grib1
                ofile_land=${output_dir}/${yyyy}${mm}/FNLland_${yyyy}${mm}${dd}${hh}.grd
                ofile_sfc=${output_dir}/${yyyy}${mm}/FNLsfc_${yyyy}${mm}${dd}${hh}.grd
                ofile_atm=${output_dir}/${yyyy}${mm}/FNLatm_${yyyy}${mm}${dd}${hh}.grd

                if [ ! -f $gfile ]; then
                    continue
                fi

                tfile="tmp.grb"
                rm -f ${tfile} ${ofile_land}
                #--Land data
                wgrib ${gfile} | grep ":LAND:" | grep ":sfc:" | wgrib -i ${gfile}         -grib -o ${tfile}
                wgrib ${gfile} | grep ":TMP:"  | grep ":sfc:" | wgrib -i ${gfile} -append -grib -o ${tfile}
                wgrib ${gfile} | grep ":TMP:"  | grep "down:" | wgrib -i ${gfile} -append -grib -o ${tfile}
                wgrib ${gfile} | grep ":SOILW:"| grep "down:" | wgrib -i ${gfile} -append -grib -o ${tfile}
                wgrib ${tfile} | wgrib -i ${tfile} -nh -ieee -o ${ofile_land}

                rm -f ${tfile} ${ofile_sfc}
                #--Surface data
                wgrib ${gfile} | grep ":PRMSL:" | grep ":MSL:"            | wgrib -i ${gfile}         -grib -o ${tfile}
                wgrib ${gfile} | grep ":PRES:"  | grep ":sfc:"            | wgrib -i ${gfile} -append -grib -o ${tfile}
                wgrib ${gfile} | grep ":UGRD:"  | grep ":10 m above gnd:" | wgrib -i ${gfile} -append -grib -o ${tfile}
                wgrib ${gfile} | grep ":VGRD:"  | grep ":10 m above gnd:" | wgrib -i ${gfile} -append -grib -o ${tfile}
                wgrib ${gfile} | grep ":TMP:"   | grep ":2 m above gnd:"  | wgrib -i ${gfile} -append -grib -o ${tfile}
                wgrib ${gfile} | grep ":RH:"    | grep ":2 m above gnd:"  | wgrib -i ${gfile} -append -grib -o ${tfile}
                wgrib ${tfile} | wgrib -i ${tfile} -nh -ieee -o ${ofile_sfc}

                rm -f ${tfile} ${ofile_atm}
                #--Upper data
                wgrib ${gfile} | grep ":HGT:"  | grep "mb" | grep -v "mb above gnd" | sort -n  | wgrib -i ${gfile}         -grib -o ${tfile}
                wgrib ${gfile} | grep ":UGRD:" | grep "mb" | grep -v "mb above gnd" | sort -n  | wgrib -i ${gfile} -append -grib -o ${tfile}
                wgrib ${gfile} | grep ":VGRD:" | grep "mb" | grep -v "mb above gnd" | sort -n  | wgrib -i ${gfile} -append -grib -o ${tfile}
                wgrib ${gfile} | grep ":TMP:"  | grep "mb" | grep -v "mb above gnd" | sort -n  | wgrib -i ${gfile} -append -grib -o ${tfile}
                wgrib ${gfile} | grep ":RH:"   | grep "mb" | grep -v "mb above gnd" | sort -n  | wgrib -i ${gfile} -append -grib -o ${tfile}
                wgrib ${tfile} | wgrib -i ${tfile} -nh -ieee -o ${ofile_atm}

            done

        iday=$((iday+1))
        done

    imonth=$((imonth+1))
    done

iyear=$((iyear+1))
done


