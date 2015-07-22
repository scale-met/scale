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
syear=2014
smonth=8
sday=10
eyear=2014
emonth=8
eday=10

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

                gfile=${dir_path_to_grib}/fnl_${yyyy}${mm}${dd}_${hh}_00.grib2
                ofile_land=${output_dir}/${yyyy}${mm}/FNLland_${yyyy}${mm}${dd}${hh}.grd
                ofile_sfc=${output_dir}/${yyyy}${mm}/FNLsfc_${yyyy}${mm}${dd}${hh}.grd
                ofile_atm=${output_dir}/${yyyy}${mm}/FNLatm_${yyyy}${mm}${dd}${hh}.grd

                tfile="tmp.grb2"
                rm -f ${tfile} ${ofile_land}
                #--Land data
                wgrib2 ${gfile} -match ":LAND:"  -match ":surface:"             -grib ${tfile}
                wgrib2 ${gfile} -match ":TMP:"   -match ":surface:"     -append -grib ${tfile}
                wgrib2 ${gfile} -match ":TMP:"   -match "below ground:" -append -grib ${tfile}
                wgrib2 ${gfile} -match ":SOILW:" -match "below ground:" -append -grib ${tfile}
                wgrib2 ${tfile} -no_header -ieee ${ofile_land} -g2clib 0

                rm -f ${tfile} ${ofile_sfc}
                #--Surface data
                wgrib2 ${gfile} -match ":PRES:"  -match ":mean sea level:"            -grib ${tfile}
                wgrib2 ${gfile} -match ":PRES:"  -match ":surface:"           -append -grib ${tfile}
                wgrib2 ${gfile} -match ":UGRD:"  -match ":10 m above ground:" -append -grib ${tfile}
                wgrib2 ${gfile} -match ":VGRD:"  -match ":10 m above ground:" -append -grib ${tfile}
                wgrib2 ${gfile} -match ":TMP:"   -match ":2 m above ground:"  -append -grib ${tfile}
                wgrib2 ${gfile} -match ":RH:"    -match ":2 m above ground:"  -append -grib ${tfile}
                wgrib2 ${tfile} -no_header -ieee ${ofile_sfc} -g2clib 0

                rm -f ${tfile} ${ofile_atm}
                #--Upper data
                wgrib2 ${gfile} -match ":HGT:"  -match "mb" -not "mb above ground" | sort -n -r | wgrib2 -i ${gfile} -grib ${tfile}
                wgrib2 ${gfile} -match ":UGRD:" -match "mb" -not "mb above ground" | sort -n -r | wgrib2 -i ${gfile} -append -grib ${tfile}
                wgrib2 ${gfile} -match ":VGRD:" -match "mb" -not "mb above ground" | sort -n -r | wgrib2 -i ${gfile} -append -grib ${tfile}
                wgrib2 ${gfile} -match ":TMP:"  -match "mb" -not "mb above ground" | sort -n -r | wgrib2 -i ${gfile} -append -grib ${tfile}
                wgrib2 ${gfile} -match ":RH:"   -match "mb" -not "mb above ground" | sort -n -r | wgrib2 -i ${gfile} -append -grib ${tfile}
                wgrib2 ${tfile} -no_header -ieee ${ofile_atm} -g2clib 0

            done

        iday=$((iday+1))
        done

    imonth=$((imonth+1))
    done

iyear=$((iyear+1))
done

