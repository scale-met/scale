#!/bin/sh
#-------------------------------------------------------
#   Sample Script to link parent data with grads format
#-------------------------------------------------------
rm -f ./*_?????.grd
#
dt=21600
#
# set start date (UTC)
#
syear=2007
smonth=7
sday=14
shour=18
#
# set end date (UTC)
#
eyear=2007
emonth=07
eday=15
ehour=12
#
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
   #echo $1 ${endday[2]}
   return
}
#---------------------------------------------------
fn=0

syyyymmdd=`printf "%4.4d%2.2d%2.2d" $syear $smonth $sday`
eyyyymmdd=`printf "%4.4d%2.2d%2.2d" $eyear $emonth $eday`
echo $syyyymmdd $eyyyymmdd

iyear=${syear}
while [ ${iyear} -le ${eyear} ] ; do

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

        iday=1
        end_day=${endday[$imonth]}
        if [ ${iyear} -eq ${syear} ]&&[ ${imonth} -eq ${smonth} ] ;then
            iday=${sday}
        fi
        if [ ${iyear} -eq ${eyear} ]&&[ ${imonth} -eq ${emonth} ] ;then
            end_day=${eday}
        fi

        while [ ${iday} -le ${end_day} ] ;do

            iyyyymmdd=`printf "%4.4d%2.2d%2.2d" $iyear $imonth $iday`

            if [ ${iyyyymmdd} -eq ${syyyymmdd} ] ; then
                hour_list=`seq -s" " ${shour} 6 18`
            elif [ ${iyyyymmdd} -eq ${eyyyymmdd} ] ; then
                hour_list=`seq -s" " 0 6 ${ehour}`
            else
                hour_list=`seq -s" " 0 6 18`
            fi

            echo ${hour_list}
            for ihour in ${hour_list} ; do

               yyyy=`printf "%4.4d" $iyear`
               mm=`printf "%2.2d" $imonth`
               dd=`printf "%2.2d" $iday`
               hh=`printf "%2.2d" $ihour`

               fmtd_fn=`printf "%05d" $fn`
#-----
  for BND in FNLatm FNLland FNLsfc ; do

    dir="../../tools/FNL_output/${yyyy}${mm}"
    gpre=${BND}
    opre=${BND}

    #echo    ${dir}/${gpre}_${yyyy}${mm}${dd}${hh}.grd  ./${opre}_${fmtd_fn}.grd
    ln -svf ${dir}/${gpre}_${yyyy}${mm}${dd}${hh}.grd  ./${opre}_${fmtd_fn}.grd

  done # BND
#-----

               fn=`expr ${fn} \+ 1`

            done # hh

        iday=$((iday+1))
        done

    imonth=$((imonth+1))
    done

iyear=$((iyear+1))
done


