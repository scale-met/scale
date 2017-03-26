#!/bin/bash
#
#  This is converter from NCEP-FNL grib1/grib2 to GrADS format for SCALE.
#
#  The way to execute the shell is as follows.
#     $ sh convert_FNL-grib2grads.sh ${START_DATE} ${END_DATE} ${RDIR} ${WDIR}
#
#  The program reads original data in ${RDIR} and then output to ${WDIR}.
#
#  GRIB1 structure: RDIR/grib1/${YEAR}/fnl_${YYYY}${MM}${DD}_${HH}_00.grib1
#  GRIB2 structure: RDIR/grib2/${YEAR}/fnl_${YYYY}${MM}${DD}_${HH}_00.grib2
#
#-------------------------------------------------------------------------------

# set date you want to convert
START_DATE=${1:-1999080100}
END_DATE=${2:-1999080100}

# path to directory with grib2 data
RDIR="${3:-.}"
# output directory
WDIR="${4:-.}"

### --- make data --- ###
START_YEAR=`echo ${START_DATE} | cut -c 1-4 | bc`
START_MON=`echo ${START_DATE} | cut -c 5-6 | bc`
START_DAY=`echo ${START_DATE} | cut -c 7-8 | bc`
START_HOUR=`echo ${START_DATE} | cut -c 9-10 | bc`
END_YEAR=`echo ${END_DATE} | cut -c 1-4 | bc`
END_MON=`echo ${END_DATE} | cut -c 5-6 | bc`
END_DAY=`echo ${END_DATE} | cut -c 7-8 | bc`
END_HOUR=`echo ${END_DATE} | cut -c 9-10 | bc`

# change file format from GRIB1 to GRIB2
CHECK_DATE1="2007-12-06 06:00:00"
# change vertical order
CHECK_DATE2="2009-12-15 06:00:00"
# change the name of soil temperature
CHECK_DATE3="2015-01-14 00:00:00"

UNIX_CHECK1=`date -u -d "$CHECK_DATE1" +%s`
UNIX_CHECK2=`date -u -d "$CHECK_DATE2" +%s`
UNIX_CHECK3=`date -u -d "$CHECK_DATE3" +%s`

YEAR=${START_YEAR}
while [ ${YEAR} -le ${END_YEAR} ]
do
  # leap year process
  if [ `expr ${YEAR} % 4`   -eq 0 -a \
       `expr ${YEAR} % 100` -ne 0 -o \
       `expr ${YEAR} % 400` -eq 0    ]; then
    MAXDAYS=( 31 29 31 30 31 30 31 31 30 31 30 31 )
  else
    MAXDAYS=( 31 28 31 30 31 30 31 31 30 31 30 31 )
  fi

  YYYY=`printf "%04d" ${YEAR}`

  if [ ${YEAR} -eq ${START_YEAR} ]; then
      MON=${START_MON}
  else
      MON=1
  fi
  if [ ${YEAR} -eq ${END_YEAR} ]; then
      MAXMON=${END_MON}
  else
      MAXMON=12
  fi

  while [ ${MON} -le ${MAXMON} ]
  do
    MM=`printf "%02d" ${MON}`

    # make output directory
    if [ ! -d ${WDIR}/${YYYY}${MM} ]; then
      mkdir -p ${WDIR}/${YYYY}${MM}
    fi

    if [ ${YEAR}${MON} -eq ${START_YEAR}${START_MON} ]; then
        DAY=${START_DAY}
    else
        DAY=1
    fi
    if [ ${YEAR}${MON} -eq ${END_YEAR}${END_MON} ]; then
        MAXDAY=${END_DAY}
    else
        MAXDAY=${MAXDAYS[`expr ${MON} - 1`]}
    fi

    while [ $DAY -le ${MAXDAY} ]
    do
      DD=`printf "%02d" $DAY`

      if [ ${YEAR}${MON}${DAY} -eq ${START_YEAR}${START_MON}${START_DAY} ]; then
          HOUR=${START_HOUR}
      else
          HOUR=0
      fi
      if [ ${YEAR}${MON}${DAY} -eq ${END_YEAR}${END_MON}${END_DAY} ]; then
          MAXHOUR=${END_HOUR}
      else
          MAXHOUR=18
      fi

      while [ $HOUR -le ${MAXHOUR} ]
      do
        HH=`printf "%02d" ${HOUR}`

        DATE="${YYYY}-${MM}-${DD} ${HH}:00:00"
        UNIX_NOW=`date -u -d "$DATE" +%s`

        if [ $UNIX_NOW -le $UNIX_CHECK1 ]; then
          FTYPE="grib1"
        else
          FTYPE="grib2"
        fi
        if [ $UNIX_NOW -le $UNIX_CHECK2 ]; then
          ZREV=""
        else
          ZREV="-r"
        fi
        if [ $UNIX_NOW -le $UNIX_CHECK3 ]; then
          TSOIL="TMP"
        else
          TSOIL="TSOIL"
        fi

        RFILE=${RDIR}/${FTYPE}/${YYYY}/fnl_${YYYY}${MM}${DD}_${HH}_00.${FTYPE}

        WFILE_ATM=${WDIR}/${YYYY}${MM}/FNL_ATM_${YYYY}${MM}${DD}${HH}.grd
        WFILE_SFC=${WDIR}/${YYYY}${MM}/FNL_SFC_${YYYY}${MM}${DD}${HH}.grd
        WFILE_LND=${WDIR}/${YYYY}${MM}/FNL_LND_${YYYY}${MM}${DD}${HH}.grd

        echo "reading ... $RFILE"

        if [ $UNIX_NOW -le $UNIX_CHECK1 ]; then
          # 3D
          wgrib ${RFILE} | grep ":HGT:"  | grep "mb:anl:" | wgrib -i ${RFILE} -nh -ieee -o HGTprs.${YYYY}${MM}${DD}${HH}.grd  > /dev/null 2>&1
          wgrib ${RFILE} | grep ":UGRD:" | grep "mb:anl:" | wgrib -i ${RFILE} -nh -ieee -o UGRDprs.${YYYY}${MM}${DD}${HH}.grd > /dev/null 2>&1
          wgrib ${RFILE} | grep ":VGRD:" | grep "mb:anl:" | wgrib -i ${RFILE} -nh -ieee -o VGRDprs.${YYYY}${MM}${DD}${HH}.grd > /dev/null 2>&1
          wgrib ${RFILE} | grep ":TMP:"  | grep "mb:anl:" | wgrib -i ${RFILE} -nh -ieee -o TMPprs.${YYYY}${MM}${DD}${HH}.grd  > /dev/null 2>&1
          wgrib ${RFILE} | grep ":RH:"   | grep "mb:anl:" | wgrib -i ${RFILE} -nh -ieee -o RHprs.${YYYY}${MM}${DD}${HH}.grd   > /dev/null 2>&1
          # 2D
          wgrib ${RFILE} | grep ":PRMSL:" | grep ":MSL:"            | wgrib -i ${RFILE} -nh -ieee -o PRMSLmsl.${YYYY}${MM}${DD}${HH}.grd >/dev/null 2>&1
          wgrib ${RFILE} | grep ":PRES:"  | grep ":sfc:"            | wgrib -i ${RFILE} -nh -ieee -o PRESsfc.${YYYY}${MM}${DD}${HH}.grd  >/dev/null 2>&1
          wgrib ${RFILE} | grep ":TMP:"   | grep ":sfc:"            | wgrib -i ${RFILE} -nh -ieee -o TMPsfc.${YYYY}${MM}${DD}${HH}.grd   >/dev/null 2>&1
          wgrib ${RFILE} | grep ":HGT:"   | grep ":sfc:"            | wgrib -i ${RFILE} -nh -ieee -o HGTsfc.${YYYY}${MM}${DD}${HH}.grd   >/dev/null 2>&1
          wgrib ${RFILE} | grep ":LAND:"  | grep ":sfc:"            | wgrib -i ${RFILE} -nh -ieee -o LANDsfc.${YYYY}${MM}${DD}${HH}.grd  >/dev/null 2>&1
          wgrib ${RFILE} | grep ":UGRD:"  | grep ":10 m above gnd:" | wgrib -i ${RFILE} -nh -ieee -o UGRD10m.${YYYY}${MM}${DD}${HH}.grd  >/dev/null 2>&1
          wgrib ${RFILE} | grep ":VGRD:"  | grep ":10 m above gnd:" | wgrib -i ${RFILE} -nh -ieee -o VGRD10m.${YYYY}${MM}${DD}${HH}.grd  >/dev/null 2>&1
          wgrib ${RFILE} | grep ":TMP:"   | grep ":2 m above gnd:"  | wgrib -i ${RFILE} -nh -ieee -o TMP2m.${YYYY}${MM}${DD}${HH}.grd    >/dev/null 2>&1
          wgrib ${RFILE} | grep ":RH:"    | grep ":2 m above gnd:"  | wgrib -i ${RFILE} -nh -ieee -o RH2m.${YYYY}${MM}${DD}${HH}.grd     >/dev/null 2>&1
          # Land: 2-layer soil until 2005/05/31-06UTC
          wgrib ${RFILE} | grep ":TMP:"   | grep "down:" | wgrib -i ${RFILE} -nh -ieee -o TSOIL.${YYYY}${MM}${DD}${HH}.grd >/dev/null 2>&1
          wgrib ${RFILE} | grep ":SOILW:" | grep "down:" | wgrib -i ${RFILE} -nh -ieee -o SOILW.${YYYY}${MM}${DD}${HH}.grd >/dev/null 2>&1
        else
          # 3D: 26-layer atmosphere until 2016/05/11-06UTC
          wgrib2 ${RFILE} -match ":HGT:"  -match "mb:anl:" | sort -n $ZREV | wgrib2 -i ${RFILE} -order "we:ns" -no_header -ieee HGTprs.${YYYY}${MM}${DD}${HH}.grd  > /dev/null 2>&1
          wgrib2 ${RFILE} -match ":UGRD:" -match "mb:anl:" | sort -n $ZREV | wgrib2 -i ${RFILE} -order "we:ns" -no_header -ieee UGRDprs.${YYYY}${MM}${DD}${HH}.grd > /dev/null 2>&1
          wgrib2 ${RFILE} -match ":VGRD:" -match "mb:anl:" | sort -n $ZREV | wgrib2 -i ${RFILE} -order "we:ns" -no_header -ieee VGRDprs.${YYYY}${MM}${DD}${HH}.grd > /dev/null 2>&1
          wgrib2 ${RFILE} -match ":TMP:"  -match "mb:anl:" | sort -n $ZREV | wgrib2 -i ${RFILE} -order "we:ns" -no_header -ieee TMPprs.${YYYY}${MM}${DD}${HH}.grd  > /dev/null 2>&1
          wgrib2 ${RFILE} -match ":RH:"   -match "mb:anl:" | sort -n $ZREV | wgrib2 -i ${RFILE} -order "we:ns" -no_header -ieee RHprs.${YYYY}${MM}${DD}${HH}.grd   > /dev/null 2>&1
          # 2D
          wgrib2 ${RFILE} -match ":PRMSL:" -match ":mean sea level:"    -order "we:ns" -no_header -ieee PRMSLmsl.${YYYY}${MM}${DD}${HH}.grd >/dev/null 2>&1
          wgrib2 ${RFILE} -match ":PRES:"  -match ":surface:"           -order "we:ns" -no_header -ieee PRESsfc.${YYYY}${MM}${DD}${HH}.grd  >/dev/null 2>&1
          wgrib2 ${RFILE} -match ":TMP:"   -match ":surface:"           -order "we:ns" -no_header -ieee TMPsfc.${YYYY}${MM}${DD}${HH}.grd   >/dev/null 2>&1
          wgrib2 ${RFILE} -match ":HGT:"   -match ":surface:"           -order "we:ns" -no_header -ieee HGTsfc.${YYYY}${MM}${DD}${HH}.grd   >/dev/null 2>&1
          wgrib2 ${RFILE} -match ":LAND:"  -match ":surface:"           -order "we:ns" -no_header -ieee LANDsfc.${YYYY}${MM}${DD}${HH}.grd  >/dev/null 2>&1
          wgrib2 ${RFILE} -match ":UGRD:"  -match ":10 m above ground:" -order "we:ns" -no_header -ieee UGRD10m.${YYYY}${MM}${DD}${HH}.grd  >/dev/null 2>&1
          wgrib2 ${RFILE} -match ":VGRD:"  -match ":10 m above ground:" -order "we:ns" -no_header -ieee VGRD10m.${YYYY}${MM}${DD}${HH}.grd  >/dev/null 2>&1
          wgrib2 ${RFILE} -match ":TMP:"   -match ":2 m above ground:"  -order "we:ns" -no_header -ieee TMP2m.${YYYY}${MM}${DD}${HH}.grd    >/dev/null 2>&1
          wgrib2 ${RFILE} -match ":RH:"    -match ":2 m above ground:"  -order "we:ns" -no_header -ieee RH2m.${YYYY}${MM}${DD}${HH}.grd     >/dev/null 2>&1
          # Land
          wgrib2 ${RFILE} -match ":$TSOIL:" -match "below ground:" -order "we:ns" -no_header -ieee TSOIL.${YYYY}${MM}${DD}${HH}.grd >/dev/null 2>&1
          wgrib2 ${RFILE} -match ":SOILW:"  -match "below ground:" -order "we:ns" -no_header -ieee SOILW.${YYYY}${MM}${DD}${HH}.grd >/dev/null 2>&1
        fi

        rm -f ${WFILE_ATM}
        rm -f ${WFILE_SFC}
        rm -f ${WFILE_LND}

        cat HGTprs.${YYYY}${MM}${DD}${HH}.grd  \
            UGRDprs.${YYYY}${MM}${DD}${HH}.grd \
            VGRDprs.${YYYY}${MM}${DD}${HH}.grd \
            TMPprs.${YYYY}${MM}${DD}${HH}.grd  \
            RHprs.${YYYY}${MM}${DD}${HH}.grd   \
        >> ${WFILE_ATM}
        cat PRMSLmsl.${YYYY}${MM}${DD}${HH}.grd \
            PRESsfc.${YYYY}${MM}${DD}${HH}.grd  \
            TMPsfc.${YYYY}${MM}${DD}${HH}.grd   \
            HGTsfc.${YYYY}${MM}${DD}${HH}.grd   \
            LANDsfc.${YYYY}${MM}${DD}${HH}.grd  \
            UGRD10m.${YYYY}${MM}${DD}${HH}.grd  \
            VGRD10m.${YYYY}${MM}${DD}${HH}.grd  \
            TMP2m.${YYYY}${MM}${DD}${HH}.grd    \
            RH2m.${YYYY}${MM}${DD}${HH}.grd     \
        >> ${WFILE_SFC}
        cat TSOIL.${YYYY}${MM}${DD}${HH}.grd \
            SOILW.${YYYY}${MM}${DD}${HH}.grd \
        >> ${WFILE_LND}

        rm -f HGTprs.${YYYY}${MM}${DD}${HH}.grd   \
              UGRDprs.${YYYY}${MM}${DD}${HH}.grd  \
              VGRDprs.${YYYY}${MM}${DD}${HH}.grd  \
              TMPprs.${YYYY}${MM}${DD}${HH}.grd   \
              RHprs.${YYYY}${MM}${DD}${HH}.grd    \
              PRMSLmsl.${YYYY}${MM}${DD}${HH}.grd \
              PRESsfc.${YYYY}${MM}${DD}${HH}.grd  \
              HGTsfc.${YYYY}${MM}${DD}${HH}.grd   \
              TMPsfc.${YYYY}${MM}${DD}${HH}.grd   \
              LANDsfc.${YYYY}${MM}${DD}${HH}.grd  \
              UGRD10m.${YYYY}${MM}${DD}${HH}.grd  \
              VGRD10m.${YYYY}${MM}${DD}${HH}.grd  \
              TMP2m.${YYYY}${MM}${DD}${HH}.grd    \
              RH2m.${YYYY}${MM}${DD}${HH}.grd     \
              TSOIL.${YYYY}${MM}${DD}${HH}.grd    \
              SOILW.${YYYY}${MM}${DD}${HH}.grd

        HOUR=`expr $HOUR + 6`
      done 
      DAY=`expr $DAY + 1`
    done
    MON=`expr $MON + 1`
  done
  YEAR=`expr $YEAR + 1`
done

echo "conversion finished ..."
