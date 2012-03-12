#!/bin/sh
#
######
#VAR=UNROLL
#COMPOPTION='-Kunroll'
#
#sed -e "s/#COMPOPTION#/${COMPOPTION}/g" Makedef.tmpl > Makedef.${VAR}
#export SCALE_SYS=${VAR}
#make clean
#make depends
#make
#mkdir -p ${VAR}
#mv ./*.lst ./${VAR}/

#####
VAR=UNROLL2
COMPOPTION='-Kunroll=2'

sed -e "s/#COMPOPTION#/${COMPOPTION}/g" Makedef.tmpl > Makedef.${VAR}
export SCALE_SYS=${VAR}
make clean
make depends
make
mkdir -p ${VAR}
mv ./*.lst ./${VAR}/
