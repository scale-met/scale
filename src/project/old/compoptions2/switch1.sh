#!/bin/sh

#####
VAR=NOSWP
COMPOPTION='-Knoswp'

sed -e "s/#COMPOPTION#/${COMPOPTION}/g" Makedef.tmpl > Makedef.${VAR}
export SCALE_SYS=${VAR}
make clean
make depends
make
mkdir -p ${VAR}
mv ./*.lst ./${VAR}/

#####
VAR=STRIPING2
COMPOPTION='-Knoswp,striping=2'

sed -e "s/#COMPOPTION#/${COMPOPTION}/g" Makedef.tmpl > Makedef.${VAR}
export SCALE_SYS=${VAR}
make clean
make depends
make
mkdir -p ${VAR}
mv ./*.lst ./${VAR}/

#####
VAR=STRIPING3
COMPOPTION='-Knoswp,striping=3'

sed -e "s/#COMPOPTION#/${COMPOPTION}/g" Makedef.tmpl > Makedef.${VAR}
export SCALE_SYS=${VAR}
make clean
make depends
make
mkdir -p ${VAR}
mv ./*.lst ./${VAR}/

#####
VAR=STRIPING4
COMPOPTION='-Knoswp,striping=4'

sed -e "s/#COMPOPTION#/${COMPOPTION}/g" Makedef.tmpl > Makedef.${VAR}
export SCALE_SYS=${VAR}
make clean
make depends
make
mkdir -p ${VAR}
mv ./*.lst ./${VAR}/
