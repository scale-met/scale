#!/bin/sh

#####
VAR=NOPREF
PREFETCH='-Knoprefetch'

sed -e "s/#PREFETCH#/${PREFETCH}/g" Makedef.tmpl > Makedef.${VAR}
export SCALE_SYS=${VAR}
make clean
make depends
make
mkdir -p ${VAR}
mv ./*.lst ./${VAR}/

#####
VAR=NOSWP
PREFETCH='-Knoswp'

sed -e "s/#PREFETCH#/${PREFETCH}/g" Makedef.tmpl > Makedef.${VAR}
export SCALE_SYS=${VAR}
make clean
make depends
make
mkdir -p ${VAR}
mv ./*.lst ./${VAR}/

#####
VAR=NOFISSION
PREFETCH='-Kloop_nofission'

sed -e "s/#PREFETCH#/${PREFETCH}/g" Makedef.tmpl > Makedef.${VAR}
export SCALE_SYS=${VAR}
make clean
make depends
make
mkdir -p ${VAR}
mv ./*.lst ./${VAR}/

#####
VAR=NOUNROLL
PREFETCH='-Knounroll'

sed -e "s/#PREFETCH#/${PREFETCH}/g" Makedef.tmpl > Makedef.${VAR}
export SCALE_SYS=${VAR}
make clean
make depends
make
mkdir -p ${VAR}
mv ./*.lst ./${VAR}/
