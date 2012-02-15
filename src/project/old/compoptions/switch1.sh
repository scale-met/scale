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
VAR=PREFINF
PREFETCH='-Kprefetch_infer,prefetch_strong_L2,prefetch_iteration_L2=50'

sed -e "s/#PREFETCH#/${PREFETCH}/g" Makedef.tmpl > Makedef.${VAR}
export SCALE_SYS=${VAR}
make clean
make depends
make
mkdir -p ${VAR}
mv ./*.lst ./${VAR}/

#####
VAR=PREFCOND
PREFETCH='-Kprefetch_conditional,prefetch_strong_L2,prefetch_iteration_L2=50'

sed -e "s/#PREFETCH#/${PREFETCH}/g" Makedef.tmpl > Makedef.${VAR}
export SCALE_SYS=${VAR}
make clean
make depends
make
mkdir -p ${VAR}
mv ./*.lst ./${VAR}/

#####
VAR=PREFSTR
PREFETCH='-Kprefetch_stride,prefetch_strong_L2,prefetch_iteration_L2=50'

sed -e "s/#PREFETCH#/${PREFETCH}/g" Makedef.tmpl > Makedef.${VAR}
export SCALE_SYS=${VAR}
make clean
make depends
make
mkdir -p ${VAR}
mv ./*.lst ./${VAR}/

#####
VAR=L1W0L2W50
PREFETCH='-Kprefetch_nostrong_L2,prefetch_iteration_L2=50'

sed -e "s/#PREFETCH#/${PREFETCH}/g" Makedef.tmpl > Makedef.${VAR}
export SCALE_SYS=${VAR}
make clean
make depends
make
mkdir -p ${VAR}
mv ./*.lst ./${VAR}/

#####
VAR=L1W0L2S50
PREFETCH='-Kprefetch_strong_L2,prefetch_iteration_L2=50'

sed -e "s/#PREFETCH#/${PREFETCH}/g" Makedef.tmpl > Makedef.${VAR}
export SCALE_SYS=${VAR}
make clean
make depends
make
mkdir -p ${VAR}
mv ./*.lst ./${VAR}/

#####
VAR=L1W2L2S50
PREFETCH='-Kprefetch_strong_L2,prefetch_iteration_L2=50,prefetch_iteration=2'

sed -e "s/#PREFETCH#/${PREFETCH}/g" Makedef.tmpl > Makedef.${VAR}
export SCALE_SYS=${VAR}
make clean
make depends
make
mkdir -p ${VAR}
mv ./*.lst ./${VAR}/

#####
VAR=L1W3L2S50
PREFETCH='-Kprefetch_strong_L2,prefetch_iteration_L2=50,prefetch_iteration=3'

sed -e "s/#PREFETCH#/${PREFETCH}/g" Makedef.tmpl > Makedef.${VAR}
export SCALE_SYS=${VAR}
make clean
make depends
make
mkdir -p ${VAR}
mv ./*.lst ./${VAR}/

#####
VAR=L1S0L2S50
PREFETCH='-Kprefetch_strong_L2,prefetch_iteration_L2=50,prefetch_strong'

sed -e "s/#PREFETCH#/${PREFETCH}/g" Makedef.tmpl > Makedef.${VAR}
export SCALE_SYS=${VAR}
make clean
make depends
make
mkdir -p ${VAR}
mv ./*.lst ./${VAR}/

#####
VAR=L1S2L2S50
PREFETCH='-Kprefetch_strong_L2,prefetch_iteration_L2=50,prefetch_strong,prefetch_iteration=2'

sed -e "s/#PREFETCH#/${PREFETCH}/g" Makedef.tmpl > Makedef.${VAR}
export SCALE_SYS=${VAR}
make clean
make depends
make
mkdir -p ${VAR}
mv ./*.lst ./${VAR}/

#####
VAR=L1S3L2S50
PREFETCH='-Kprefetch_strong_L2,prefetch_iteration_L2=50,prefetch_strong,prefetch_iteration=2'

sed -e "s/#PREFETCH#/${PREFETCH}/g" Makedef.tmpl > Makedef.${VAR}
export SCALE_SYS=${VAR}
make clean
make depends
make
mkdir -p ${VAR}
mv ./*.lst ./${VAR}/
