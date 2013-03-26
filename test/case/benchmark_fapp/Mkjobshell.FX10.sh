#! /bin/bash -x

# Arguments
BINDIR=${1}
INITNAME=${2}
BINNAME=${3}
INITCONF=${4}
RUNCONF=${5}
TPROC=${6}

# System specific
MPIEXEC="mpiexec"

array=( `echo ${TPROC} | tr -s 'x' ' '`)
x=${array[0]}
y=${array[1]:-1}
let xy="${x} * ${y}"

if [ ${xy} -gt 480 ]; then
   rscgrp="x-large"
elif [ ${xy} -gt 372 ]; then
   rscgrp="large"
elif [ ${xy} -gt 216 ]; then
   rscgrp="medium"
elif [ ${xy} -gt 12 ]; then
   rscgrp="small"
else
   rscgrp="short"
fi

# Generate run.sh

cat << EOF1 > ./run.sh
#! /bin/bash -x
################################################################################
#
# for OAKLEAF-FX
#
################################################################################
#PJM --rsc-list "rscgrp=${rscgrp}"
#PJM --rsc-list "node=${TPROC}"
#PJM --rsc-list "elapse=00:30:00"
#PJM -j
#PJM -s
export PARALLEL=8
export OMP_NUM_THREADS=8

ln -sv ${BINDIR}/${INITNAME} .
ln -sv ${BINDIR}/${BINNAME}  .

rm -rf ./pa[1-7]
mkdir -p ./pa1
mkdir -p ./pa2
mkdir -p ./pa3
mkdir -p ./pa4
mkdir -p ./pa5
mkdir -p ./pa6
mkdir -p ./pa7

# run
                      ${MPIEXEC} ./${INITNAME} ${INITCONF} || exit
fapp -C -d pa1 -Hpa=1 ${MPIEXEC} ./${BINNAME}  ${RUNCONF}  || exit
fapp -C -d pa2 -Hpa=2 ${MPIEXEC} ./${BINNAME}  ${RUNCONF}  || exit
fapp -C -d pa3 -Hpa=3 ${MPIEXEC} ./${BINNAME}  ${RUNCONF}  || exit
fapp -C -d pa4 -Hpa=4 ${MPIEXEC} ./${BINNAME}  ${RUNCONF}  || exit
fapp -C -d pa5 -Hpa=5 ${MPIEXEC} ./${BINNAME}  ${RUNCONF}  || exit
fapp -C -d pa6 -Hpa=6 ${MPIEXEC} ./${BINNAME}  ${RUNCONF}  || exit
fapp -C -d pa7 -Hpa=7 ${MPIEXEC} ./${BINNAME}  ${RUNCONF}  || exit

################################################################################
EOF1

cat << EOF2 > ./fapppx.sh
#! /bin/bash -x
################################################################################
#
# for K computer
#
################################################################################

fapppx -A -p0,1 -l0 -d pa1 -o output_prof_1.csv -tcsv -Hpa
fapppx -A -p0,1 -l0 -d pa2 -o output_prof_2.csv -tcsv -Hpa
fapppx -A -p0,1 -l0 -d pa3 -o output_prof_3.csv -tcsv -Hpa
fapppx -A -p0,1 -l0 -d pa4 -o output_prof_4.csv -tcsv -Hpa
fapppx -A -p0,1 -l0 -d pa5 -o output_prof_5.csv -tcsv -Hpa
fapppx -A -p0,1 -l0 -d pa6 -o output_prof_6.csv -tcsv -Hpa
fapppx -A -p0,1 -l0 -d pa7 -o output_prof_7.csv -tcsv -Hpa

################################################################################
EOF2
