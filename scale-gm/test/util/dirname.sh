#! /bin/bash -x

GLEV=${1:-5}
RLEV=${2:-0}
NMPI=${3:-5}
GL=`printf %02d ${GLEV}`
RL=`printf %02d ${RLEV}`
if   [ ${NMPI} -ge 10000 ]; then
	NP=`printf %05d ${NMPI}`
elif [ ${NMPI} -ge 1000 ]; then
	NP=`printf %04d ${NMPI}`
elif [ ${NMPI} -ge 100 ]; then
	NP=`printf %03d ${NMPI}`
else
	NP=`printf %02d ${NMPI}`
fi
ZL=${4:-40}
output_selector=${5:-"dir"}


if [ ${output_selector} == "dir" ]; then
   dir=gl${GL}rl${RL}z${ZL}pe${NP}
   echo ${dir}
elif [ ${output_selector} == "res" ]; then
   res=GL${GL}RL${RL}z${ZL}
   echo ${res}
fi