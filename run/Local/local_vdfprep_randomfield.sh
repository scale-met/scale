#! /bin/bash -x

export HMDIR=~/GCMresults/sol/latest
export BIN=~/Dropbox/Inbox/scale3/sbin/${SCALE_SYS}
export EXE=vdfprep

export OUTDIR=${HMDIR}/output/randomfield_check

cd ${OUTDIR}

${BIN}/${EXE}