#!/bin/sh
################################################################################
#
# ---- Preparation script for easy start of projects.
# ---- H.Yashiro 11/08/19
#
################################################################################

projname=${1}

if [ -z ${projname} ]; then
   echo "usage: sh mkproject.sh projectname "
   exit 1
fi

if [ ! -d ${projname} ]; then
   mkdir -p ./${projname}
fi

if [ ! -f ${projname}/${projname}.f90 ]; then
   cp ../scale3.f90 ./${projname}/${projname}.f90
fi

sed -e "s/#PROJECT#/${projname}/g" Makefile.tmpl > ./${projname}/Makefile

################################################################################
