#!/bin/bash
DNUM=1
while [ $DNUM -le $NUM_DOMAIN ]
do
  D=`expr $DNUM - 1`

  FNUM=`printf "%02d" $DNUM`

cat << EOF > ${OUTPUT_CONFIGDIR}/pp/Makefile
################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = \$(shell pwd)
TOPDIR      = \$(abspath ../../../../../..)
TESTDIR     = ../../../..


# user-defined source files
CODE_DIR    = .
ORG_SRCS    =

# parameters for run
INITCONF    = pp.d${FNUM}.conf

TPROC       = `expr ${PRC_NUM_X[$D]} \* ${PRC_NUM_Y[$D]}`

# required data (parameters,distributed files)
DATDIR      =
DATPARAM    =
DATDISTS    =



all: run

# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include \$(TESTDIR)/Makefile.tutorial.common

EOF

cat << EOF > ${OUTPUT_CONFIGDIR}/init/Makefile
################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = \$(shell pwd)
TOPDIR      = \$(abspath ../../../../../..)
TESTDIR     = ../../../..


# user-defined source files
CODE_DIR    = .
ORG_SRCS    =

# parameters for run
INITCONF    = init.d${FNUM}.conf

TPROC       = 4

# required data (parameters,distributed files)
DATDIR      = ../../data
DATPARAM    = ${BASENAME_ORG}
DATDISTS    =



all: run

# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include \$(TESTDIR)/Makefile.tutorial.common

EOF

cat << EOF > ${OUTPUT_CONFIGDIR}/run/Makefile
################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = \$(shell pwd)
TOPDIR      = \$(abspath ../../../../../..)
TESTDIR     = ../../../..


# user-defined source files
CODE_DIR    = .
ORG_SRCS    =

# parameters for run
RUNCONF     = run.d${FNUM}.conf

TPROC       = `expr ${PRC_NUM_X[$D]} \* ${PRC_NUM_Y[$D]}`

# required data (parameters,distributed files)
DATDIR      =
DATPARAM    =
DATDISTS    =



all: run

# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include \$(TESTDIR)/Makefile.tutorial.common

EOF

cat << EOF > ${OUTPUT_CONFIGDIR}/net2g/Makefile
################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = \$(shell pwd)
TOPDIR      = \$(abspath ../../../../../..)
TESTDIR     = ../../../..


# user-defined source files
CODE_DIR    = .
ORG_SRCS    =

# parameters for run


N2GCONF     = net2g.2D.d${FNUM}.conf #,net2g.3D.d${FNUM}.conf
TPROC       = `expr ${PRC_NUM_X[$D]} \* ${PRC_NUM_Y[$D]}`

# required data (parameters,distributed files)
DATDIR      =
DATPARAM    =
DATDISTS    =



all: run

# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include \$(TESTDIR)/Makefile.tutorial.common

EOF

  DNUM=`expr $DNUM + 1`
done

