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
PPCONF      = ${LIST_PP_CONF_FILES}

TPROC       = ${LIST_PRC_DOMAINS}

# required data (parameters,distributed files)
DATPARAM    =
DATDISTS    =


# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include \$(TESTDIR)/Makefile.common


all: run

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
INITCONF    = ${INIT_CONF_FILE}

TPROC       = ${TPROC}

# required data (parameters,distributed files)
DATPARAM    = ${DATPARAM}
DATDISTS    = ${DATDISTS}


# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include \$(TESTDIR)/Makefile.common


all: run

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
RUNCONF     = ${RUN_CONF_FILE}

TPROC       = ${TPROC}

# required data (parameters,distributed files)
DATPARAM    =
DATDISTS    =


# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include \$(TESTDIR)/Makefile.common


all: run

EOF

cat << EOF > ${OUTPUT_CONFIGDIR}/sno/Makefile
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
SNOCONF     = ${LIST_SNO_CONF_FILES}

TPROC       = ${PRC_SNO}

# required data (parameters,distributed files)
DATDIR      =
DATPARAM    =
DATDISTS    =

# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include \$(TESTDIR)/Makefile.common


all: run

EOF

  DNUM=`expr $DNUM + 1`
done

