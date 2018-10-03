################################################################################
#
# Makefile for each test program
#
################################################################################

PWD         = $(shell pwd)
TOPDIR      = $(abspath ../../../../..)
TESTDIR     = ../../..

ifndef PREPDIR
   PREPDIR = ../pp
endif

ifndef INITDIR
   INITDIR = ../init_FNL
endif

ifndef HISTDIR
   HISTDIR = ../run_FNL
endif

# user-defined source files
CODE_DIR    = .
ORG_SRCS    =

# parameters for run
INITCONF    = init_d02_history.conf,init_d02_restart.conf

TPROC       = 6,6

# required data (parameters,distributed files)
DATPARAM    = " [${TESTDIR}/data/land/param.bucket.conf param.bucket.conf ] \
                [$(PREPDIR)/latlon_domain_catalogue_d01.txt latlon_domain_catalogue_d01.txt] "
DATDISTS    = " [6 $(PREPDIR)/topo_d02                     topo_d02                    ] \
                [6 $(PREPDIR)/landuse_d02                  landuse_d02                 ] \
                [3 $(INITDIR)/init_d01_20100501-000000.000 init_d01_20100501-000000.000] \
                [3 $(HISTDIR)/history_d01_5min             history_d01_5min            ] "



# build, makedir, run, jobshell, allclean, clean is inside of common Makefile
include $(TESTDIR)/Makefile.common

# special execute
run: jobshell
	$(JOBSUB) run.sh
