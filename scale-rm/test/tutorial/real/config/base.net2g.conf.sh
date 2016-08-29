#!/bin/bash

str=""; for v in ${POPSCA_PLEV[*]}; do str="${str}${v}, "; done; LIST_POPSCA_PLEV=${str%, }

str=""; for v in ${POPSCA_2D[*]}; do str="${str}\"${v}\", "; done; LIST_POPSCA_2D=${str%, }
str=""; for v in ${POPSCA_3D[*]}; do str="${str}\"${v}\", "; done; LIST_POPSCA_3D=${str%, }

cat << EOF > base.net2g.2D.conf

#################################################
#
# post process configuration: POPSCA
#
#################################################

&LOGOUT
 LOG_BASENAME   = "${NET2G_2D_IO_LOG_BASENAME}",
 LOG_ALL_OUTPUT = .false.,
 LOG_LEVEL      = 1,
/

&INFO
 TIME_STARTDATE = ${TIME_STARTDATE},
 START_TSTEP    = 1,
 END_TSTEP      = ${MAXSTEP_2D},
 INC_TSTEP      = 1,
 DOMAIN_NUM     = ${DNUM},
 ZSTART         = 1,
 ZCOUNT         = 1,
 CONFFILE       = "./${RUNDIR}/run.d${FNUM}.conf",
 IDIR           = "./${RUNDIR}",
 ODIR           = ".",
 Z_LEV_TYPE     = "original",
 Z_MERGE_OUT    = .false.,
 T_MERGE_OUT    = .true.,
/

&VARI
 VNAME = ${LIST_POPSCA_2D},
/
EOF

cat << EOF > base.net2g.3D.conf

#################################################
#
# post process configuration: POPSCA
#
#################################################

&LOGOUT
 LOG_BASENAME   = "${NET2G_3D_IO_LOG_BASENAME}",
 LOG_ALL_OUTPUT = .false.,
 LOG_LEVEL      = 1,
/

&INFO
 TIME_STARTDATE = ${TIME_STARTDATE},
 START_TSTEP    = 1,
 END_TSTEP      = ${MAXSTEP_3D},
 INC_TSTEP      = 1,
 DOMAIN_NUM     = ${DNUM},
 ZSTART         = 1,
 ZCOUNT         = ${#POPSCA_PLEV[*]},
 CONFFILE       = "./${RUNDIR}/run.d${FNUM}.conf",
 IDIR           = "./${RUNDIR}",
 ODIR           = ".",
 Z_LEV_TYPE     = "plev",
 Z_MERGE_OUT    = .true.,
 T_MERGE_OUT    = .true.,
/

&EXTRA
 EXTRA_TINTERVAL = ${TIME_DT_HISTORY_3D},
 EXTRA_TUNIT     = "${TIME_DT_UNIT}",
/

&VARI
 VNAME       = ${LIST_POPSCA_3D},
 TARGET_ZLEV = ${LIST_POPSCA_PLEV},
/
EOF
