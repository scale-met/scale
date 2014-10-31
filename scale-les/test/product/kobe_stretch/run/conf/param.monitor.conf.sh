#!/bin/bash

cat << EOF > conf/param.monitor.conf

#################################################
#
# model configuration: monitor
#
#################################################

&PARAM_MONITOR
 MONITOR_OUT_BASENAME  = "${MONITOR_OUT_BASENAME}",
 MONITOR_STEP_INTERVAL = 120,
/

&MONITITEM item="QDRY" /
&MONITITEM item="QTOT" /
&MONITITEM item="ENGT" /
&MONITITEM item="ENGP" /
&MONITITEM item="ENGK" /
&MONITITEM item="ENGI" /

&MONITITEM item="ENGFLXT" /

&MONITITEM item="ENGSFC_SH" /
&MONITITEM item="ENGSFC_LH" /
&MONITITEM item="ENGSFC_RD" /
&MONITITEM item="ENGTOA_RD" /
EOF
