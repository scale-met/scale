#!/bin/bash

cat << EOF > base.launch.conf

#################################################
#
# launcher configuration
#
#################################################

&PARAM_LAUNCHER
 NUM_DOMAIN  = ${NUM_DOMAIN},
 PRC_DOMAINS = ${LIST_PRC_DOMAINS},
 CONF_FILES  = ${LIST_CONF_FILES},
/
EOF
