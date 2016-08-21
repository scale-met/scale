#!/bin/bash

cat << EOF > base.init.launch.conf

#################################################
#
# launcher configuration
#
#################################################

&PARAM_LAUNCHER
 NUM_DOMAIN  = ${NUM_DOMAIN},
 PRC_DOMAINS = ${LIST_PRC_DOMAINS},
 CONF_FILES  = "${LIST_INIT_CONF_FILES}",
/
EOF

cat << EOF > base.run.launch.conf

#################################################
#
# launcher configuration
#
#################################################

&PARAM_LAUNCHER
 NUM_DOMAIN  = ${NUM_DOMAIN},
 PRC_DOMAINS = ${LIST_PRC_DOMAINS},
 CONF_FILES  = "${LIST_RUN_CONF_FILES}",
/
EOF
