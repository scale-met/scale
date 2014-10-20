#!/bin/bash

cat << EOF > conf/param.region.conf

#################################################
#
# model configuration: process
#
#################################################

&PARAM_PRC
 PRC_NUM_X      = 3,
 PRC_NUM_Y      = 3,
 PRC_PERIODIC_X = .false.,
 PRC_PERIODIC_Y = .false.,
/

#################################################
#
# model configuration: region
#
#################################################

&PARAM_INDEX
 KMAX = 30,
 IMAX = 20,
 JMAX = 20,
/

&PARAM_LAND_INDEX
 LKMAX = 5,
/

&PARAM_URBAN_INDEX
 UKMAX = 4,
/

&PARAM_LAND_GRID
 LDZ = 0.05D0, 0.15D0, 0.30D0, 0.50D0, 1.00D0,
/

&PARAM_URBAN_GRID
 UDZ = 0.05D0, 0.15D0, 0.30D0, 0.50D0,
/

&PARAM_GRID
 DX = 2500.D0,
 DY = 2500.D0,
 FZ(:) =    500.0000D0,  1000.0000D0,  1500.0000D0,  2000.0000D0,  2270.0000D0,
           2576.4500D0,  2924.2708D0,  3319.0474D0,  3767.1187D0,  4275.6797D0,
           4852.8965D0,  5508.0376D0,  6251.6226D0,  7095.5918D0,  7995.5918D0,
           8895.5918D0,  9795.5918D0, 10695.5918D0, 11595.5918D0, 12495.5918D0,
          13083.2021D0, 13656.1221D0, 14214.7188D0, 14759.3506D0, 15290.3672D0,
          15808.1084D0, 16312.9062D0, 16805.0840D0, 17284.9570D0, 17752.8340D0,
 BUFFER_DZ =  5000.D0,
 BUFFER_DX = 10000.D0,
 BUFFER_DY = 10000.D0,
/

&PARAM_MAPPROJ
 MPRJ_basepoint_lon = 134.25D0,
 MPRJ_basepoint_lat =  34.00D0,
/
EOF