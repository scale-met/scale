#! /bin/bash -x

~/Dropbox/Inbox/scale3/sbin/MacOSX-ifort/spd2bin GRID_DXYZ=1000 \
                                                 GRID_KMAX=19   \
                                                 GRID_IMAX=126  \
                                                 GRID_JMAX=126  \
                                                 PRC_nmax=6     \
                                                 PRC_NUM_X=3    \
                                                 PRC_NUM_Y=2    \
                                                 gridfile="grid_1000m_19x126x126" \
                                                 init_coldbubble_63072000000.000
