#! /bin/bash -x

~/Dropbox/Inbox/scale3/sbin/MacOSX-ifort/spd2bin GRID_IMAX=100 \
                                                 GRID_JMAX=100 \
                                                 GRID_KMAX=40  \
                                                 GRID_DX=200   \
                                                 GRID_DY=200   \
                                                 GRID_DZ=200   \
                                                 PRC_nmax=4    \
                                                 PRC_NUM_X=2   \
                                                 PRC_NUM_Y=2   \
                                                 restart_out_63072000030.000