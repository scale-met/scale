#! /bin/bash -x

~/Dropbox/Inbox/scale3/sbin/MacOSX-ifort/spd2bin GRID_DXYZ=500 \
                                                 GRID_IMAX=120 \
                                                 GRID_JMAX=120 \
                                                 GRID_KMAX=25  \
                                                 PRC_nmax=1    \
                                                 PRC_NUM_X=1   \
                                                 PRC_NUM_Y=1   \
                                                 gridfile=grid_500m_120x120x25 \
                                                 history
