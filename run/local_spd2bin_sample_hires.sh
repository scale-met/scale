#! /bin/bash -x

~/Dropbox/Inbox/scale3/sbin/MacOSX-ifort/spd2bin GRID_DXYZ=40  \
                                                 GRID_IMAX=70  \
                                                 GRID_JMAX=100 \
                                                 GRID_KMAX=220 \
                                                 PRC_nmax=1    \
                                                 PRC_NUM_X=1   \
                                                 PRC_NUM_Y=1   \
                                                 gridfile=grid_40m_70x100x220 \
                                                 history
