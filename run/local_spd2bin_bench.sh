#! /bin/bash -x

~/Dropbox/Inbox/scale3/sbin/MacOSX-ifort/spd2bin GRID_DXYZ=20  \
                                                 GRID_KMAX=336 \
                                                 GRID_IMAX=63  \
                                                 GRID_JMAX=63  \
                                                 PRC_nmax=1    \
                                                 PRC_NUM_X=1   \
                                                 PRC_NUM_Y=1   \
                                                 gridfile=grid_20m_336x63x63 \
                                                 step_str=1    \
                                                 step_end=10   \
                                                 selectvar='DENS,U,V,W,PT' \
                                                 history
