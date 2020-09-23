#! /bin/bash -x

### Lat-Lon conversion ###
echo "+convert by sno"

TOPDIR=`pwd`/../../../../../..

cd ..
${TOPDIR}/bin/sno convert/sno_d01_latlon.conf || exit 1
${TOPDIR}/bin/sno convert/sno_d02_latlon.conf || exit 1
