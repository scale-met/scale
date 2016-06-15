#!/bin/sh
#-----------------------------------------
#   Test Script for XXX  (Ver.0.1)
#   2009/10/26 --- Ryuji Yoshida.
#   2016/06/15  modified  S.Adachi
#-----------------------------------------
#
FNAME="scale_users_guide"
#
#---- Replace Version
scale_version=`gawk '{print $0}' ../../../src/VERSION`
#echo "SCALE VERSION" $scale_version
#sed 3s/_VERSION_/$scale_version/ template/scale_users_guide.tex > scale_users_guide.tex
#sed 4s/_VERSION_/$scale_version/ template/11_introduction.tex   > 11_introduction.tex
#
#
#
# first time
platex ${FNAME}.tex
# bibtex first time
bibtex ${FNAME}
# second time
platex ${FNAME}.tex
# bibtex second time
bibtex ${FNAME}
# final time
platex ${FNAME}.tex
#
# dvi to pdf
dvipdfmx -p a4 ${FNAME}.dvi
#
if [ -e /etc/vine-release ]; then
 evince ${FNAME}.pdf
fi
#
#-----------------------------------------
#EOF
