#!/bin/sh
#-----------------------------------------
#   Test Script for XXX  (Ver.0.1)
#   2009/10/26 --- Ryuji Yoshida.
#-----------------------------------------
#
FNAME="scale_users_guide"
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
