#!/bin/sh
#-----------------------------------------
#   Test Script for XXX  (Ver.0.1)
#   2009/10/26 --- Ryuji Yoshida.
#-----------------------------------------
# first time
platex scale_users_guide.tex
# second time
platex scale_users_guide.tex
# bibtex
bibtex scale_users_guide
# final time
platex scale_users_guide.tex
# dvi to pdf
dvipdfmx scale_users_guide.dvi
#
#acroread scale_users_guide.pdf
#
#-----------------------------------------
#EOF
