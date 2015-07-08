#!/bin/bash -x

# mkdir -p html
# for f in $(ls *.tex)
# do
#    ff=`basename ${f}`
#
#    pandoc ${f} --to=html >  ./html/${ff%.tex}.html
# done

pandoc scale_users_guide.tex --to=html >  scale_users_guide.html
