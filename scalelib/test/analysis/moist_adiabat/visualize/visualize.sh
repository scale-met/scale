#!/bin/bash -x

echo "Skew-T log-P diagram"
ruby visualize/skewTlogP.rb output.dat
mv dcl_0001.png skewTlogP.png
