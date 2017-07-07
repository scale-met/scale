#! /bin/bash -x

ruby ./evaluate/diff.rb restartA_d01_*.pe000000.nc restartB2_d01_*.pe000000.nc
ruby ./evaluate/diff.rb restartA_d02_*.pe000000.nc restartB2_d02_*.pe000000.nc
