#!/bin/usr/env bash
export R_LIBS=~/Rlibs2
export R_LIBS_USER=~/Rlibs2
for S in 4 4 4 4
do
  for ss in 1000 2000 5000 10000
  do
    sbatch  --export=s=$s,ss=ss ~/calibration/scenario2.sbatch
  done
done


