#!/bin/usr/env bash
export R_LIBS=~/Rlibs2
export R_LIBS_USER=~/Rlibs2
for S in 4 11
do
  for ss in 1000 2000 5000 10000
  do
    sbatch  --export=S=$S,ss=$ss ~/calibration/scenario.sbatch
  done
done


