#!/bin/bash

#SBATCH --array=1
#SBATCH --nodes=1
#SBATCH --output=Rout/scenario2-%J.out
#SBATCH --error=Rout/scenario2-%J.err
#SBATCH --cpus-per-task=1
echo "LOADING R"
module load R
echo "R LOADED"
Rscript -e 'source("R_setup.R"); S = Sys.getenv("S"); ss = Sys.getenv("ss"); source("~/calibration/main-script-alg2.R")'
