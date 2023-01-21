#!/bin/bash

#SBATCH --array=0-28
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4 --mem-per-cpu=20G 
#SBATCH --time=14-0

ml R/4.1.2-foss-2020b

currDir='/home/eulloape/CATE'
cd $currDir

S=(11 11 11 11)
ss=(1000 2000 5000 10000)

srun R --vanilla "--args nall=${ss[${SLURM_ARRAY_TASK_ID}]} ntest=10000 scen=${S[${SLURM_ARRAY_TASK_ID}]}" <main-script-alg2.R> job.${SLURM_ARRAY_TASK_ID}.Rout





