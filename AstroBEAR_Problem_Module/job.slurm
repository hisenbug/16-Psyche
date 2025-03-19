#!/bin/bash

###=================================================
## AstroBEAR job submission script on bluehive2
## can be submitted with: sbatch job.slurm
##  by Baowei Liu 04.16.2014
##  modified 03.19.2025 by Atma Anand
## =================================================

## Reserve multiple cores/RAM/time
# SBATCH -N 3 --mem-per-cpu=2000m --ntasks-per-node=36 -p standard -t 5-00:00:00 --mail-type=END,FAIL,BEGIN

## Use Single core for making output of points on the 10 km surface layer (Uncomment this, comment previous)
## SBATCH -N 1 --ntasks-per-node=1 -p debug -t 01:00:00

source /software/modules/init/bash

module purge
module load astrobear/b1
module load slurm

ulimit -s unlimited

mkdir -p out
## IMPORTANT: Enter path to your AstroBEAR Makefile location on your local disk in this line
## Or the directory should already contained the compiled "astrobear"

mpirun -np $SLURM_NTASKS ./astrobear
