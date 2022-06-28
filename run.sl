#!/bin/sh
#SBATCH -J oasis
#SBATCH -A safe     # Project Account
#SBATCH --time=48:00:00   # Walltime
#SBATCH --mem-per-cpu=10G # memory/cpu 

ml Armadillo/10.5.3-foss-2020b
ml OpenFAST/3.0.0-foss-2020b
/home/projects/safe/oasis/bin/oasis "$1"