#!/bin/sh
#SBATCH -J oasis
#SBATCH -A corewind     # Project Account
#SBATCH --time=2:00:00   # Walltime
#SBATCH --mem-per-cpu=4G # memory/cpu 

ml Armadillo/10.5.3-foss-2020b
ml OpenFAST/3.0.0-foss-2020b
./bin/oasis "$1"