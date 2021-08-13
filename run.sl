#!/bin/sh
#SBATCH -J oasis
#SBATCH -A corewind     # Project Account
#SBATCH --time=2:00:00   # Walltime
#SBATCH --mem-per-cpu=4G # memory/cpu 

module purge
ml Armadillo/10.5.3-foss-2020b
ml OpenFAST/3.0.0-foss-2020b
ml OpenMPI/4.0.5-GCC-10.2.0
mpiexec -np 1 ./bin/oasis "$1"