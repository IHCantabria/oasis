#!/bin/sh
#SBATCH -J oasis
#SBATCH -A fconcrete     # Project Account
#SBATCH --time=4:00:00   # Walltime
#SBATCH --mem-per-cpu=10G # memory/cpu 

ml HDF5/1.8.19-foss-2017a
ml Armadillo/9.800.4-foss-2017a
ml CMake/3.7.2-foss-2017a
chmod 777 ./bin/oasis
./bin/oasis "$1"