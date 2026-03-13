#!/bin/sh
#SBATCH -J oasis_elastic_anchor
#SBATCH -A oasys     # Project Account
#SBATCH --time=8:00:00   # Walltime
#SBATCH --mem-per-cpu=4G # memory/cpu

module purge
ml Armadillo/10.5.3-foss-2020b
ml FASTurbine_wrapper/0.1.1-foss-2020b
../../../bin/oasis "$1"
rm -f core*
