#!/bin/sh
#SBATCH -J oasis_morison_wind_drag_yaml
#SBATCH -A oasys     # Project Account
#SBATCH --time=8:00:00   # Walltime
#SBATCH --mem-per-cpu=32G # memory/cpu

module purge
ml Armadillo/10.5.3-foss-2020b
ml FASTurbine_wrapper/0.1.1-foss-2020b
../../../bin/oasis "$1"
rm -f core*
