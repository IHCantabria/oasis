#!/bin/sh
#SBATCH -J oasis
#SBATCH -A corewind     # Project Account
#SBATCH --time=2:00:00   # Walltime
#SBATCH --mem-per-cpu=4G # memory/cpu 

module purge
ml Armadillo/10.5.3-foss-2020b
ml FASTurbine_wrapper/0.0-foss-2020b
/home/projects/energia/developments/oasis_wind_turbine/oasis/bin/oasis "$1"