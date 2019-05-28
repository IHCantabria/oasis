#!/bin/bash
#SBATCH -J SEMLines
#SBATCH -A telwind        # Project Account
#SBATCH --time=50:00:00   # Walltime
#SBATCH --mem-per-cpu=1GB  # memory/cpu 
# -------------------------------------------------------

# Cargo el modulo de Openblas
ml OpenBLAS/0.2.19-GCC-6.3.0-2.27-LAPACK-3.7.0

# Doy permisos al ejecutable
chmod 777 ../../src/SEM_Lines.exe

# Ejecuto el caso
../../src/SEM_Lines.exe
