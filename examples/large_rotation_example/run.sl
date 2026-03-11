#!/bin/sh
#SBATCH -J oasis_merge
#SBATCH -A oasys     # Project Account
#SBATCH --time=8:00:00   # Walltime
#SBATCH --mem-per-cpu=4G # memory/cpu

module purge
ml Armadillo/10.5.3-foss-2020b
ml FASTurbine_wrapper/0.1.1-foss-2020b

echo "====== Case 1: Simplified rotation (rotSimpFlag=1) ======"
../../bin/oasis case_simplified
RESULT1=$?

echo "====== Case 2: Full rotation dynamics (rotSimpFlag=0) ======"
../../bin/oasis case_full
RESULT2=$?

rm -f core*

if [ $RESULT1 -ne 0 ] || [ $RESULT2 -ne 0 ]; then
    echo "FAILED: one or both cases failed"
    exit 1
fi
echo "Both cases completed successfully."
