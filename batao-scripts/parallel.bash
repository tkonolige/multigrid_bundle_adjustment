#!/bin/bash
#SBATCH -A m1489
#SBATCH -p regular
#SBATCH -C haswell
#SBATCH --ntasks-per-node=4
#SBATCH --hint=memory_bound
#SBATCH --mail-type=END,FAIL

export JULIA_DEPOT_PATH="$SCRATCH/julia-pkg"
export PATH="$SCRATCH/julia/bin:$PATH"
echo "$@"
srun $SCRATCH/bundle_adjustment/ba-tao/build/bin/TaoBundleAdjustment -options_file $SCRATCH/bundle_adjustment/ba-tao/default_options.txt "$@"
