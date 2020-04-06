#!/bin/bash
#SBATCH -A m1489
#SBATCH -p regular
#SBATCH -C haswell
#SBATCH --ntasks-per-node=1
#SBATCH -N 1
#SBATCH --hint=memory_bound
#SBATCH --time 04:00:00

export JULIA_DEPOT_PATH="$SCRATCH/julia-pkg"
export PATH="$SCRATCH/julia/bin:$PATH"
echo "$@"
$SCRATCH/bundle_adjustment/ba-tao/build/bin/CeresBundleAdjustment "$@"
