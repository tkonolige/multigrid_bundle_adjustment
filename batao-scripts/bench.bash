#!/bin/bash
#SBATCH -A m1489
#SBATCH -p regular
#SBATCH -C haswell
#SBATCH --ntasks-per-node=4
#SBATCH --hint=memory_bound

srun shifter --image=registry.services.nersc.gov/tkonolig/batao:latest TaoBundleAdjustment "$@"
