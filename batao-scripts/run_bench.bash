#!/bin/bash -e

# What number of machines to benchmark on
nps="2 4"

probs=("100K_whole_city.problem")
solvers=("pbjacobi" "gamg")
flags=("" "-fieldsplit_camera_pc_type gamg")

for np in $nps; do
  for prob in "${probs[@]}"; do
    for i in "${!solvers[@]}"; do
      probname="${prob%.*}"
      d=$(date "+%y-%m-%d")
      set -x
      sbatch -N $np --time=01:00:00 --output "slurm-$d-$probname-${solvers[$i]}-np$np.out" -J "$probname-${solvers[$i]}" batao-scripts/parallel.bash -bal "ba-problems/$prob/problem_noised.bbal" -csv "ba-problems/$prob/parallel_${solvers[$i]}_np$np.csv" -robust ${flags[$i]} -log_view "ascii:ba-problems/$prob/parallel_${solvers[$i]}_np$np.log" "$@"
      set +x
    done
  done
done
