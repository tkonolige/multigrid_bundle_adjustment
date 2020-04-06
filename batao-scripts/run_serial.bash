#!/bin/bash

#probs=("26K_29K_whole_city_hard.problem" "26K_29K_whole_city.problem" "17K_25K_whole_city.problem" "15K_6K_inner_city.problem")
probs=("5K_long_path.problem" "10K_long_path.problem" "15K_long_path.problem" "20K_long_path.problem" "25K_long_path.problem" "30K_long_path.problem")
# probs=("26K_29K_whole_city_hard.problem" "26K_29K_whole_city" "17K_25K_whole_city")
# solvers=("visibility")
# flags=("--cluster_jacobi")
solvers=("pbjacobi" "multigrid" "visibility")
flags=("" "--multigrid" "--cluster_jacobi")

for prob in "${probs[@]}"; do
  for i in "${!solvers[@]}"; do
    set -x
    sbatch --time 02:00:00 -J "$prob-${solvers[$i]}" --mail-type=END,FAIL batao-scripts/serial.bash --bal "ba-problems/$prob/problem_noised.bbal" "${flags[$i]}" --csv "ba-problems/$prob/ceres_benchmark_${solvers[$i]}_robust01.csv" "$@"
    set +x
  done
done
