# Multigrid for Bundle Adjustment

This project contains the implementation of a serial multigrid solver for bundle adjustment in addition to a distributed memory Levenberg-Marquardt implementation in PETSc.

## Problem file format

Bundle adjustment problems are stored in either the Bundle Adjustment in the Large format ([http://grail.cs.washington.edu/projects/bal/](http://grail.cs.washington.edu/projects/bal/)), or in a binary version of it.

## Important Directories

- BALUtils: julia code for interacting with bundle adjustment problems.
- ba-problems: test data and results.
- ba-tao: Distributed memory bundle adjuster. Also contains the test harness for the serial solver.
- ba_eig: tools for computing the eigenvalues/eigenvectors of bundle adjustment problems.
- bamg: the multigrid solver (written in Julia).
- ceres-solver: a fork of Ceres Solver with hooks to use the serial multigrid solver.
- city2ba: tools for generating synthetic problems.
- models: 3D city models used to generate synthetic problems.
- multigrid_analysis: tools to apply compatible relaxation analysis to the multigrid solver.
- notebooks: jupyter notebooks for result analysis and plotting.
- options: various configuration used for testing bundle adjustment solvers.
- petsc: fork of PETSc with a parallel Levenberg-Marquardt implementation.

## How to Build

This project requires:

- cmake
- [julia](https://julialang.org)
- [embree](https://www.embree.org)
- [ghc](https://www.haskell.org/ghc/)/[cabal](https://www.haskell.org/cabal/)
- [rust](https://www.rust-lang.org)

With these installed, run `setup.bash` to build everything.

## Collecting Results

`Build.hs` is a shake build file that can build everything and run tests. It requires `ghc` (the Haskell compiler) and the `shake` package to be installed (you can install it with `cabal install shake`). There are a variety of targets that it can build:

- ba-problems/*/ceres_benchmark_pbjacobi.csv: collect results for the point block Jacobi preconditioner.
- ba-problems/*/ceres_benchmark_visibility.csv: collect results for the visibility preconditioner.
- ba-problems/*/ceres_benchmark_multigrid.csv: collect results for the multigrid preconditioner.

You can list all targets by running `./build.sh -h`.
