#!/bin/bash -ex

git submodule update --init --recursive

export PETSC_DIR=$(pwd)/petsc
export PETSC_ARCH=arch-opt
export SLEPC_DIR=$(pwd)/slepc

echo "Building PETSc"
cd petsc
./configure --with-debugging=0 COPTFLAGS='-O3 -march=native -mtune=native -g -ggdb' CXXOPTFLAGS='-O3 -march=native -mtune=native -g -ggdb' FOPTFLAGS='-O3 -march=native -mtune=native -g -ggdb' --download-ptscotch --with-hdf5
make
cd ..

echo "Building SLEPc"
cd slepc
./configure
make
cd ..

echo "Building Ceres Solver"
mkdir -p ceres-solver/cmake-build
cmake -S ceres-solver -B ceres-solver/cmake-build -DCMAKE_BUILD_TYPE=Release -DEXPORT_BUILD_DIR=On
cmake --build ceres-solver/cmake-build -- -j $(nproc)

echo "Building BA-TAO"
mkdir -p ba-tao/build
cmake -S ba-tao -B ba-tao/build -DCMAKE_BUILD_TYPE=Release -DCeres_DIR=$(pwd)/ceres-solver/cmake-build
cmake --build ba-tao/build -- -j $(nproc) VERBOSE=1

echo "Building BA-EIG"
mkdir -p ba_eig/build
cmake -S ba_eig -B ba_eig/build -DCMAKE_BUILD_TYPE=Release
cmake --build ba_eig/build -- -j $(nproc) VERBOSE=1

echo "Done!"
