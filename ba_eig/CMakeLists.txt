cmake_minimum_required(VERSION 3.1)

project(BA_Eig C CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

find_package(MPI REQUIRED)
find_package(PETSc REQUIRED)
find_package(SLEPc REQUIRED)
find_package(HDF5 REQUIRED)
find_package(Julia REQUIRED)
find_package(LAPACK REQUIRED)

add_executable(ba_eig ba_eig.cpp)
target_compile_features(ba_eig PRIVATE cxx_std_11)
target_link_libraries(ba_eig PRIVATE PkgConfig::SLEPC PkgConfig::PETSC MPI::MPI_CXX "${Julia_LIBRARIES}" "${LAPACK_LIBRARIES}")
target_include_directories(ba_eig PRIVATE "${HDF5_C_INCLUDE_DIRS}" "${Julia_INCLUDE_DIRS}")

install(TARGETS ba_eig # bundle_adjuster openmvg2bal
        RUNTIME DESTINATION bin
        )
