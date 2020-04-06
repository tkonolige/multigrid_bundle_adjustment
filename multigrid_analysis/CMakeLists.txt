cmake_minimum_required( VERSION 3.12 )

project( Multigrid_Analysis
         VERSION 3.12
         DESCRIPTION "Various techniques for analyzing Multigrid performance." )

# Put FindPackage modules and others in cmake/
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

find_package(PETSc REQUIRED)
find_package(SLEPc REQUIRED)
# cmake doesn't include /usr/local/include in its pkgconfig file, cmake doesn't look there either
find_package(HDF5 REQUIRED)

add_executable(cr_analysis)
target_sources(cr_analysis PRIVATE "cr-analysis.cpp")
target_link_libraries(cr_analysis PRIVATE PETSc::PETSc SLEPc::SLEPc)
target_include_directories(cr_analysis PRIVATE ${HDF5_C_INCLUDE_DIRS})
