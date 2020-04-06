find_package(PkgConfig REQUIRED)
find_package(MPI REQUIRED)

if(DEFINED ENV{PETSC_DIR})
  set(ENV{PKG_CONFIG_PATH} "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}" )
  set(ENV{PKG_CONFIG_PATH} "$ENV{PETSC_DIR}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}" )
endif()
pkg_check_modules(PETSC IMPORTED_TARGET GLOBAL PETSc)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (PETSc
  FAIL_MESSAGE "PETSc could not be found. Please set PETSC_ARCH and PETSC_DIR in your environment."
  REQUIRED_VARS PETSC_FOUND)

if(PETSC_FOUND)
  add_library(PETSc::PETSc INTERFACE IMPORTED GLOBAL)
  target_link_libraries(PETSc::PETSc INTERFACE MPI::MPI_C PkgConfig::PETSC)
endif()
