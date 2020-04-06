find_package(PkgConfig REQUIRED)
if(DEFINED ENV{PETSC_DIR})
  set(ENV{PKG_CONFIG_PATH} "$ENV{PKG_CONFIG_PATH}:$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/lib/pkgconfig" )
endif()
pkg_check_modules(PETSC IMPORTED_TARGET PETSc)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (PETSC
  "PETSc could not be found."
  PETSC_FOUND)