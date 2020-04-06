#pragma once

#include "julia.h"

void check_error() {
  if (jl_exception_occurred()) {
    const char *p = (const char *)jl_unbox_voidpointer(jl_eval_string("pointer(sprint(showerror, ccall(:jl_exception_occurred, Any, ())))"));

    throw std::runtime_error(p);
  }
}

jl_value_t* eval_string(const std::string& str) {
  auto ret = jl_eval_string(str.data());
  check_error();
  return ret;
}

jl_function_t* get_function(const std::string& fnname, const std::string& module) {
  auto mod = (jl_module_t*)jl_get_global(jl_main_module, jl_symbol(module.c_str()));
  if(mod == nullptr) {
    throw("Could not get module " + module);
  }
  check_error();
  jl_function_t* func = jl_get_function(mod, fnname.c_str());
  if(func == nullptr) {
    throw("Could not get function " + fnname + " from module " + module);
  }
  return func;
}

jl_value_t* wrap_array(std::vector<double>& vec) {
  jl_value_t* array_type = jl_apply_array_type((jl_value_t*)jl_float64_type, 1);
  return (jl_value_t*)jl_ptr_to_array_1d(array_type, vec.data(), vec.size(), 0);
}

jl_value_t* wrap_array(std::vector<int32_t>& vec) {
  jl_value_t* array_type = jl_apply_array_type((jl_value_t*)jl_int32_type, 1);
  return (jl_value_t*)jl_ptr_to_array_1d(array_type, vec.data(), vec.size(), 0);
}

jl_value_t* wrap_array(const int32_t* vec, int len) {
  jl_value_t* array_type = jl_apply_array_type((jl_value_t*)jl_int32_type, 1);
  return (jl_value_t*)jl_ptr_to_array_1d(array_type, (void*)vec, len, 0);
}

jl_value_t* wrap_array(const double* vec, int len) {
  jl_value_t* array_type = jl_apply_array_type((jl_value_t*)jl_float64_type, 1);
  return (jl_value_t*)jl_ptr_to_array_1d(array_type, (void*)vec, len, 0);
}

PetscErrorCode bamg_apply(PC pc, Vec x, Vec y) {
  PetscErrorCode ierr;
  PetscFunctionBegin;

  PetscInt n;
  ierr = VecGetSize(x, &n);CHKERRQ(ierr);
  const PetscScalar *x_ary;
  PetscScalar *y_ary;
  ierr = VecGetArrayRead(x, &x_ary);CHKERRQ(ierr);
  ierr = VecGetArray(y, &y_ary);CHKERRQ(ierr);

  auto mg = jl_get_global(jl_main_module, jl_symbol("mg"));
  auto lmul = get_function("ldiv!", "LinearAlgebra");
  jl_call3(lmul, wrap_array(y_ary, n), mg, wrap_array(x_ary, n));
  check_error();

  ierr = VecRestoreArrayRead(x, &x_ary);CHKERRQ(ierr);
  ierr = VecRestoreArray(y, &y_ary);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

void bamg_create(PC pc, Mat A, const double* solution, const double* scale, std::string bafile) {
  // TODO: use scaled visibility
  jl_init();
  eval_string("using bamg");
  eval_string("using LinearAlgebra");
  eval_string("using BALUtils");
  eval_string("ba = readbal(\"" + bafile + "\")");
  eval_string("mg = create_multigrid(ba, Options())");
  auto mg = jl_get_global(jl_main_module, jl_symbol("mg"));
  auto update = get_function("update!", "bamg");

  PetscInt n;
  PetscBool done;
  const PetscInt *ia,*ja;
  MatGetRowIJ(A, 0, PETSC_FALSE, PETSC_FALSE, &n, &ia, &ja, &done);
  if (!done) {
    throw("PETSc could not give CSR matrix indices");
  }
  PetscScalar *data;
  MatSeqAIJGetArray(A, &data);

  int nnz = ia[n];
  jl_value_t* args[8] = { mg
                        , wrap_array(ia, n+1)
                        , wrap_array(ja, nnz)
                        , wrap_array(data, nnz)
                        , wrap_array(scale, n)
                        , wrap_array(solution, n)
                        , jl_box_voidpointer(NULL)
                        , jl_box_int64(nnz)
                        };
  jl_call(update, args, 8);
  check_error();

  MatSeqAIJRestoreArray(A, &data);
  MatRestoreRowIJ(A, 0, PETSC_TRUE, PETSC_FALSE, &n, &ia, &ja, &done);


  PCSetType(pc, PCSHELL);
  // PCSetUp(pc);
  PCShellSetName(pc, "bamg");
  PCShellSetApply(pc, bamg_apply);
}
