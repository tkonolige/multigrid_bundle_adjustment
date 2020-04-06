// Find eigenvectors/values of SᵀA'S where A' is A normalized by columns and rows
// From Brannick, James J., and Robert D. Falgout. "Compatible relaxation and coarsening in algebraic multigrid." SIAM Journal on Scientific Computing 32.3 (2010): 1393-1416.
#include <petsc.h>
#include <slepc.h>
#include <iostream>
#include <vector>
#include "petscviewerhdf5.h"

using namespace std;

Mat load_from_file(const std::string filename, int index) {
  Mat A;
  PetscViewer viewer;
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_READ,
                        &viewer);
  MatCreate(PETSC_COMM_WORLD, &A);
  MatSetFromOptions(A);
  for(int i = 0; i <= index; i++) {
    MatLoad(A, viewer);
  }
  PetscViewerDestroy(&viewer);
  return A;
}

PetscErrorCode ksp_solve(Mat m, Vec x, Vec y) {
  KSP ctx;
  MatShellGetContext(m, &ctx);
  KSPSolve(ctx, x, y);
  return 0;
}

int main(int argv, char** argc) {
  char a_file[PETSC_MAX_PATH_LEN];
  a_file[0] = '\0';
  char s_file[PETSC_MAX_PATH_LEN];
  s_file[0] = '\0';
  char out_file[PETSC_MAX_PATH_LEN];
  out_file[0] = '\0';
  PetscErrorCode ierr;
  PetscBool use_composite = PETSC_FALSE;
  PetscInt index = 0;

  ierr = SlepcInitialize(&argv, &argc, (char*)0, NULL);CHKERRQ(ierr);

  ierr = PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "Compatible Relaxation Analysis", NULL);CHKERRQ(ierr);
  ierr = PetscOptionsString("-A", "A matrix", "", a_file, a_file, PETSC_MAX_PATH_LEN, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsString("-S", "S matrix", "", s_file, s_file, PETSC_MAX_PATH_LEN, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsString("-out", "output.h5", "", out_file, out_file, PETSC_MAX_PATH_LEN, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-use_composite", "Use composite matrix", NULL, use_composite, &use_composite, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsInt("-index", "Index of A to use", NULL, index, &index, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();CHKERRQ(ierr);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  bool leader = rank == 0;

  Mat S = load_from_file(s_file, 0);
  Mat A = load_from_file(a_file, index);
  Mat St;
  ierr = MatTranspose(S, MAT_INITIAL_MATRIX, &St);CHKERRQ(ierr);

  PetscInt M, N, m, n;
  ierr = MatGetSize(A, &M, &N);CHKERRQ(ierr);
  ierr = MatGetLocalSize(A, &m, &n);CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD, "A: %d x %d\n", M, N);

  Vec diag;
  ierr = MatCreateVecs(A, NULL, &diag);CHKERRQ(ierr);
  ierr = MatGetDiagonal(A, diag);CHKERRQ(ierr);
  // VecSqrtAbs(diag);

  Mat M_;
  ierr = MatCreateAIJ(MPI_COMM_WORLD, m, m, M, M, 1, NULL, 0, NULL, &M_);CHKERRQ(ierr);
  ierr = MatDiagonalSet(M_, diag, INSERT_VALUES);CHKERRQ(ierr);

  ierr = VecReciprocal(diag);CHKERRQ(ierr);

  Mat Minv;
  ierr = MatCreateAIJ(MPI_COMM_WORLD, m, n, M, N, 1, NULL, 0, NULL, &Minv);CHKERRQ(ierr);
  ierr = MatInvertBlockDiagonalMat(A, Minv);CHKERRQ(ierr);
  // ierr = MatDiagonalSet(Minv, diag, INSERT_VALUES);CHKERRQ(ierr);
  // MatDiagonalScale(A, diag, diag);



  ierr = MatGetSize(S, &M, &N);CHKERRQ(ierr);
  ierr = MatGetLocalSize(S, &m, &n);CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD, "S: %d x %d\n", M, N);

  // create (SᵀM⁻¹S)SᵀAS
  // TODO: why do we use SᵀM⁻¹S instead of (SᵀMS)⁻¹
  Mat C;
  if(use_composite) {
    // Mat mats[6] = {St, Minv, S, St, A, S};
    ierr = MatCreate(MPI_COMM_WORLD, &C);CHKERRQ(ierr);
    ierr = MatSetSizes(C, n, n, N, N);CHKERRQ(ierr);
    ierr = MatSetType(C, MATCOMPOSITE);CHKERRQ(ierr);
    ierr = MatCompositeSetType(C, MAT_COMPOSITE_MULTIPLICATIVE);CHKERRQ(ierr);
    ierr = MatCompositeAddMat(C, S);CHKERRQ(ierr);
    ierr = MatCompositeAddMat(C, A);CHKERRQ(ierr);

    // TODO: these two sometime cause issues
    ierr = MatCompositeAddMat(C, St);CHKERRQ(ierr);
    ierr = MatCompositeAddMat(C, S);CHKERRQ(ierr);

    ierr = MatCompositeAddMat(C, Minv);CHKERRQ(ierr);
    ierr = MatCompositeAddMat(C, St);CHKERRQ(ierr);

    Mat StMS;
    ierr = MatPtAP(M_, S, MAT_INITIAL_MATRIX, PETSC_DECIDE, &StMS);CHKERRQ(ierr);

    Mat StMS_inv;
    KSP ksp;
    ierr = KSPCreate(MPI_COMM_WORLD, &ksp);CHKERRQ(ierr);
    ierr = KSPSetOperators(ksp, StMS, StMS);CHKERRQ(ierr);
    ierr = MatCreateShell(MPI_COMM_WORLD, n, n, N, N, ksp, &StMS_inv);CHKERRQ(ierr);
    ierr = MatShellSetOperation(StMS_inv, MATOP_MULT, (void(*)())ksp_solve);CHKERRQ(ierr);

    // MatCompositeAddMat(C, StMS_inv);

    // MatCompositeAddMat(C, S);
    ierr = MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatCompositeMerge(C);CHKERRQ(ierr);
  } else {
    Mat StMinvS;
    ierr = MatPtAP(Minv, S, MAT_INITIAL_MATRIX, PETSC_DECIDE, &StMinvS);CHKERRQ(ierr);
    Mat StAS;
    ierr = MatPtAP(A, S, MAT_INITIAL_MATRIX, PETSC_DECIDE, &StAS);CHKERRQ(ierr);
    ierr = MatMatMult(StMinvS, StAS, MAT_INITIAL_MATRIX, PETSC_DECIDE, &C);CHKERRQ(ierr);

    // Mat MinvA;
    // ierr = MatMatMult(Minv, A, MAT_INITIAL_MATRIX, PETSC_DECIDE, &MinvA);CHKERRQ(ierr);
    // ierr = MatPtAP(MinvA, S, MAT_INITIAL_MATRIX, PETSC_DECIDE, &C);CHKERRQ(ierr);
  }
  if(leader) std::cout << "matrix created" << std::endl;
  ierr = PetscViewerPushFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_INFO);CHKERRQ(ierr);
  ierr = MatView(C, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscViewerPopFormat(PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  // set up eigensolver
  EPS eps;
  ierr = EPSCreate(PETSC_COMM_WORLD, &eps);CHKERRQ(ierr);
  ierr = EPSSetOperators(eps, C, NULL);CHKERRQ(ierr);
  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);

  // solve
  ierr = EPSSolve(eps);CHKERRQ(ierr);

  ierr = EPSReasonView(eps, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  EPSConvergedReason reason;
  ierr = EPSGetConvergedReason(eps, &reason);CHKERRQ(ierr);
  if(reason < 0) {
    SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "EPS failed to converge");
  }

  // converged pairs
  PetscInt nconv;
  ierr = EPSGetConverged(eps, &nconv);CHKERRQ(ierr);
  if(leader) std::cout << "# converged eigenvalues " << nconv << std::endl;

  std::vector<PetscScalar> eigs;
  for (int i = 0; i < nconv; i++) {
    PetscScalar vr, vi;
    ierr = EPSGetEigenvalue(eps, i, &vr, &vi);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "Eigenvalue %d: % 06.1e + % 06.1ei\n", i, vr, vi);CHKERRQ(ierr);
    if(leader)
      eigs.push_back(vr);
  }

  PetscViewer viewer;
  ierr = PetscViewerHDF5Open(MPI_COMM_WORLD, out_file, FILE_MODE_WRITE, &viewer);CHKERRQ(ierr);

  Vec eigs_vec;
  ierr = VecCreateMPIWithArray(MPI_COMM_WORLD, 1, eigs.size(), nconv, eigs.data(), &eigs_vec);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)eigs_vec, "eigenvalues");CHKERRQ(ierr);
  ierr = VecView(eigs_vec, viewer);CHKERRQ(ierr);

  // write out eigenvectors
  Vec x;
  Vec corrected;
  ierr = MatCreateVecs(C, &x, NULL);CHKERRQ(ierr);
  ierr = MatCreateVecs(S, NULL, &corrected);CHKERRQ(ierr);
  for (int i = 0; i < nconv; i++) {
    ierr = EPSGetEigenvector(eps, i, x, NULL);CHKERRQ(ierr);

    // Eigenvectors are in the fine space only; we want them in the full space.
    // S: fine space -> full space
    ierr = MatMult(S, x, corrected);CHKERRQ(ierr);

    ierr = PetscObjectSetName((PetscObject)corrected, ("eigenvector" + std::to_string(i)).c_str());CHKERRQ(ierr);
    ierr = VecView(corrected, viewer);CHKERRQ(ierr);

    ierr = PetscObjectSetName((PetscObject)x, ("eigenvector_uncorrected" + std::to_string(i)).c_str());CHKERRQ(ierr);
    ierr = VecView(x, viewer);CHKERRQ(ierr);
  }
}

