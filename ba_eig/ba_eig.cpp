// Spectral analysis of Bundle Adjustment problems
// Usage: ./ba_eig -H H.petsc -poses number_of_poses -pose_size pose_size -out outfile.h5 -st_type sinvert -st_shift 0 -eps_target 0
// A.petsc and b.petsc come from the Tao bundle adjuster.
// Output file structure:
// "eigenvalues" => array of eigenvalues
// "eigenvector$(n)" => nth eigenvector real values as pose_size x number_of_poses array
// "eigenvector$(n)i" => nth eigenvector imaginary values as pose_size x number_of_poses array

// #include "bundle_adjustment_problem.h"
#include <iostream>
#include <petscviewerhdf5.h>
#include <petscblaslapack.h>
#include <petsc.h>
#include <petscoptions.h>
#include <slepc.h>
#include <string>
#include <vector>
#include "bamg_preconditioner.hpp"

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

void write_to_file(Mat A, const std::string filename) {
  PetscViewer viewer;
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_WRITE,
                        &viewer);
  MatView(A, viewer);
  PetscViewerDestroy(&viewer);
}

Vec load_from_file_vec(const std::string filename, int index) {
  Vec A;
  PetscViewer viewer;
  PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_READ,
                        &viewer);
  VecCreate(PETSC_COMM_WORLD, &A);
  VecSetFromOptions(A);
  for(int i = 0; i <= index; i++) {
    VecLoad(A, viewer);
  }
  PetscViewerDestroy(&viewer);
  return A;
}

int main(int argc, char **args) {
  SlepcInitialize(&argc, &args, (char *)0, NULL);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  bool leader = rank == 0;

  char jacobian_path[PETSC_MAX_PATH_LEN];
  jacobian_path[0] = '\0';
  PetscBool use_jacobian = PETSC_FALSE;

  char hessian_path[PETSC_MAX_PATH_LEN];
  hessian_path[0] = '\0';
  PetscBool use_hessian = PETSC_FALSE;

  char hdf5_path[PETSC_MAX_PATH_LEN];
  hdf5_path[0] = '\0';
  PetscBool use_hdf5 = PETSC_FALSE;

  char outfile[PETSC_MAX_PATH_LEN];
  outfile[0] = '\0';
  PetscBool use_output = PETSC_FALSE;
  char rhs_path[PETSC_MAX_PATH_LEN];
  int load_index = 0;
  PetscBool use_nullspace = PETSC_FALSE;
  PetscBool use_bamg = PETSC_FALSE;

  char bafile[PETSC_MAX_PATH_LEN];
  bafile[0] = '\0';
  PetscBool use_bafile = PETSC_FALSE;
  PetscBool blockdiag_scale = PETSC_FALSE;
  PetscBool bamg_mat = PETSC_FALSE;

  PetscOptionsBegin(MPI_COMM_WORLD, NULL, NULL, NULL);
  PetscOptionsString("-J", "Jacobian matrix", "", jacobian_path, jacobian_path, PETSC_MAX_PATH_LEN, &use_jacobian);
  PetscOptionsString("-H", "hessian matrix", "", hessian_path, hessian_path, PETSC_MAX_PATH_LEN, &use_hessian);
  PetscOptionsString("-hdf5", "Hessian matrix from HDF5", "", hdf5_path, hdf5_path, PETSC_MAX_PATH_LEN, &use_hdf5);
  PetscOptionsString("-out", "output file", "", outfile, outfile, PETSC_MAX_PATH_LEN, &use_output);
  PetscOptionsInt("-index", "which # mat to load from the file", "", load_index, &load_index, NULL);
  PetscOptionsBool("-nullspace", "Project out nullspace", "", use_nullspace, &use_nullspace, NULL);
  PetscOptionsBool("-bamg", "Use bamg as the PC", "", use_bamg, &use_bamg, NULL);
  PetscOptionsString("-bafile", "bundle adjustment file", "", bafile, bafile, PETSC_MAX_PATH_LEN, &use_bafile);
  PetscOptionsBool("-bdscale", "Scale problem by block diagonal", "", blockdiag_scale, &blockdiag_scale, NULL);
  PetscOptionsBool("-bamg_mat", "Find preconditioned operator eigenvectors", "", bamg_mat, &bamg_mat, NULL);
  PetscOptionsEnd();

  if (use_bamg && !use_bafile) {
    std::cerr << "Must specify -bafile with -bamg" << std::endl;
    return 1;
  }

  Mat S;
  Vec b;
  Vec nullspace[7];
  Vec scale;
  Vec poses;
  if(use_jacobian) {
    S = load_from_file(jacobian_path, load_index);
  } else if (use_hessian) {
    S = load_from_file(hessian_path, load_index);
    MatSetOption(S, MAT_SYMMETRIC, PETSC_TRUE);
  } else if(use_hdf5) {
    PetscViewer viewer;
    PetscViewerHDF5Open(PETSC_COMM_WORLD, hdf5_path, FILE_MODE_READ, &viewer);
    MatCreate(PETSC_COMM_WORLD, &S);
    MatSetFromOptions(S);
    MatSetBlockSize(S, 9);
    std::string name = std::to_string(load_index) + "/A";
    PetscObjectSetName((PetscObject)S,name.c_str());
    MatLoad(S, viewer);
    MatSetOption(S, MAT_SPD, PETSC_TRUE);

    if(use_nullspace) {
      for(int i = 0; i < 7; i++) {
        VecCreate(PETSC_COMM_WORLD,&nullspace[i]);
        std::string name = std::to_string(load_index) + "/nullspace" + std::to_string(i+1);
        PetscObjectSetName((PetscObject)nullspace[i],name.c_str());
        VecLoad(nullspace[i], viewer);
      }
    }

    VecCreate(PETSC_COMM_WORLD,&scale);
    PetscObjectSetName((PetscObject)scale,(std::to_string(load_index) + "/scale").c_str());
    VecLoad(scale, viewer);

    VecCreate(PETSC_COMM_WORLD,&poses);
    PetscObjectSetName((PetscObject)poses,(std::to_string(load_index) + "/poses").c_str());
    VecLoad(poses, viewer);

    PetscViewerDestroy(&viewer);
  } else {
    std::cout << "-J, -H, or -hdf5 must be specified" << std::endl;
    std::exit(1);
  }

  if(blockdiag_scale) {
    PetscInt m, n, M, N, bs;
    MatGetLocalSize(S, &m, &n);
    MatGetSize(S, &M, &N);
    MatGetBlockSize(S, &bs);
    PetscScalar *diag;
    PetscMalloc1(m*bs, &diag);
    PetscInt *inds;
    PetscMalloc1(bs, &inds);
    PetscFPTrapPush(PETSC_FP_TRAP_OFF);
    PetscBLASInt bbs;
    PetscBLASIntCast(bs,&bbs);
    for(PetscInt i = 0; i < m/bs; i++) {
      for(int j = 0; j < bs; j++) inds[j] = i*bs+j;
      MatGetValues(S, bs, inds, bs, inds, diag+i*bs*bs);
      // zero the lower half of the block
      for(int j = 0; j < bs; j++) {
        for(int k = j+1; k < bs; k++) {
          diag[i*bs*bs+j*bs+k] = 0;
        }
      }
      PetscBLASInt bierr;
      PetscStackCallBLAS("LAPACKpotrf",LAPACKpotrf_("U",&bbs,diag+i*bs*bs,&bbs,&bierr));
      if (bierr) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_LIB,"Error in POTRF Lapack routine %d",(int)bierr);
      PetscStackCallBLAS("LAPACKtrtri",LAPACKtrtri_("U","N",&bbs,diag+i*bs*bs,&bbs,&bierr));
      if (bierr) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_LIB,"Error in TFTRI Lapack routine %d",(int)bierr);
    }
    PetscFPTrapPop();

    Mat BD;
    MatCreate(MPI_COMM_WORLD,&BD);
    MatType type;
    MatGetType(S, &type);
    MatSetType(BD, type);
    MatSetSizes(BD,m,n,M,N);
    MatSetBlockSize(BD,bs);
    PetscInt *dnnz;
    PetscInt rstart,rend,i,j;
    PetscMalloc1(m/bs,&dnnz);
    for (j = 0; j < m/bs; j++) dnnz[j] = 1;
    MatXAIJSetPreallocation(BD,bs,dnnz,NULL,NULL,NULL);
    PetscFree(dnnz);
    MatGetOwnershipRange(BD,&rstart,&rend);
    MatSetOption(BD,MAT_ROW_ORIENTED,PETSC_FALSE);
    for (i = rstart/bs; i < rend/bs; i++) {
      MatSetValuesBlocked(BD,1,&i,1,&i,&diag[(i-rstart/bs)*bs*bs],INSERT_VALUES);
    }
    MatAssemblyBegin(BD,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(BD,MAT_FINAL_ASSEMBLY);
    MatSetOption(BD,MAT_ROW_ORIENTED,PETSC_TRUE);
    write_to_file(S, "S.petsc");

    Mat S_;
    MatPtAP(S, BD, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &S_);
    MatDestroy(&S);
    S = S_;
  }

  PetscInt N, M;
  MatGetSize(S, &N, &M);
  PetscPrintf(PETSC_COMM_WORLD, "Loaded %d x %d matrix\n", N, M);

  PetscInt nconv;

  EPS eps;
  SVD svd;
  if(use_jacobian) {
    SVDCreate(PETSC_COMM_WORLD, &svd);
    SVDSetOperator(svd, S);
    SVDSetFromOptions(svd);

    SVDSolve(svd);

    SVDGetConverged(svd, &nconv);
  } else {
    // Set up eigen solver
    EPSCreate(PETSC_COMM_WORLD, &eps);
    EPSSetOperators(eps, S, NULL);
    EPSSetFromOptions(eps);
    EPSSetProblemType(eps, EPS_HEP);
    if(use_hdf5 && use_nullspace) {
      EPSSetDeflationSpace(eps,7,nullspace);
    }
    EPSSetUp(eps);

    if(use_hdf5 && use_bamg && use_bafile) {
      if(bamg_mat) {
        PC pc;
        const PetscScalar *poses_ary;
        const PetscScalar *scale_ary;
        VecGetArrayRead(poses, &poses_ary);
        VecGetArrayRead(scale, &scale_ary);
        PCCreate(PETSC_COMM_WORLD,&pc);
        bamg_create(pc, S, poses_ary, scale_ary, std::string(bafile));
        VecRestoreArrayRead(poses, &poses_ary);
        VecRestoreArrayRead(scale, &scale_ary);
        Mat Sc;
        Mat bamgmat;
        MatCreateShell(PETSC_COMM_WORLD, M, N, M, N, pc, &bamgmat);
        MatShellSetOperation(bamgmat,MATOP_MULT,(void(*)(void))bamg_apply);
        Mat mats[2];
        mats[1] = bamgmat;
        mats[0] = S;
        S = NULL;
        MatCreateComposite(PETSC_COMM_WORLD, 2, mats, &S);
      } else {
        ST st;
        KSP ksp;
        PC pc;
        EPSGetST(eps, &st);
        STGetKSP(st, &ksp);
        const PetscScalar *poses_ary;
        const PetscScalar *scale_ary;
        VecGetArrayRead(poses, &poses_ary);
        VecGetArrayRead(scale, &scale_ary);
        KSPGetPC(ksp, &pc);
        bamg_create(pc, S, poses_ary, scale_ary, std::string(bafile));
        VecRestoreArrayRead(poses, &poses_ary);
        VecRestoreArrayRead(scale, &scale_ary);
      }
    }

    EPSSolve(eps);

    EPSGetConverged(eps, &nconv);
  }

  // print eigenvalues
  std::vector<PetscScalar> eigs;
  if (leader) {
    std::cout << std::endl << "Eigenvalues" << std::endl;
    for (int i = 0; i < nconv; i++) {
      PetscScalar vr, vi;
      vi = 0;
      if(use_jacobian) {
        SVDGetSingularTriplet(svd, i, &vr, NULL, NULL);
      } else {
        EPSGetEigenvalue(eps, i, &vr, &vi);
      }
      std::cout << i << ": " << vr << " + " << vi << "i" << std::endl;
      eigs.push_back(vr);
    }
  }

  if(use_output) {
    // write eigenvalues
    PetscViewer viewer;
    PetscViewerHDF5Open(MPI_COMM_WORLD, outfile, FILE_MODE_WRITE,
        &viewer);
    Vec eigs_vec;
    VecCreateMPIWithArray(MPI_COMM_WORLD, 1, eigs.size(), nconv, eigs.data(),
        &eigs_vec);
    PetscObjectSetName((PetscObject)eigs_vec, "eigenvalues");
    VecView(eigs_vec, viewer);

    // write eigenvectors
    if(leader) std::cout << std::endl << "Writing to " << outfile << std::endl;
    for (int i = 0; i < nconv; i++) {
      Vec x;
      Vec x_i;
      if(use_jacobian) {
        MatCreateVecs(S, &x_i, &x);
        SVDGetSingularTriplet(svd, i, NULL, x, x_i);
      } else {
        MatCreateVecs(S, &x, NULL);
        MatCreateVecs(S, &x_i, NULL);
        EPSGetEigenvector(eps, i, x, x_i);
      }
      PetscObjectSetName((PetscObject)x,
          ("eigenvector" + std::to_string(i)).c_str());
      VecView(x, viewer);
      PetscObjectSetName((PetscObject)x_i,
          ("eigenvector" + std::to_string(i) + "i").c_str());
      VecView(x_i, viewer);
    }
  }
}
