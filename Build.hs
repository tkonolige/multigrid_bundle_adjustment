import Development.Shake
import Development.Shake.Command
import Development.Shake.FilePath
import Development.Shake.Util
import System.Directory
import System.Info.Extra

bafiles = ["S-solution.petsc", "S-rhs.petsc", "S.petsc", "dscale.petsc", "solution.petsc"]
balfile dir = (dropFileName dir) </> "problem_noised.bbal"

baloutput f = [ "-fieldsplit_camera_ksp_view_solution"
              , "binary:" <> f </> "S-solution.petsc"
              , "-fieldsplit_camera_ksp_view_rhs"
              , "binary:" <> f </> "S-rhs.petsc"
              , "-fieldsplit_camera_ksp_view_pmat"
              , "binary:" <> f </> "S.petsc"
              , "-ksp_view_diagonal_scale"
              , "binary:" <> f </> "dscale.petsc"
              , "-tao_view_solution"
              , "binary:" <> f </> "solution.petsc"
              ]

bamgProject = do
  getDirectoryFiles "" ["bamg//*.jl", "bamg/Project.toml"] >>= need
  need ["bamg/Manifest.toml"]

main :: IO ()
main = shakeArgs shakeOptions{shakeFiles="_shake", shakeChange=ChangeModtimeAndDigest} $ do
    want []

    let libpetsc = if isMac then "petsc/arch-darwin-opt-deb/lib/libpetsc.dylib"
                   else "petsc/arch-linux-c-opt/lib/libpetsc.so"
    let libceres = "ceres-solver/cmake-build/lib/libceres.a"

    libpetsc %> \out -> do
      let arch = takeDirectory1 $ dropDirectory1 out
      srcs <- getDirectoryFiles "" ["petsc/include//*.h", "petsc/src//*.h", "petsc/src//*.c", "petsc/src//makefile"]
      need srcs
      cwd <- liftIO $ getCurrentDirectory
      cmd_ "make" ["PETSC_DIR=" <> cwd </> "petsc", "PETSC_ARCH=" <> arch] "-C petsc"

    libceres %> \out -> do
      srcs <- getDirectoryFiles "" ["ceres-solver/include//*.h", "ceres-solver/internal//*.h", "ceres-solver/internal//*.cc"]
      need srcs
      bamgProject
      cmd_ "cmake --build" [takeDirectory $ takeDirectory out] "-- -j 4"

    "ba-problems/*.problem/problem_noised.bbal" %> \out -> do
      let infile = (dropFileName out) </> "problem.bbal"
      srcs <- getDirectoryFiles "" ["city2ba/src//*.rs", "city2ba/Cargo.toml"]
      need srcs
      need [infile]
      lines <- readFileLines $ (dropFileName out) </> "noise.txt"
      cmd_ (Cwd "city2ba") "cargo run --release --bin city2ba noise" (map (\x -> ".." </> x) [infile, out]) lines

    (map (\f -> "ba-problems/*.problem" </> f) bafiles) &%> \outfiles -> do
      let bin = "ba-tao/build/bin/CeresBundleAdjustment"
      let balproblem = balfile $ head outfiles
      let outdir = dropFileName $ head outfiles
      let options = "ba-tao/default_options.txt"
      need [bin, balproblem, options]
      cmd_ bin "-bal" [balproblem] "-options_file" [options] (baloutput outdir)

    "ba-tao/build/bin/CeresBundleAdjustment" %> \out -> do
      srcs <- getDirectoryFiles "" ["ba-tao/cmake/*.cmake", "ba-tao/src/*hpp", "ba-tao/src/*.h", "ba-tao/src/*.cpp", "ba-tao/CMakeLists.txt"]
      need srcs
      need [libpetsc, libceres]
      cmd_ "cmake --build" [takeDirectory $ takeDirectory out] "-- -j 4"

    "bamg/Manifest.toml" %> \out -> do
      need ["bamg/Project.toml"]
      let command = "using Pkg; Pkg.add(PackageSpec(path=\"./BALUtils\")); Pkg.add(PackageSpec(path=\"./StaticArrays.jl\")); Pkg.add(PackageSpec(path=\"./IterativeSolvers.jl\")); Pkg.add(PackageSpec(path=\"./Arpack.jl\")); Pkg.instantiate()"
      cmd_ "julia --project=bamg -e " [command]

    "ba-problems/*.problem/benchmark_*_*.csv" %> \out -> do
      let balproblem = (dropFileName out) </> "problem_noised.bbal"
      bamgProject
      let Just [_, nlopts, lopts] = filePattern "ba-problems/*.problem/benchmark_*_*.csv" out
      llines <- readFileLines $ "options/linear_solver" </> lopts <.> "txt"
      nllines <- readFileLines $ "options/nonlinear_solver" </> nlopts <.> "txt"
      need [balproblem]
      cmd_ "julia -O3 --project=bamg bamg/bench/bench.jl" [balproblem] "--csv" [out] nllines llines

    "ba-problems/*.problem/nash_log_*_*.csv" %> \out -> do
      let balproblem = (dropFileName out) </> "problem_noised.bbal"
      bamgProject
      let Just [_, nlopts, lopts] = filePattern "ba-problems/*.problem/nash_log_*_*.csv" out
      llines <- readFileLines $ "options/linear_solver" </> lopts <.> "txt"
      nllines <- readFileLines $ "options/nonlinear_solver" </> nlopts <.> "txt"
      need [balproblem]
      cmd_ "julia -O3 --project=bamg bamg/bench/bench.jl" [balproblem] "--nash-log" [out] nllines llines

    "ba_eig/build/ba_eig" %> \out -> do
      need [libpetsc]
      need ["ba_eig/CMakeLists.txt", "ba_eig/ba_eig.cpp"]
      cmd_ "cmake --build ba_eig/build"

    "ba-problems/*.problem/eig_*_*.h5" %> \out -> do
      need ["ba_eig/build/ba_eig"]
      let Just [problem, index, slargest] = filePattern "ba-problems/*.problem/eig_*_*.h5" out
      let largest = case slargest of
            "largest" -> True
            "smallest" -> False
      let problem = dropDirectory1 $ dropFileName out
      let matrix = (dropFileName out) </> "default_dump.h5"
      need [matrix]
      let options = "ba_eig/options_" <> slargest <.> "txt"
      need [options]
      let bafile = (dropFileName out) </> "problem_noised.bbal"
      need [bafile]
      bamgProject
      let smopts = if not largest then ["-bafile", bafile] else []
      cmd_ "ba_eig/build/ba_eig -hdf5" [matrix] "-out" out "-index" index "-options_file" [options] smopts "-bdscale"

    "ba-problems/*.problem/condition_number_*.csv" %> \out -> do
      let dir = dropFileName out
      let Just [problem, index] = filePattern "ba-problems/*.problem/condition_number_*.csv" out
      let largest = dir </> "eig_" <> index <> "_largest.h5"
      let smallest = dir </> "eig_" <> index <> "_smallest.h5"
      need [largest, smallest]
      cmd_ "julia cond.jl" [smallest, largest, out]

    "ba-problems/*.problem/gt_condition_number.csv" %> \out -> do
      let dir = dropFileName out
      let Just [problem] = filePattern "ba-problems/*.problem/gt_condition_number.csv" out
      let largest = dir </> "gt_eig_largest.h5"
      let smallest = dir </> "gt_eig_smallest.h5"
      need [largest, smallest]
      cmd_ "julia cond.jl" [smallest, largest, out]

    ["ba-problems/*.problem/bamg_*_*_P.petsc", "ba-problems/*.problem/bamg_*_*_S.petsc"] &%> \[out, strength] -> do
      bamgProject
      let problem = dropFileName out
      let Just [_, index, opts] = filePattern "ba-problems/*.problem/bamg_*_*_P.petsc" out
      lines <- readFileLines $ "options/linear_solver" </> opts <.> "txt"
      need $ map (problem </>) bafiles
      cmd_ "julia --project=bamg bamg/scripts/build_prolongation.jl" [problem, index, out, strength] lines

    "ba-problems/*.problem/fine_grid_*_*.petsc" %> \out -> do
      let Just [_, index, opts] = filePattern "ba-problems/*.problem/fine_grid_*_*.petsc" out
      let matrix = replaceFileName out $ "bamg_" <> index <> "_" <> opts <> "_P.petsc"
      need ["multigrid_analysis/P_to_S.jl", matrix]
      cmd_ "julia multigrid_analysis/P_to_S.jl" [matrix, out]

    "multigrid_analysis/build/cr_analysis" %> \out -> do
      srcs <- getDirectoryFiles "" ["multigrid_analysis/*.cpp", "multigrid_analysis/cmake/*.cmake", "multigrid_analysis/CMakeLists.txt"]
      need srcs
      cmd_ "cmake --build multigrid_analysis/build"

    "ba-problems/*.problem/cr_eigs_*_*.h5" %> \out -> do
      let Just [_, index, opts] = filePattern "ba-problems/*.problem/cr_eigs_*_*.h5" out
      let a = replaceFileName out "S.petsc"
      let f = replaceFileName out $ "fine_grid_" <> index <> "_" <> opts <> ".petsc"
      let options = "multigrid_analysis/default_options.txt"
      need ["multigrid_analysis/build/cr_analysis", a, f, options]
      cmd_ "multigrid_analysis/build/cr_analysis -A" [a] "-S" [f] "-index" [index] "-out" [out] "-options_file" [options]

    "ba-problems/*.problem/convergence_test_study_*.csv" %> \out -> do
      let balproblem = (dropFileName out) </> "problem_noised.bbal"
      bamgProject
      let Just [_, opts] = filePattern "ba-problems/*.problem/convergence_test_study_*.csv" out
      lines <- readFileLines $ "options/linear_solver" </> opts <.> "txt"
      need [balproblem]
      withTempFile $ \feisen -> do
        withTempFile $ \fnash -> do
          cmd_ "julia -O3 --project=bamg bamg/bench/bench.jl" [balproblem] "--csv" [feisen] lines "--convergence eisenstat"
          cmd_ "julia -O3 --project=bamg bamg/bench/bench.jl" [balproblem] "--csv" [fnash] lines "--convergence nash --tol 0.01"
          cmd_ Shell "csvstack" [feisen, fnash] ">" [out]

    "ba-problems/*.problem/ceres_benchmark_*.csv" %> \out -> do
      let balproblem = (dropFileName out) </> "problem_noised.bbal"
      let bin = "ba-tao/build/bin/CeresBundleAdjustment"
      let Just [_, opts] = filePattern "ba-problems/*.problem/ceres_benchmark_*.csv" out
      lines <- readFileLines $ "options/ceres" </> opts <.> "txt"
      bamgProject
      need [bin, balproblem]
      cmd_ [bin] "--bal" [balproblem] lines "--csv" [out]

    "ba-problems/*.problem/gt_benchmark_*.csv" %> \out -> do
      let balproblem = (dropFileName out) </> "problem.bbal"
      let bin = "ba-tao/build/bin/CeresBundleAdjustment"
      let Just [_, opts] = filePattern "ba-problems/*.problem/gt_benchmark_*.csv" out
      lines <- readFileLines $ "options/ceres" </> opts <.> "txt"
      bamgProject
      need [bin, balproblem]
      cmd_ [bin] "--bal" [balproblem] lines "--csv" [out] "--no-damping --num_iter 1 --abstol -1 --random_rhs --rtol 1e-3"

    "ba-problems/*.problem/multigrid_benchmark_*.csv" %> \out -> do
      let balproblem = (dropFileName out) </> "problem_noised.bbal"
      let bin = "ba-tao/build/bin/CeresBundleAdjustment"
      let Just [_, opts] = filePattern "ba-problems/*.problem/multigrid_benchmark_*.csv" out
      lines <- readFileLines "options/ceres/multigrid_robust01.txt"
      let optsfile = "options/ceres/multigrid" </> opts <.> "txt"
      bamgProject
      need [bin, balproblem, optsfile]
      cmd_ [bin] "--bal" [balproblem] lines "--csv" [out] "--options_file" [optsfile]

    "ba-problems/*.problem/*.ply" %> \out -> do
      let balproblem = (dropExtension out) <.> ".bbal"
      srcs <- getDirectoryFiles "" ["city2ba/src//*.rs", "city2ba/Cargo.toml"]
      need srcs
      need [balproblem]
      cmd_ (Cwd "city2ba") "cargo run --release --bin city2ba ply" (map (\x -> ".." </> x) [balproblem, out])

    "ba-problems/*.problem/*_dump.h5" %> \out -> do
      let balproblem = (dropFileName out) </> "problem_noised.bbal"
      let bin = "ba-tao/build/bin/CeresBundleAdjustment"
      lines <- readFileLines "options/ceres/multigrid.txt"
      let Just [_, opts] = filePattern "ba-problems/*.problem/*_dump.h5" out
      let optsfile = "options/ceres/multigrid" </> opts <.> "txt"
      need [bin, balproblem, optsfile]
      bamgProject
      cmd_ [bin] "--bal" [balproblem] lines "--options_file" [optsfile] "--dump_file" [out] "--robust --num_iter 10"

    "ba-problems/*.problem/dump_gt.h5" %> \out -> do
      let balproblem = (dropFileName out) </> "problem.bbal"
      let bin = "ba-tao/build/bin/CeresBundleAdjustment"
      lines <- readFileLines "options/ceres/multigrid.txt"
      let optsfile = "options/ceres/multigrid/default.txt"
      need [bin, balproblem, optsfile]
      bamgProject
      cmd_ [bin] "--bal" [balproblem] lines "--options_file" [optsfile] "--dump_file" [out] "--num_iter 1 --no_damping"

    "ba-problems/*.problem/gt_eig_*.h5" %> \out -> do
      need ["ba_eig/build/ba_eig"]
      let Just [problem, slargest] = filePattern "ba-problems/*.problem/gt_eig_*.h5" out
      let largest = case slargest of
            "largest" -> True
            "smallest" -> False
      let problem = dropDirectory1 $ dropFileName out
      let matrix = (dropFileName out) </> "dump_gt.h5"
      need [matrix]
      let options = "ba_eig/options_" <> slargest <.> "txt"
      need [options]
      let bafile = (dropFileName out) </> "problem.bbal"
      need [bafile]
      bamgProject
      let smopts = if not largest then ["-bafile", bafile] else []
      cmd_ "ba_eig/build/ba_eig -hdf5" [matrix] "-out" out "-index 1" "-options_file" [options] smopts "-nullspace"

    "ba-problems/*.problem/diameter.csv" %> \out -> do
      let balproblem = (dropFileName out) </> "problem_noised.bbal"
      let bin = "diameter.jl"
      need [bin, balproblem]
      cmd_ "julia" [bin, balproblem, out]

    "ba-problems/bal-*/ceres_benchmark_*.csv" %> \out -> do
      let balproblem = (dropFileName out) </> "problem.bal"
      let bin = "ba-tao/build/bin/CeresBundleAdjustment"
      let Just [_, opts] = filePattern "ba-problems/*/ceres_benchmark_*.csv" out
      lines <- readFileLines $ "options/ceres" </> opts <.> "txt"
      need [bin, balproblem]
      bamgProject
      cmd_ [bin] "--bal" [balproblem] lines "--csv" [out]

    "ba-problems/*_block_drift.problem/problem.bbal" %> \out -> do
      srcs <- getDirectoryFiles "" ["city2ba/src//*.rs", "city2ba/Cargo.toml"]
      let Just [blocks] = filePattern "ba-problems/*_block_drift.problem/problem.bbal" out
      need srcs
      cmd_ (Cwd "city2ba") "cargo run --release --bin city2ba synthetic" (map (\x -> ".." </> x) [out]) "--blocks" [blocks]
