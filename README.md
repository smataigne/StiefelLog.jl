# StiefelLog.jl
This repository is associated with the publication "Mataigne S., Zimmermann R., Miolane N., An efficient algorithm to compute the Riemannian logarithm on the Stiefel manifold for a family of Riemannian metrics"

To use this repository, the user must have a Julia installation and have installed the following packages for computations: `LinearAlgebra`, `SkewLinearAlgebra`, `MatrixEquations`. And also the following packages to make the plots: `Plots`, `Colors`, `LaTeXStrings`, `XSLX`, `Distributions`, `BenchmarkTools`.

It is made of three principal folders.
* The folder `src` contains the files:
  * `Manifold.jl` that must be run in local before any use of the other functions. If you want to use the functions implemented in `src`, simply add `ìmport .Manifold as mfd` to the top of your personal file.
  * `log.jl` contains an implementation of the routines introduced in the paper, notably Algorithm 3.1, 4.1 and 4.2 in their different versions (pure forward, pseudo-backward and accelerated forward).
  * The other files contain supplementary material, not directly related to the paper.
    * `stiefel.jl` contains the `StiefelVector` type utilized in `log.jl` and implements basic operations on the Stiefel manifold.
    * `orthonormal.jl` contains the `Orthonormal` type to deal with non-square orthogonal matrices and perform basic operations on these matrices.
    * and other files not related to the publication.
* The folder `make_figures_and_tables` provides the code to reperform every numerical experiment that is presented in the paper, being Figures 2 to 5 and Tables 1 to 3. Run `Manifold.jl` in local before running any of these files.
  * The folder `benchmark` contains `benchmark.jl`. Choose, $n, p$ and $\delta$ and run the file to output the table of performances.
  * The folder `logistic_regression` contains `Logistic_data.jl` to create the xlsx data files in the folder `logistic_data`. In `Logistic_data.jl`, you can choose $n$, $p$ and the number `N`of samples.  `Logistic_regression.jl` creates the figures in the folder `logistic_figures`.
  * The folder `make_plots` contains the files to create the figures 2 to 5 in the folder `figures`.
* The folder `bin`can be ignored.
      
