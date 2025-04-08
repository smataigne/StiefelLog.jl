# StiefelLog.jl
This repository is associated with the publication ["Mataigne S., Zimmermann R., Miolane N., An Efficient Algorithm for the Riemannian Logarithm on the Stiefel Manifold for a Family of Riemannian Metrics, SIAM Journal on Matrix Analysis and Applications, 46(2), 879-905, 2025"](https://epubs.siam.org/doi/10.1137/24M1647801)

## Requirements
To use this repository, the user must have a Julia installation and have installed the following packages for computations: `LinearAlgebra`, `SkewLinearAlgebra`, `MatrixEquations`. And also the following packages to make the plots: `Plots`, `Colors`, `LaTeXStrings`, `XSLX`, `Distributions`, `BenchmarkTools`. These packages are easily obtained from the package installation environment as follows. In Julia REPL, press `]` to access the installation environment and for each package, do
```julia
(@v1.6) pkg> add Name_of_Package
```


## Use

This repository contains three principal folders.
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

## Bibtex
If you use the content of this repository, please cite
```
@article{doi:10.1137/24M1647801,
author = {Mataigne, Simon and Zimmermann, Ralf and Miolane, Nina},
title = {An Efficient Algorithm for the Riemannian Logarithm on the Stiefel Manifold for a Family of Riemannian Metrics},
journal = {SIAM Journal on Matrix Analysis and Applications},
volume = {46},
number = {2},
pages = {879-905},
year = {2025},
doi = {10.1137/24M1647801},
URL = {https://doi.org/10.1137/24M1647801},
eprint = {https://doi.org/10.1137/24M1647801}
}
```
      
