# TMB_contrib_R

This repository contains user contributed R packages.
Install with ```devtools:::install_github(...)```.

### TMBhelper -- Dumping ground for useful contributed R functions
Includes
* AICTMB -- calculate AIC based on model output
* Check_Identifiable -- automatically check for non-identiable fixed effects

Install using
```R
devtools::install_github("kaskr/TMB_contrib_R/TMBhelper")
```

### TMBdebug -- Debugging tool for Windows users
* Designed to prevent crash in R terminal for windows users when the program overruns a vector or matrix
* Install using
```R
devtools::install_github("kaskr/TMB_contrib_R/TMBdebug")
```

### TMBphase -- Wrapper for fitting models in TMB with phases (a la ADMB)
* Designed to fit TMB models with phases, to allow sequential estimation of model parameters
* Install using
```R
devtools::install_github("kaskr/TMB_contrib_R/TMBphase")
```

### [TMBtools](https://github.com/mlysy/TMBtools) -- Tools for developing R packages interfacing with TMB

Provides helper functions for creating packages which contain TMB source code such that:

* Size of compiled code for multiple TMB models is minimized

* The package can contain non TMB-related C++ source code, e.g., using [Rcpp](http://www.rcpp.org/) 

* TMB compile chain has been cross-platform tested to pass `R CMD --as-cran check`

Install using
```R
devtools::install_github("mlysy/TMBtools")
```

Please see quickstart instructions [here](https://github.com/mlysy/TMBtools) and package vignette [here](http://htmlpreview.github.io/?https://github.com/mlysy/TMBtools/master/doc/TMBtools.html).
