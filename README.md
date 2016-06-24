# TMB_contrib_R

This repository contains user contributed R packages.
Install with ```devtools:::install_github(...)```.

### TMBhelper -- Dumping ground for useful contributed R functions
Includes
* AICTMB -- calculate AIC based on model output

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
