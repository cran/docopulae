# docopulae 0.4

* new function `wDsensitivity`
* new function `wDefficiency`
* inconsistencies due to odd implementation of `duplicated` fixed
* new function `rowsduplicated`

# docopulae 0.3.5

* `DerivLogf`, `Deriv2Logf` and `update.param` print their progress
* new function `grow.grid`
* `update.param` handles interrupts by user and returns partially updated model object
* minor bug in C code fixed

# docopulae 0.3.4

* `buildf` takes an additional argument `continuous` to build joint probability mass functions

# docopulae 0.3.3

* `expr2f` removed in favor of package `Deriv`
  * `buildf` returns a function in every use case
  * arguments changed for `DerivLogf` and `Deriv2Logf`
  * `buildf`, `DerivLogf` and `Deriv2Logf` use package `Deriv` to simplify and cache
* `nint_transform` completely rewritten
  * takes a list of transformations
  * transforms interval and function dimensions
  * handles function dimensions correctly
  * implements transformation of function dimensions to interval dimensions
  * builtin transformations replaced by `nint_tanTransform`
* refocus on D_A-optimality, arguments changed for `Dsensitivity` and `Defficiency`
* argument `names` renamed for `buildf`, `DerivLogf`, `Deriv2Logf`, `fisherI`, `Dsensitivity`, `Defficiency`
* `FedorovWynn` renamed to `Wynn`
* defaults for `method.args` for `numDeriv2Logf` fixed
* handling of `nint_funcDim` in `nint_integrateNFunc` fixed
* performance increased for `Dsensitivity` and `Defficiency`
* `rowmatch` uses C-code and requires matrices of doubles
* `plot.desigh` recursively finds base design
* `plot.desigh` won't draw second axis anymore if `sensArgs$axes == F`
* typo fixed and default sensitivity label redefined for `plot.desigh`

# docopulae 0.3.2

* `numDerivLogf` and `numDeriv2Logf` optionally take a log likelihood function
* sensitivity function redefined to a positive function, concerns
  * `Dsensitivity`
  * argument `tol` to `FedorovWynn`
  * new argument `sensTol` to `plot.desigh`

# docopulae 0.3.1

* `buildf` is more general

# docopulae 0.3.0

* class `desigh` redefined, constructor `design` added
* major change in workflow
  * sensitivity function is 'liberated'
  * `Dsensitivity` builds a sensitivity function for D- and Ds-optimality
  * `FedorovWynn` takes a sensitivity function
  * related functions and their arguments are adjusted accordingly
* component `dDim` removed from model definition
* `update.param` also takes designs therefore replacing `update.desigh` and `update_reference`
* new defaults for arguments `dsNames` and `names`

# docopulae 0.2.0

* initial CRAN submission
* start package versioning
* all previous commits are considered irrelevant
