# docopulae

A direct approach to optimal designs for copula models based on the Fisher information.
Provides flexible functions for building joint PDFs, evaluating the Fisher information and finding Ds-optimal designs.
It includes an extensible solution to summation and integration called `nint`, functions for transforming, plotting and comparing designs, as well as a set of tools for common low-level tasks.

## Goals

*docopulae* strives to provide functions which allow the user to
* define a wide variety of models, by the joint probability density function (PDF)
* evaluate the Fisher information, using a convenient interface
* find Ds-optimal designs
* visualize designs
* compare Ds-optimal designs

## Workflow

See [./misc/workflow.dot.pdf](./misc/workflow.dot.pdf) and [./misc/nint.dot.pdf](./misc/nint.dot.pdf). Vignettes explaining these steps in detail are in production.

## Quick Start

TODO. For the moment see and follow the example for function `param` on the corresponding help page. Have fun :)

## Install

* from CRAN (hopefully soon)
  * `install.packages('docopulae')`
* from GitHub with devtools
  * `devtools::install_github('arappold/docopulae')`
* from GitHub without devtools
  * download https://github.com/arappold/docopulae/archive/master.zip
  * extract archive
  * `install.packages('/path/to/docopulae-master', repos=NULL, type='source')`

## Bugs

If you are absolutely certain that you found a bug, please let me know by creating an issue at https://github.com/arappold/docopulae/issues. Explain how to reproduce the bug, best by attaching a small script, and I will investigate as soon as I got time to.

(just a) **Warning**: *docopulae* allows complex/complicated scripts. And even though we might think we know what it tells *R* what to do, we most often don't.

## Feature Requests

If you feel unhappy about certain aspects of *docopulae* and perhaps have an adequate solution, please create an issue and lets discuss about it.
