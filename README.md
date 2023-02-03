<!-- omit from toc -->
# Stochastic Block Model Prior with Ordering Constraints for Gaussian Graphical Models

![GitHub last commit](https://img.shields.io/github/last-commit/teobucci/CNN-Plants-Classifier?logo=github)

This project was developed for the course of **Bayesian Statistics** for the MSc. in Mathematical Engineering at Politecnico di Milano, A.Y. 2022/2023.

<!-- omit from toc -->
# Table of contents

- [Installation](#installation)
  - [How to clone the repository](#how-to-clone-the-repository)
  - [How to build the `FGM` package](#how-to-build-the-fgm-package)
  - [How to install the packages](#how-to-install-the-packages)
  - [How to compile the PDF files](#how-to-compile-the-pdf-files)
- [Running the analysis](#running-the-analysis)
- [Authors](#authors)
- [Developer notes](#developer-notes)
  - [How to update the submodule](#how-to-update-the-submodule)
  - [How to format R code](#how-to-format-r-code)
  - [Some help for debugging on macOS](#some-help-for-debugging-on-macos)


# Installation

## How to clone the repository

```
git clone https://github.com/teobucci/bayesian-statistics-project
git submodule update --init
git submodule update --recursive
```

## How to build the `FGM` package

Open the `./FGM/FGM.Rproj` in RStudio and type:

- `Ctrl+Shift+B` on Windows
- `CMD+Shift+B` on macOS

On macOS on M1 chip you may get an error involving `gfortran`, in which case proceed as follows according to [this](https://stackoverflow.com/a/72997915/16222204):

1. Install `gcc` which includes `gfortran` with
   ```
   brew install gcc
   ```
2. Create a file `~/.R/Makevars` (if it does not exist yet). For example running with a terminal
    ```
    mkdir -p ~/.R
    touch ~/.R/Makevars
    ```
3. Add the following lines to `~/.R/Makevars`
    ```
    FC = /opt/homebrew/Cellar/gcc/11.3.0_2/bin/gfortran
    F77 = /opt/homebrew/Cellar/gcc/11.3.0_2/bin/gfortran
    FLIBS = -L/opt/homebrew/Cellar/gcc/11.3.0_2/lib/gcc/11
    ```
    This can be done by opening it in a normal text editor such as VSCode (`code ~/.R/Makevars`) or SublimeText (`subl ~/.R/Makevars`).
    
    Note that you might have to change gcc version `11.3.0_2` to whatever your gcc version is.

## How to install the packages

Install the required packages from CRAN

```
packages_list <- c("tidyverse","mvtnorm","salso","logr","gmp","mcclust","igraph","ggraph","tidygraph", "uuid", "dittodb", "latex2exp", "kableExtra", "doSNOW", "doParallel")
install.packages(packages_list)
```

and install the custom utilities by Alessandro Colombi and `mcclust.ext`

```
devtools::install_github("alessandrocolombi/ACutils")
devtools::install_github("sarawade/mcclust.ext")
```

## How to compile the PDF files

To compile the presentations, run the following in the root of the repo

```
make prese1
make prese2
make prese3
```

To compile the report, run

```
make report
```

To compile everything, run

```
make pdf
```

To remove temporary `LaTeX` files, run

```
make clean
```

To remove both temporary and pdf files, run

```
make distclean
```

# Running the analysis

The repository contains different files to perform the analysis

- `01_simulations_basic.Rmd` is a notebook containing a vanilla implementation for running a single simulation, meant to be used a playground for on-the-go configurations.
- The second block of files implements a gridsearch approach to run different simulations varying parameters to see how well the MCMC behaves and how robust it is:
  - `02a_simulation_grid_generation.R` generates the grid of required configurations.
  - `02b_simulation_grid_parallel_execution.R` runs and saves all the simulations from the grid by the previous file, which is sourced here. The execution is run in parallel and on a MacBook Pro M1 14" takes about 10 minutes.
  - `02c_simulation_grid_visualization_notebook.Rmd` reads all the simulations generated from the previous files and, when knitted, produces a PDF where for each section there is a simulation. At the beginning there is a comprehensive table with all the relevant indexes across the grid.
  - `02d_simulation_grid_visualization_export_files.R` is a script that, when sourced, reads all the simulations generated by `02b_simulation_grid_parallel_execution.R` and saves all the relevant plots and figures to file, useful for embedding in presentations and report.
  - `02e_simulation_grid_kl_comparison.Rmd` reads all the simulations generated from the previous files and, when knitted, produces a PDF comparing the evolution of the KL distance across iterations for different configurations. For a better understanding, it is advised to run the simulations without burn-in in this case. 
- `03_simulations_real_dataset.Rmd` is a notebook where the algorithm is run on a real dataset, which is meant to be stored in `dataset`. It is essentially a copy of `01_simulations_basic.Rmd` but without the knowledge of the true graph and partition.
- `04_execution_time_regression.R` is a script that implements a polynomial linear regression of the execution time against the number of nodes, taken from the simulation grid results.

# Authors

Supervisor: Alessandro Colombi ([@alessandrocolombi](https://github.com/alessandrocolombi))

- Teo Bucci ([@teobucci](https://github.com/teobucci))
- Filippo Cipriani ([@SmearyTundra](https://github.com/SmearyTundra))
- Filippo Pagella ([@effefpi2](https://github.com/effefpi2))
- Flavia Petruso ([@fl-hi1](https://github.com/fl-hi1))
- Andrea Puricelli ([@apuri99](https://github.com/apuri99))
- Giulio Venturini ([@Vinavil334](https://github.com/Vinavil334))

# Developer notes

## How to update the submodule

```
git pull origin master
git submodule foreach git pull
git commit -a -m 'Update submodule with fixes'
git push origin master
```

## How to format R code

Using [this guide](https://bookdown.org/dli/rguide/r-style-guide.html).

In R you can use the commands: `Code` > `Reformat Code` to format the selected chunk of code.

## Some help for debugging on macOS

Open an R console in this way
```
R -d lldb
```
then
```
run
```
then
```
source("src/main.R")
```

More [stuff here](https://blog.davisvaughan.com/posts/2019-04-05-debug-r-package-with-cpp/).
