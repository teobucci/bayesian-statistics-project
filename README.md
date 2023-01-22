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

1. Install gcc which includes gfortran with
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
packages_list <- c("tidyverse","mvtnorm","salso","logr","gmp","mcclust","igraph","ggraph","tidygraph", "uuid", "dittodb")
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
```

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
