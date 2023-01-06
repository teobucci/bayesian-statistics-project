![GitHub last commit](https://img.shields.io/github/last-commit/teobucci/CNN-Plants-Classifier?logo=github)

# Stochastic Block Model Prior with Ordering Constraints for Gaussian Graphical Models

This project was developed for the course of **Bayesian Statistics** for the MSc. in Mathematical Engineering at Politecnico di Milano, A.Y. 2022/2023.

# Instructions

## How to clone the repository

```
git clone https://github.com/teobucci/bayesian-statistics-project
git submodule update --init
git submodule update --recursive
```

## How to build the `FGM` package

Open the `FGM.rproj` and type:

- `Ctrl+Shift+B` on Windows
- `CMD+Shift+B` on macOS

## How to install the packages

devtools::install_github("alessandrocolombi/ACutils")
packages_list <- c("tidyverse","mvtnorm","salso")
install.packages(packages_list)

## How to update the submodule

```
git pull origin master
git submodule foreach git pull
git commit -a -m 'Update submodule with fixes'
git push origin master
```

# Output

## How to compile the PDF files

To compile the presentations, run the following in the root of the repo

```
make prese1
make prese2
```

## How to format R code

Using [this guide](https://bookdown.org/dli/rguide/r-style-guide.html).

In R you can use the commands: `Code` > `Reformat Code` to format the selected chunk of code.

## Authors

Supervisor: Alessandro Colombi ([@alessandrocolombi](https://github.com/alessandrocolombi))

- Teo Bucci ([@teobucci](https://github.com/teobucci))
- Filippo Cipriani ([@SmearyTundra](https://github.com/SmearyTundra))
- Filippo Pagella ([@effefpi2](https://github.com/effefpi2))
- Flavia Petruso ([@fl-hi1](https://github.com/fl-hi1))
- Andrea Puricelli ([@apuri99](https://github.com/apuri99))
- Giulio Venturini ([@Vinavil334](https://github.com/Vinavil334))
