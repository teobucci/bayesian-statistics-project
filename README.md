# Stochastic Block Model Prior with Ordering Constraints for Gaussian Graphical Models

## How to clone the repository

```
git clone https://github.com/teobucci/bayesian-statistics-project
git submodule update --init
git submodule update --recursive
```

## How to update the submodule

```
git pull origin master
git submodule foreach git pull
git commit -a -m 'Update submodule with fixes'
git push origin master
```

## How to compile the PDF files

To compile the presentations, run the following in the root of the repo

```
make prese1
make prese2
```

## How to format R code

Using [this guide](https://bookdown.org/dli/rguide/r-style-guide.html).

In R you can use the commands: `Code` > `Reformat Code` to format the selected chunk of code.

