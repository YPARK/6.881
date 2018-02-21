# Helper functions for 6.881 lab

## Installation using `devtools`

First install packages needed, `devtools`, `Rcpp`, `RcppProgress`

```
install.packages('devtools')
install.packages('Rcpp')
install.packages('RcppProgress')
```

Then install this package
```
> library(devtools)
> devtools::install_github('YPARK/6.881/util')
```

## Installation from local repository (a simple version)

Just do this in `R`
```
> install.packages('util6881_1.0.0.tar.gz')
```

## Installation (optional)

Prerequisite: `Rcpp`, `RcppEigen`, `RcppProgress` packages

Make sure your R development environment support `C++14` by including
`-std=c++14` to `CFLAGS` and `CXXFLAGS` in `~/.R/Makevars` file.
For instance,
```
CXX = g++-6
CXXFLAGS = -O3 -std=c++14
CFLAGS = -O3 -std=c++14
```

Build package locally.
```
$ R CMD build .
```

You will have `'util6881_x.x.x.tar.gz` gzipped file in current directory
(`x.x.x` can be any version of the package).  Install package within
R:

```
> install.packages('util6881_x.x.x.tar.gz')
```

# Bug reports

Yongjin Park `ypp@csail.mit.edu`

