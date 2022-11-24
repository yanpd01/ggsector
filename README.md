<!-- README.md is generated from README.Rmd. Please edit that file -->

# ggsector: Easily draw sectors with grid and ggplot2

<!-- badges: start -->

[![Project Status: Active - The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
<!-- [![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/ggsector)](https://cran.r-project.org/package=ggsector) -->
[![](https://img.shields.io/badge/devel%20version-1.5.3-green.svg)](https://github.com/yanpd01/ggsector)

![](https://img.shields.io/badge/platforms-all-green.svg)
![](https://img.shields.io/badge/Windows-passing-green.svg)
![](https://img.shields.io/badge/Linux-passing-green.svg)
<!-- badges: end -->

The “ggsector” is an R package that can use “grid” or “ggplot2” to
easily draw sectors. At the same time, “ggsector” provides some cases of
how to call `grid.sector` in “ComplexHeatmap”, and provides a function
to interface with “Seurat” to facilitate the display of cell markers.

## :writing_hand: Authors

Pengdong Yan

## :arrow_double_down: Installation

Get the released version from Cran:

``` r
## The new version of the R package is still under review.
## If you want to try it out, you can use the following github method to install it
install.packages("ggsector")
```

Or the development version from github or gitee:

``` r
## install.packages("remotes")
# from github
remotes::install_github("yanpd01/ggsector")
# from gitee
remotes::install_git("https://gitee.com/yanpd01/ggsector")
```

## :books: Usage

For the usage of this R package, please type `vignette("ggsector")`
after installation to view it. There are detailed case descriptions in
it.

## :sparkling_heart: Acknowledgments

The code of this R package refers to
[jjplot](https://github.com/junjunlab/jjPlot) of
[JunJunLao](https://github.com/junjunlab) and
[ggplot2](https://github.com/tidyverse/ggplot2) of
[Hadley](https://github.com/hadley).

The Description, vignette, and readme of this R package refer to
[clusterProfiler](https://github.com/YuLab-SMU/clusterProfiler),
[ggfun](https://github.com/YuLab-SMU/ggfun), and
[treeio](https://github.com/YuLab-SMU/treeio) of [Guangchuang
YU](https://github.com/YuLab-SMU/).

Here, I would like to express my highest respect to thank the big guys
for their open source spirit.
