---
editor_options: 
  markdown: 
    wrap: 72
---

# BatChef

<img src="https://drive.google.com/uc?export=view&amp;id=1_XyW_NBxIw2el05dzVzXCapwr2ktgL4k" align="left" width="200"/>

`BatChef` is an R package that provides a common interface for batch
effect correction methods in single-cell RNA-seq data.

It allows users to compute various correction methods through a generic
function by specifying the desired approach. The package contains
multiple methods implemented in both R and Python, each based on
different mathematical approaches.

Also, it provides metrics to evaluate the performance of the correction
methods, such as the Wasserstein distance, Local Inverse Simpson's Index
(LISI), Average Silhouette Width (ASW), and Adjusted Rand Index (ARI).

## Installation

``` r
# install.packages("devtools")
devtools::install_github("zuinelena3/BatChef")
```
