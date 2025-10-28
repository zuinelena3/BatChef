# BatChef

<img src="https://drive.google.com/uc?export=view&amp;id=1_XyW_NBxIw2el05dzVzXCapwr2ktgL4k" align="left" width="200"/>

`BatChef` is an R package that:

-   implements a variety of correction methods in R and Python. It allows users to compute the desired approach through a generic function;

-   provides metrics to evaluate the performance of the correction methods, such as the Wasserstein distance, Local Inverse Simpson's Index (LISI), Average Silhouette Width (ASW), and Adjusted Rand Index (ARI);

-   can be used as a guideline to identify the optimal batch effects correction method based on the data features.

## Installation

``` r
# install.packages("devtools")
devtools::install_github("zuinelena3/BatChef")
```

## Usage

See the vignette for more details.

``` r
sce <- simulate_data(n_genes = 2000, batch_cells = c(1200, 800),
                     group_prob = c(0.6, 0.4), n_hvgs = 1000,
                     compute_pca = TRUE, pca_ncomp = 10, 
                     output_format = "SingleCellExperiment")

# To compute the desired method
sce <- batchCorrect(input = sce, batch = "Batch", params = HarmonyParams())

# To compute evaluation metrics
metrics <- metrics(input = sce, batch = "Batch",
                   group = "Group", reduction = "harmony", rep = 5)

# To predict the suggested method
pred <- suggested_method(input = sce, batch = "Batch")
```
