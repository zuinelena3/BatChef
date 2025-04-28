# BatChef

<!-- badges: start -->
  <!-- badges: end -->
  
<img src="https://drive.google.com/uc?export=view&id=1GYlIpgOopyva1zj3CTAYJ-9V6ts7qZwr" alt="image" width="25%" height="auto">
    
## Overview
    
`BatChef` package contains a common interface for single-cell batch correction methods in R and Python environments. This includes methods based on different strategies and a few metrics to evaluate the performance.
  
- `batchCorrect()` function allows the batch effects correction specifying the desired integration method with the `params` parameter. The values of this parameter are:
    
  - `LimmaParams`
  - `CombatParams`
  - `Seuratv3Params`
  - `Seuratv5Params`
  - `FastMNNParams` 
  - `HarmonyParams` 
  - `ScanoramaParams` 
  - `ScVIParams` 
  - `ScMergeParams` 
  - `BBKNNParams`
  - `LigerParams`
  
  The input can be a `SingleCellExperiment`, `Seurat`, and `AnnData` objects.
  The output can be a `SingleCellExperiment`, `Seurat`, and `AnnData` objects, depending on the class of the input.
  
- `metrics()` function computes the evaluation metrics, such as Normalized Mutual Information (NMI), Adjusted Range Index (ARI), Local Inverse Simpson's Index (LISI), Average Silhouette Width (ASW), and Wasserstein distance. The output is a `data.frame` object;

## Installation

``` r
# install.packages("devtools")
devtools::install_github("zuinelena3/BatChef")
```

## Usage

``` r
library(BatChef)

# Create a simulated dataset with two batches anche three cell types

data <- simulated_data(nGenes = 1000, batchCells = c(1000, 2000), group.prob = c(0.3, 0.7), ncomp = 10, n_hvgs = 1000)
```
