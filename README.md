# BatChef

<img align="left" width="200" src="https://drive.google.com/uc?export=view&id=1_XyW_NBxIw2el05dzVzXCapwr2ktgL4k">

  `BatChef` is an R package that provides a common interface for batch effect 
correction methods in single-cell RNA-seq data.

It allows users to compute various correction methods through a generic 
function by specifying the desired approach. The package contains multiple 
methods implemented in both R and Python, each based on different mathematical approach. 

Also, it provides metrics to evaluate the performance of the correction methods,
such as the Wasserstein distance, Local Inverse Simpson's Index (LISI),
Average Silhouette Width (ASW), and Adjusted Rand Index (ARI).

## Installation

``` r
# install.packages("devtools")
devtools::install_github("zuinelena3/BatChef")
```

## Overview

- `batchCorrect(input, batch, params)` function allows to specify the desired correction method:

  - `input` can be a `SingleCellExperiment`, `Seurat`, and `AnnData` objects;
  
  - `batch` is a string specifying the batch variable
  
  - `parmas` named arguments to pass to individual methods upon dispatch. The values of this parameter are:
  
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
    
    The output can be a `SingleCellExperiment`, `Seurat`, and `AnnData` objects, depending on the class of the input.

- `metrics()` function computes the evaluation metrics, such as Normalized Mutual Information (NMI), Adjusted Rand Index (ARI), Local Inverse Simpson's Index (LISI), Average Silhouette Width (ASW), and Wasserstein distance. The output is a `data.frame` object;

- `simulated_data()` function allows to simulated single-cell RNA-seq data using `Splatter` package, normalize data, select highly variable genes and compute Principal Component Analysis.

## Usage

``` r
library(BatChef)

# Create a simulated dataset with two batches anche two cell types
data <- simulated_data(nGenes = 2000, batchCells = c(1000, 2000), group.prob = c(0.3, 0.7), ncomp = 10, n_hvgs = 1000)

# Batch effects correction with Harmony
out <- batchCorrect(input = data, batch = "Batch", params = HarmonyParams())

# Batch effects correction with Scanorama
out <- batchCorrect(input = out, batch = "Batch", params = ScanoramaParams(assay_type = "logcounts", return_dimred = TRUE))

# Performance evaluation metrics
red <- SingleCellExperiment::reducedDimNames(out)
metrics <- lapply(red, function(x) metrics(input = out, batch = "Batch", group = "Group", reduction = x, rep = 10))
metrics <- do.call(rbind, metrics)
```
