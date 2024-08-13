# BatChef

<!-- badges: start -->
  <!-- badges: end -->
  
  <img src="https://drive.google.com/uc?export=view&id=1GYlIpgOopyva1zj3CTAYJ-9V6ts7qZwr" alt="image" width="25%" height="auto">
    
## Overview
    
`BatChef` package contains a common interface for single-cell batch correction methods in R and Python environments. This includes methods based on different strategies and a few metrics to evaluate the performance.
  
- `scIntegration` function allows the batch effects correction specifying the desired integration method with the `METHOD` parameter. The values of this parameter are:
    
  - `limmaMethod`
  - `combatMethod`
  - `seuratv3Method`
  - `seuratv5Method`
  - `fastMNNMethod` 
  - `harmonyMethod` 
  - `scanoramaMethod` 
  - `bbknnMethod` 
  - `scVIMethod` 
  
  The output can be a `SingleCellExperiment` or `AltOutput` object, containing three slots: corrected expression matrix, embedding, and metadata;
  
  - `clustering` function computes the Leiden clustering. The resolution parameter is chosen using the Normalized Mutual Information (NMI) metric. The output is a `AnnData` object;
  
  - `metrics` function computes the common quantitative metrics, such as Normalized Mutual Information (NMI), Adjusted Range Index (ARI), Local Inverse Simpson's Index (LISI), Average Silhouette Width (ASW), and Homogeneity. The output is a `data.frame` object;

- `wasserstein_metric` function computes the Wasserstein distance metrics between two batches. We can specify the number of random selection cells and replicates. The output is a `data.frame` object;

## Installation

``` r
# install.packages("pak")
pak::pak("zuinelena3/BatChef")
```

## Usage

``` r
library(BatChef)

# Create a simulated dataset with two batches anche three cell types
data <- simulated_data(nGenes = 1000, batchCells = c(155, 150), group.prob = c(0.3, 0.5, 0.2), ncomp = 10)

# Correct batch effects with Harmony method
harmony_out <- scIntegration(data, batch = "Batch", METHOD = harmonyMethod())

# Compute clustering after Harmony integration
clust_harmony <- clustering(harmony_out, cell_type = "Group", reduction = "harmony")

# Compute common quantitative metrics to evaluate the integration
metrics_harmony <- metrics(anndata = clust_harmony, cell_type = "Group", batch = "Batch", reduction = "harmony")

# Compute the Wasserstein metric by selecting 10 cells randomly twice
wass_harmony <- wasserstein_metric(harmony_out, batch = "Batch", reduction = "HARMONY", n = 10, rep = 2)
```
