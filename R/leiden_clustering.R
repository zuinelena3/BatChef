#' Leiden clustering
#'
#' Leiden clustering algorithm.
#'
#' The clustering algorithm can be executed by specifying either a single
#' resolution parameter or range of resolution parameters (from 0.1 to 2).
#' In the case of multiple resolutions, the clustering outcome that corresponds
#' to the highest Normalized Mutual Information (NMI) score is selected. Finding
#' the optimal clustering can be useful to compute the performance evaluation of
#' batch correction methods.
#'
#' @param input A \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object can be supplied.
#' @param label_true A string specifying the ground truth label.
#' @param reduction A string specifying the dimensional reduction
#' on which the clustering analysis will be performed.
#' @param nmi_compute A Boolean value indicating NMI metric calculation to
#' identify the optimal clustering is to be performed (Default: TRUE).
#' @param resolution A numeric value specifying the resolution parameter.
#' @param k An integer scalar specifying the number of nearest neighbors.
#' @param store A Boolean value indicating whether cluster labels are stored
#' within the input object (Default: FALSE).
#'
#' @export
#' @importFrom scran buildSNNGraph
#' @importFrom igraph cluster_leiden membership
#' @importFrom aricode NMI
#'
#' @return A A \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object that contains the cluster labels.
#' @examples
#'
#' sim <- simulate_data(n_genes = 1000, batch_cells = c(150, 50),
#'                      group_prob = c(0.5, 0.5), n_hvgs = 500,
#'                      compute_pca = TRUE, output_format = "SingleCellExperiment")
#' clust <- leiden_clustering(input = sim, reduction = "PCA",
#'                            nmi_compute = FALSE, resolution = 1)
#'
leiden_clustering <- function(input, label_true = NULL, reduction,
                       nmi_compute = TRUE, resolution = NULL, k = 10,
                       store = FALSE) {
  sce <- clustInput(input = input, reduction = reduction)

  neighbors <- buildSNNGraph(x = sce, use.dimred = reduction, k = k)

  if (nmi_compute == FALSE) {
    stopifnot("Specify the resolution parameter" = !is.null(resolution))

    clust <- cluster_leiden(graph = neighbors, resolution = resolution)
    clust <- membership(clust)
  }

  else {
    stopifnot("Specify the true label" = !is.null(label_true))

    res <- seq(0.1, 2, by = 0.1)
    max <- 0

    for (i in seq_along(res)) {
      clust <- cluster_leiden(graph = neighbors, resolution = res[i])

      nmi <- NMI(c1 = as.vector(colData(sce)[, label_true]),
                 c2 = as.vector(membership(clust)),
                 variant = "sum")

      if (nmi > max) {
        max <- nmi
        names(max) <- res[i]
      }
    }
    clust <- cluster_leiden(graph = neighbors, resolution = max)
    clust <- membership(clust)
  }

  if (store == TRUE) {
    ifelse(inherits(input, "AnnDataR6"),
           input$obs$cluster <- clust, input$cluster <- clust)
    return(input)
  }

  else return(clust)
}
