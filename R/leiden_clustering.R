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
#' @param n_iter Number of iterations of the Leiden clustering algorithm to perform.
#' Positive values above 2 define the total number of iterations to perform,
#' -1 has the algorithm run until it reaches its optimal clustering
#'
#' @export
#' @importFrom scran buildSNNGraph
#' @importFrom leidenAlg leiden.community
#' @importFrom aricode NMI
#'
#' @return A A \link[SingleCellExperiment]{SingleCellExperiment}
#' \link[Seurat]{Seurat} or `AnnData` object that contains the cluster labels.
#' @examples
#'
#' sim <- simulate_data(
#'   n_genes = 1000, batch_cells = c(150, 50),
#'   group_prob = c(0.5, 0.5), n_hvgs = 500,
#'   compute_pca = TRUE, output_format = "SingleCellExperiment"
#' )
#' clust <- leiden_clustering(
#'   input = sim, reduction = "PCA",
#'   nmi_compute = FALSE, resolution = 1, n_iter = 2
#' )
#'
leiden_clustering <- function(input, label_true = NULL, reduction,
                              nmi_compute = TRUE, resolution = NULL, k = 15,
                              store = FALSE, n_iter = -1) {
  sce <- clustInput(input = input, reduction = reduction)
  neighbors <- buildSNNGraph(x = sce, use.dimred = reduction, k = k)

  if (nmi_compute == FALSE) {
    stopifnot("Specify the resolution parameter" = !is.null(resolution))

    clust <- leiden.community(neighbors,
      resolution = resolution,
      n.iterations = n_iter
    )
    clust <- clust$membership
  } else {
    stopifnot("Specify the true label" = !is.null(label_true))

    res <- seq(0.1, 2, by = 0.1)
    max <- 0

    for (i in seq_along(res)) {
      clust <- leiden.community(neighbors,
        resolution = res[i],
        n.iterations = n_iter
      )

      nmi <- NMI(
        c1 = as.vector(colData(sce)[, label_true]),
        c2 = as.vector(clust$membership),
        variant = "sum"
      )

      if (nmi > max) {
        max <- nmi
        names(max) <- res[i]
      }
    }
    clust <- leiden.community(neighbors,
      resolution = as.numeric(names(max)),
      n.iterations = n_iter
    )
    clust <- clust$membership
  }

  if (store == TRUE) {
    ifelse(inherits(input, "AnnDataR6"),
      input$obs$cluster <- clust, input$cluster <- clust
    )
    return(input)
  } else {
    return(clust)
  }
}
