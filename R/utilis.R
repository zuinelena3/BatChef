#' Retrive knn index and distances from BBKNN output
#
#' @param dist A dgRMatrix of pairwise distances
#'
knn_index_dist <- function(dist) {
  ncol <- max(dist@p[-1] - dist@p[-length(x = dist@p)])
  nrow = ncol(dist)
  knn_index <- knn_dist <- matrix(0, nrow = nrow, ncol = ncol)
  for (i in seq_len(nrow)) {
    idx <- ((i - 1) * ncol + 1):(i * ncol)
    knn_dist[i, ] <- dist@x[idx][order(dist@x[idx])]
    knn_index[i, ] <- dist@j[idx][order(dist@x[idx])] + 1L
  }
  knn_dist <- cbind(rep(x = 0, nrow), knn_dist)
  knn_index <- cbind(1:nrow, knn_index)
  return(list(idx = knn_index, dist = knn_dist))
}
