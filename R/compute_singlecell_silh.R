#' Computes silhouette width for a sequence of single-cell clustering
#' @param clustering.list list of clustering result with names of the field being number of clusters and value being clusterting result
#' @param d cell-to-cell distance (dissimilarity)
#'
#' @return names vector with names being number of clusters and values being average silhouette coefficient
#' @export

compute_singlecell_silh <- function(
  clustering.list,
  d
){

  k.seq <- as.numeric(names(clustering.list))
  res   <- c()

  for(k.ch in names(clustering.list)){
    cur.cl <- clustering.list[[k.ch]]

    print(k.ch)
    s <- summary(cluster::silhouette(
      x = cur.cl,
      dist = d
    ))$avg.width

    res <- c(res, s)
  }

  names(res) <- names(clustering.list)
  return(res)
}
