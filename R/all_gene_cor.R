#' Computes gene-gene correlations and p values
#'
#' @param ge super-cell gene expression matrix
#' @param cell.size super-cell size
#'
#' @return list with fields:
#'  \code{cell.size} - super-cell size, \code{cor} unweighted correlation (\code{if(return.unwt.cor)}), \code{wt.cor} - weighted correlation betweeen genes, \code{wt.pval} p.values of weighted correlations
#' @export

all_gene_cor_with_pval <- function(ge, cell.size = NULL, return.unwt.cor = FALSE,  mean1 = FALSE){
  if(return.unwt.cor){
    cor          <- cor(Matrix::t(ge))
  } else {
    cor          <- NULL
  }

  w.cor              <- weights::wtd.cor(Matrix::t(ge), weight = cell.size, mean1 = mean1)

  w.cor$correlation[is.na(w.cor$correlation)] <- 0
  w.cor$p.value[is.na(w.cor$p.value)] <- 1

  pval.adj           <- matrix(p.adjust(w.cor$p.value), ncol = ncol(w.cor$p.value))
  colnames(pval.adj) <- colnames(w.cor$p.value)
  rownames(pval.adj) <- rownames(w.cor$p.value)

  res <- list(cell.size    = cell.size,
              cor          = cor,
              wt.cor       = w.cor$correlation,
              wt.pval      = w.cor$p.value,
              wt.adj.pval  = pval.adj
              )
  return(res)
}



#' Computes gene-gene correlations
#'
#' @param ge super-cell gene expression matrix
#' @param cell.size super-cell size
#'
#' @return list with fields:
#'  \code{cell.size} - super-cell size, \code{cor} unweighted correlation (\code{if(return.unwt.cor)}), \code{wt.cor} - weighted correlation betweeen genes, \code{wt.pval} p.values of weighted correlations
#' @export

all_gene_cor <- function(ge, cell.size = NULL, return.unwt.cor = FALSE){
  if(return.unwt.cor){
    cor          <- cor(Matrix::t(ge))
  } else {
    cor <- NULL
  }

  if(!is.null(cell.size)){
    w.cor          <- cov.wt(as.matrix(Matrix::t(ge)), wt = cell.size, cor = TRUE)$cor
  } else {
    w.cor          <- cor(Matrix::t(ge))
  }

  w.cor[is.na(w.cor)] <- 0

  res <- list(cell.size    = cell.size,
              cor          = cor,
              wt.cor       = w.cor)
  return(res)
}
