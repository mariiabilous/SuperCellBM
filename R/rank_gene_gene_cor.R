#' rank gene-gene correlation
#'
#' @param gene.gene.cor gene-gene correlation matrix (square)
#' @param N.top A number of top correlate genes to rank (by default all pairs are ranked)
#'
#' @return top.cor.genes -- ranking of pairwise correlation.
#'
#' @export

rank_gene_gene_cor <- function(gene.gene.cor, N.top = NULL, p.val.thresh = 0.05, rank.absolute = FALSE){

  gg.cor  <- gene.gene.cor$wt.cor
  #gg.pval <- gene.gene.cor$wt.adj.pval

  N.g     <- ncol(gene.gene.cor)

  if(ncol(gg.cor) != nrow(gg.cor)){
    stop("gg.cor has to be a square matrix")
  }


  #print(paste("N.top =", N.top))

  gg.cor[lower.tri(gg.cor, diag = TRUE)] <- NA
  #gg.cor[gg.pval >= p.val.thresh] <- NA # keep statistically significant correlations
  gg.cor <- reshape2::melt(gg.cor, na.rm = TRUE)

  if(is.null(N.top)){
    N.top <- length(gg.cor$value) # all available pairs with statistically significant pvals
  }

  top.cor.genes <- gg.cor[order((gg.cor$value), decreasing = TRUE),][1:N.top,]

  #print(which(is.na(top.cor.genes)))

  top.cor.genes$rank     <- 1:nrow(top.cor.genes)
  top.cor.genes$rank_rev <- nrow(top.cor.genes):1

  if(rank.absolute)
    top.cor.genes$rank.abs <- order(abs(top.cor.genes$value), decreasing = TRUE)

  top.cor.genes$genegene <- paste(top.cor.genes$Var1,
                                  top.cor.genes$Var2, sep = "_")
  rownames(top.cor.genes) <- top.cor.genes$genegene

  return(top.cor.genes)
}

#' computes p value of gene gene correlations
#'
#' @param gg.rank output of \link{rank_gene_gene_cor}
#' @param ge super-cell gene expression matrix
#' @param membership super-cell memberrship vector
#'
#' @return add p value and adjusted p value columns to \code{gg.rank}
#'
#' @export

add_cor_pval <- function(
  gg.rank,
  ge,
  membership = NULL,
  mean1 = FALSE
){
  if(is.null(membership)){

    membership <- 1:ncol(ge)
  }
  pval <- apply(gg.rank, 1, function(x){
    x1 <- ge[ x["Var1"], membership]
    x2 <- ge[ x["Var2"], membership]
    res <- weights::wtd.cor(x = x1, y = x2, mean1 = mean1)[4]
    return(res)
  })

  N.g <- nrow(ge)
  pval.adj <- p.adjust(pval, n = N.g*N.g/2) # adjust for the total number of comparisons

  gg.rank$wt.pval     <- pval
  gg.rank$wt.adj.pval <- pval.adj

  return(gg.rank)
}
