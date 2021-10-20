
#' computes jaccard index of GO annotation of a pair of genes, used for an assessment of gene similarity (my approach)
#'
#' @param g.g.rank output of \link{rank_gene_gene_cor} or \link{add_cor_pval}
#' @param martDS \code{dataset} parameter of \code{useMart} funcrion of package "biomaRt"
#' @param go_domain atr which GO domain (out of 'biological_process','molecular_function','cellular_component') to use, default all
#'
#' @return g.g.rank data frame with additional fields \code{jaccard} and \code{shared}, where:
#' - \code{jaccard} jaccard index of GO annotation of a pair of genes
#' - \code{shared} absolute number of shared GO ids between a pair of genes
#'
#' @export


GO_match_top_cor_genes <- function(
  g.g.rank,
  bm,
  go_domain = c('biological_process','molecular_function','cellular_component')
  )
{

  if(go_domain == "all" | is.null(go_domain)){
    go_domain_name <- "all"
    go_domain <- c('biological_process','molecular_function','cellular_component')
  } else {
    go_domain_name <- paste(go_domain, collapse = ",")
  }

  res <- list()
  res[["jaccard"]] <- c()
  res[["shared"]]  <- c()


  g.pairs <- rownames(g.g.rank)
  genes   <- unique(unlist(strsplit(g.pairs, split = "_")))

  bm.BP   <- bm[bm$namespace_1003 %in% go_domain,]

  N <- length(g.pairs)
  for(i in 1:N){
    cur.gene <- strsplit(g.pairs[i], split = "_")[[1]]
    gene_x <- cur.gene[1]
    gene_y <- cur.gene[2]

    GO_gene_x   <- bm.BP$go_id[bm.BP$hgnc_symbol == gene_x]
    GO_gene_y   <- bm.BP$go_id[bm.BP$hgnc_symbol == gene_y]

    if(is.nan(jaccard(GO_gene_x, GO_gene_y)) | is.null(GO_gene_x) | is.null(GO_gene_y)){
      res[["jaccard"]] <- c(res[["jaccard"]], NA)

      res[["shared"]]  <- c(res[["shared"]], NA)
    } else {
      res[["jaccard"]] <- c(res[["jaccard"]], jaccard(GO_gene_x, GO_gene_y))

      res[["shared"]]  <- c(res[["shared"]], length(intersect(GO_gene_x, GO_gene_y)))
    }
  }

  g.g.rank[[paste0("GO_jaccard_", go_domain_name)]]  <- res[["jaccard"]]
  g.g.rank[[paste0("GO_shared_", go_domain_name)]]   <- res[["shared"]]
  return(g.g.rank)
}
