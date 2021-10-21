
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

    jaccard <- function(x, y){
      return(length(intersect(x, y))/length(union(x, y)))
    }

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



#' Computes mean and median values of jaccard GO match score within N.top statistically significant top correlated gene pairs
#'
#' @param GO_match_SC list of \link{GO_match_top_cor_genes} of super-cells
#' @param GO_match_sc list of \link{GO_match_top_cor_genes} of single cells
#' @param N.top number of statistically significant top cor genes pais
#' @param pval.thresh statistical significance
#' @param go_domain GO domain to consider
#' @param na_to_0 whether NA GO match replace with 0
#'
#' @return data frame with mean and median values of jaccard GO match score within N.top statistically significant top correlated gene pairs
#'
#' @export

GO_match_summary <- function(
  GO_match_SC,
  GO_match_sc,
  N.top = 1000,
  pval.thresh = 0.05,
  go_domain = "all",
  seed = 12345,
  na_to_0 = TRUE
){

  go_jaccard_name <- paste0("GO_jaccard_", go_domain)
  go_shared_name  <- paste0("GO_shared_", go_domain)
  res <- list()

  for(cl in names(GO_match_SC)){

    res[[cl]] <- data.frame()

    ### single-cell ###
    GO_match_sc_cur <- GO_match_sc[[cl]]

    if(na_to_0){
      GO_match_sc_cur[[go_jaccard_name]][is.na(GO_match_sc_cur[[go_jaccard_name]])] <- 0
      GO_match_sc_cur[[go_shared_name]][is.na(GO_match_sc_cur[[go_shared_name]])] <- 0
    }

    GO_match_sc_cur.jacc.signif     <- GO_match_sc_cur[[go_jaccard_name]][GO_match_sc_cur[["wt.adj.pval"]] < pval.thresh]
    GO_match_sc_cur.jacc.signif.top <- GO_match_sc_cur.jacc.signif[1:N.top]

    GO_match_sc_cur.shar.signif     <- GO_match_sc_cur[[go_shared_name]][GO_match_sc_cur[["wt.adj.pval"]] < pval.thresh]
    GO_match_sc_cur.shar.signif.top <- GO_match_sc_cur.shar.signif[1:N.top]

    if(na_to_0){
      GO_match_sc_cur.jacc.signif.top[is.na(GO_match_sc_cur.jacc.signif.top)] <- 0
      GO_match_sc_cur.shar.signif.top[is.na(GO_match_sc_cur.shar.signif.top)] <- 0
    }

    s_jaccard <- as.array(summary(GO_match_sc_cur.jacc.signif.top))

    s_shared  <- as.array(summary(GO_match_sc_cur.shar.signif.top))


    res.cur <- data.frame(
      Meth = "Exact",
      Gamma = 1,
      Gamma_actual = 1,
      Seed = seed,
      GO_jaccard_Mean = unname(s_jaccard[4]),
      GO_jaccard_Median = unname(s_jaccard[3]),
      GO_jaccard_N_NAs = unname(ifelse(length(s_jaccard)>=7, s_jaccard[7], 0)),
      GO_shared_Mean = unname(s_shared[4]),
      GO_shared_Median = unname(s_shared[3]),
      GO_shared_N_NAs = unname(ifelse(length(s_shared)>=7, s_shared[7], 0)),
      stringsAsFactors = FALSE
    )

    res[[cl]] <- rbind(res[[cl]], res.cur)

    ### super-cell ###
    for(meth in names(GO_match_SC[[cl]])){

      for(gamma.ch in names(GO_match_SC[[cl]][[meth]])){

        gamma <- as.numeric(gamma.ch)
        for(seed.i.ch in names(GO_match_SC[[cl]][[meth]][[gamma.ch]])){

          seed.i <- as.numeric(seed.i.ch)

          GO_match_SC_cur <- GO_match_SC[[cl]][[meth]][[gamma.ch]][[seed.i.ch]]


          if(na_to_0){
            GO_match_SC_cur[[go_jaccard_name]][is.na(GO_match_SC_cur[[go_jaccard_name]])] <- 0
            GO_match_SC_cur[[go_shared_name]][is.na(GO_match_SC_cur[[go_shared_name]])] <- 0
          }

          GO_match_SC_cur.jacc.signif     <- GO_match_SC_cur[[go_jaccard_name]][GO_match_SC_cur[["wt.adj.pval"]] < pval.thresh]
          GO_match_SC_cur.jacc.signif.top <- GO_match_SC_cur.jacc.signif[1:N.top]

          GO_match_SC_cur.shar.signif     <- GO_match_SC_cur[[go_shared_name]][GO_match_SC_cur[["wt.adj.pval"]] < pval.thresh]
          GO_match_SC_cur.shar.signif.top <- GO_match_SC_cur.shar.signif[1:N.top]

          if(na_to_0){
            GO_match_SC_cur.jacc.signif.top[is.na(GO_match_SC_cur.jacc.signif.top)] <- 0
            GO_match_SC_cur.shar.signif.top[is.na(GO_match_SC_cur.shar.signif.top)] <- 0
          }


          s_jaccard <- as.array(summary(GO_match_SC_cur.jacc.signif.top))
          s_jaccard[is.na(s_jaccard) | is.nan(s_jaccard)] <- 0
          s_shared  <- as.array(summary(GO_match_SC_cur.shar.signif.top))
          s_shared[is.na(s_shared) | is.nan(s_shared)] <- 0

          res.cur <- data.frame(
            Meth = meth,
            Gamma = gamma,
            Gamma_actual = gamma,
            Seed = seed.i,
            GO_jaccard_Mean = unname(s_jaccard[4]),
            GO_jaccard_Median = unname(s_jaccard[3]),
            GO_jaccard_N_NAs = unname(ifelse(length(s_jaccard)>=7, s_jaccard[7], 0)),
            GO_shared_Mean = unname(s_shared[4]),
            GO_shared_Median = unname(s_shared[3]),
            GO_shared_N_NAs = unname(ifelse(length(s_shared)>=7, s_shared[7], 0)),
            stringsAsFactors = FALSE
          )

          res[[cl]] <- rbind(res[[cl]], res.cur)
        }
      }
    }
  }
  return(res)
}
