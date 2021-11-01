
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

    if(length(GO_match_sc_cur.jacc.signif)>=N.top){
      GO_match_sc_cur.jacc.signif.top <- GO_match_sc_cur.jacc.signif[1:N.top]
    } else {
      GO_match_sc_cur.jacc.signif.top <- NA
    }


    GO_match_sc_cur.shar.signif     <- GO_match_sc_cur[[go_shared_name]][GO_match_sc_cur[["wt.adj.pval"]] < pval.thresh]
    if(length(GO_match_sc_cur.shar.signif) >= N.top){
      GO_match_sc_cur.shar.signif.top <- GO_match_sc_cur.shar.signif[1:N.top]
    } else {
      GO_match_sc_cur.shar.signif.top <- NA
    }

    s_jaccard <- (as.array(summary(GO_match_sc_cur.jacc.signif.top)))
    s_shared  <- (as.array(summary(GO_match_sc_cur.shar.signif.top)))


    res.cur <- data.frame(
      Meth = "Exact",
      Gamma = 1,
      Gamma_actual = 1,
      Seed = seed,
      GO_jaccard_Mean = as.numeric(unname(s_jaccard[4])),
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

          if(length(GO_match_SC_cur.jacc.signif) >= N.top){
            GO_match_SC_cur.jacc.signif.top <- GO_match_SC_cur.jacc.signif[1:N.top]
          } else {
            GO_match_SC_cur.jacc.signif.top <- NA
          }


          GO_match_SC_cur.shar.signif     <- GO_match_SC_cur[[go_shared_name]][GO_match_SC_cur[["wt.adj.pval"]] < pval.thresh]
          if(length(GO_match_SC_cur.shar.signif) >= N.top){
            GO_match_SC_cur.shar.signif.top <- GO_match_SC_cur.shar.signif[1:N.top]
          } else {
            GO_match_SC_cur.shar.signif.top <- NA
          }


          #if(na_to_0){
          #  GO_match_SC_cur.jacc.signif.top[is.na(GO_match_SC_cur.jacc.signif.top)] <- 0
          #  GO_match_SC_cur.shar.signif.top[is.na(GO_match_SC_cur.shar.signif.top)] <- 0
          #}


          s_jaccard <- (as.array(summary(GO_match_SC_cur.jacc.signif.top)))
          #s_jaccard[is.na(s_jaccard) | is.nan(s_jaccard)] <- 0
          s_shared  <- (as.array(summary(GO_match_SC_cur.shar.signif.top)))
          #s_shared[is.na(s_shared) | is.nan(s_shared)] <- 0

          res.cur <- data.frame(
            Meth = meth,
            Gamma = gamma,
            Gamma_actual = gamma,
            Seed = seed.i,
            GO_jaccard_Mean = as.numeric(unname(s_jaccard[4])),
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


#' Plot (normalized) GO match score
#'
#' @param GO_summary a field of the output of \link{GO_match_summary}
#' @param score.name name of the score (one of the colnames of \code{GO_summary})
#' @param error_bars error bar types fo multiple testing simplifications (i.e., Approx, Subsampling, Random grouping)
#' @param do.normalize whether to normalize GO match scores with a single-cell GO match score
#'
#' @return plots GO match score and returns dataset used for plotting
#'
#' @export
#'

plot_GO_match_score <- function(
  GO_summary,
  score.name = "GO_jaccard_Mean", # other options would be c("GO_jaccard_Median", "GO_shared_Mean", "GO_shared_Median")
  error_bars = c('extr', 'quartiles', 'sd')[1],
  do.normalize = TRUE,
  to.save.plot = TRUE,
  to.save.plot.raw = FALSE,
  asp = 0.5,
  fig.folder = './plots',
  fig.name = "",
  ignore.gammas = c(),
  ignore.methods = c(),
  .shapes = c("Exact"=1, "Approx"=0, "Subsampling"=2, "Random"=3,
              "Metacell_default_fp"=4, "Metacell_default_av" = 8,
              "Metacell_SC_like_fp"=4, "Metacell_SC_like_av" = 8),
  .colors = c("Exact"="darkred", "Approx"="royalblue", "Subsampling"="black", "Random"="gray",
              "Metacell_default_fp"="forestgreen", "Metacell_default_av" = "forestgreen",
              "Metacell_SC_like_fp"="gold", "Metacell_SC_like_av" = "gold"),
  verbose = FALSE,
  ...
){

  `%>%` <- dplyr::`%>%`

  if(!(score.name %in% colnames(GO_summary))){
    stop(paste("score.name:", score.name, "is not in GO_summary data frame, please, provide a correct score name!"))
  }

  GO_summary$Score  <- GO_summary[[score.name]]
  sc.value          <- GO_summary$Score[GO_summary$Gamma == 1][1]

  # noralize with sinle-cell score
  if(do.normalize){
    GO_summary$Score  <- GO_summary$Score/sc.value
    sc.value          <- 1
  }

  GO_summary <- GO_summary %>%
    dplyr::group_by(Meth, Gamma) %>%
    dplyr::summarize(
      meanScore    = unname(summary(Score)[4]),
      firstQScore  = unname(summary(Score)[2]),
      thirdQScore  = unname(summary(Score)[5]),
      medianScore  = unname(summary(Score)[3]),
      sdScore      = ifelse(!is.na(sd(Score)), sd(Score), 0),
      minScore     = min(Score),
      maxScore     = max(Score),
      medianPsd    = min(median(Score)+ifelse(!is.na(sd(Score)), sd(Score), 0), 1),
      medianMsd    = max(median(Score)-ifelse(!is.na(sd(Score)), sd(Score), 0), 0)
    )

  if(is.null(error_bars)){
    error_bars <- 'extr'
  }
  if(!(error_bars %in% c('extr', 'quartiles', 'sd'))){
    stop(paste("Error bar name:", error_bars, "not known. Available names are 'extr', 'quartiles', 'sd' "))
  }


  switch (error_bars,

          extr = {
            min_err_name <- 'minScore'
            max_err_name <- 'maxScore'
          },

          quartiles = {
            min_err_name <- 'firstQScore'
            max_err_name <- 'thirdQScore'
          },

          sd = {
            min_err_name <- 'medianMsd'
            max_err_name <- 'medianPsd'
          },

          {
            min_err_name <- 'minScore'
            max_err_name <- 'maxScore'
          }
  )

  df.to.plot <- GO_summary %>%
    dplyr::filter(
      !(Meth %in% ignore.methods) & !(Gamma %in% ignore.gammas))

  df.to.plot[['min_err_bar']] <- df.to.plot[[min_err_name]]
  df.to.plot[['max_err_bar']] <- df.to.plot[[max_err_name]]

  .colors <- .colors[unique(df.to.plot$Meth)]
  .shapes <- .shapes[unique(df.to.plot$Meth)]

  ymin <- min(df.to.plot$medianScore, na.rm = TRUE) - 0.09
  ymin <- round(ymin, 1)
  ymax <- max(df.to.plot$medianScore, na.rm = TRUE) + 0.09
  ymax <- round(ymax, 1)
  breaks <- seq(ymin, ymax, 0.2)
  if(!(sc.value %in% breaks)) breaks <- breaks + 0.1

  df.to.plot <- df.to.plot[!is.na(df.to.plot$medianScore), ]
  g <- ggplot2::ggplot(df.to.plot,
                       ggplot2::aes(x = Gamma, y = medianScore, color = factor(Meth), shape = factor(Meth))) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin=min_err_bar, ymax=max_err_bar), width=.0,
      position = ggplot2::position_dodge(0.02)) +
    ggplot2::geom_hline(yintercept = 1, color = "gray", linetype = 2) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_continuous(breaks = breaks) +
    ggplot2::scale_color_manual(values = .colors) +
    ggplot2::scale_shape_manual(values = .shapes) +
    ggplot2::theme(asp = 0.5, legend.position = "none") +
    ggplot2::labs(y= score.name, x = "Graining level", title = fig.name)

  plot(g)

  if(to.save.plot){
    fig.folder.save = file.path(fig.folder, 'save', "GO")
    if(!dir.exists(fig.folder.save))
      dir.create(fig.folder.save, recursive = TRUE)

    filename = paste0(score.name, '_errbar_', error_bars, "_", fig.name)
    SCBM_saveplot(p = g, folder = fig.folder.save, name = filename, save.raw.ggplot = FALSE, asp = asp, ...)
  }

  return(invisible(df.to.plot))

}

