#' Computes differential expression analysis for super-cells
#'
#' @param SC.list list of super-cells and other simplifications (output of \link{compute_supercell})
#' @param SC.GE.list list of super-cell gene expressions
#' @param cluster.name name of clusterig result field
#' @param ident.1 name(s) of cluster for which markers are computed, if \code{NULL} (default), markers for all clusters will be computed (running \link[SuperCell]{supercell_FindAllMarkers} instead of \link[SuperCell]{supercell_FindMarkers}
#' @param ident.2 name(s) of clusters for comparison. If \code{NULL} (defauld), then all the other clusters used
#' @param seed sandom seed to compute logFC_AUC score
#' @param ... other parameters of \link[SuperCell]{supercell_FindMarkers} or \link[SuperCell]{supercell_FindAllMarkers}
#'
#' @return list of DEA results for each super-cell like structure
#'
#' @export


compute_supercells_DEA <- function(
  SC.list,
  SC.GE.list,
  cluster.name,
  ident.1 = NULL,
  ident.2 = NULL,
  seed = 12345,
  verbose = FALSE,
  pval.thresh = 0.05,
  ...
){

  res <- list()

  method.seq <- names(SC.list)
  for(meth in method.seq){
    res[[meth]] <- list()
    cur.gamma.seq <- names(SC.list[[meth]])

    for(gamma.ch in cur.gamma.seq){
      res[[meth]][[gamma.ch]] <- list()
      cur.seed.seq <- names(SC.list[[meth]][[gamma.ch]])

      for(seed.i.ch in cur.seed.seq){
        if(verbose) print(paste("Method:", meth, "Gamma:", gamma.ch, "Seed:", seed.i.ch))

        cur.SC     <- SC.list[[meth]][[gamma.ch]][[seed.i.ch]]
        cur.GE     <- SC.GE.list[[meth]][[gamma.ch]][[seed.i.ch]]

        if(cluster.name %in% names(cur.SC)){
          clusters   <- cur.SC[[cluster.name]]
        } else {
          stop(paste("Field (cluster.name)", cluster.name, "not found within SC fields, please provide valid clustering name (cluster.name)"))
        }

        if(is.null(ident.1)){ # find all genes
          cur.res <- SuperCell::supercell_FindAllMarkers(
            ge = cur.GE,
            supercell_size = cur.SC$supercell_size,
            clusters = clusters,
            ...
          )
        } else {
          cur.res <- SuperCell::supercell_FindMarkers(
            ge = cur.GE,
            supercell_size = cur.SC$supercell_size,
            clusters = clusters,
            ident.1 = ident.1,
            ident.2 = ident.2,
            ...
          )
        }

        cur.res <- lapply(
          cur.res,
          FUN = function(x){
            r <- x
            mlt <- ifelse(as.numeric(r$adj.p.value) < pval.thresh, 1, 0)

            print(sum(mlt)/length(mlt))

            set.seed(seed)
            r$logFC_AUC <- mlt * r$logFC +
              (1 - mlt)*runif(min = 0, max = 1e-1, n = length(mlt)) * sign(r$logFC)

            return(r)
          })

        ## save actual gamma
        if("gamma.actial" %in% names(cur.SC)){
          gamma.actual <- cur.SC$gamma.actual
        } else {
          gamma.actual <- as.numeric(gamma.ch)
        }

        cur.res$gamma.actual <- gamma.actual

        res[[meth]][[gamma.ch]][[seed.i.ch]] <- cur.res
      }
    }
  }

  return(res)
}



#' Compute consistency of super-cell (meta-cell) DEA with the GT DEA (DEA for GT annotation)
#'
#'@param DEA.list list of DEA result (output of \link{compute_supercells_DEA})
#'@param GT.DEA.bool ground truth differentially expressed genes (named boolean with names beeing genes)
#'@param sc.DEA if NULL (dafault) equalt to \code{GT.DEA.bool}
#'
#'@return TPR, F1, AUC between GT and super-cell DEA
#'
#'@export

compute_consistency_of_supercell_DEA <- function(
  DEA.list,
  GT.DEA.bool,
  sc.DEA = NULL,
  verbose = FALSE
){

  tpr.dea <- data.frame()

  N.clusters    <- length(names(GT.DEA)) # GT number of clusters (cell types)

  for(meth in names(DEA.list)){

    for(gamma.ch in names(DEA.list[[meth]])){
      for(seed.i.ch in names(DEA.list[[meth]][[gamma.ch]])){
        if(verbose) print(paste(meth, "Gamma:", gamma.ch, "Seed:", seed.i.ch))

        cur.DEA           <- DEA.list[[meth]][[gamma.ch]][[seed.i.ch]]

        for(cl in names(cur.DEA)){
          cur.DEA.cl      <- cur.DEA[[cl]]
          cur.GT.DEA.cl   <- GT.DEA[[cl]]
          cur.N.genes     <- length(cur.GT.DEA.cl)

          cur.markers.logFC_AUC     <- cur.DEA.cl$logFC_AUC
          cur.markers.logFC_AUC     <- cur.markers.logFC_AUC[names(GT.DEA.bool)]

          cur.markers.logFC_AUC[is.na(cur.markers.logFC_AUC)] <- -1000
          names(cur.markers.logFC_AUC)    <- names(GT.DEA.bool)

          cur.rocit <- prediction(as.numeric(cur.markers.logFC_AUC), as.numeric(GT.DEA.bool))
          cur.perf  <- performance(cur.rocit,"tpr","fpr")

          cur.TPR   <- cur.perf@y.values[[1]]
          cur.FPR   <- cur.perf@x.values[[1]]

          cur.AUC   <- performance(cur.rocit,"auc")
          cur.AUC   <- cur.AUC@y.values[[1]]

          cur.F1    <- performance(cur.rocit,"f")
          cur.F1    <- cur.F1@y.values[[1]]

          cur.n.pos <- cur.rocit@n.pos.pred[[1]]

          cur.df <- data.frame(
            TPR = cur.TPR,
            FPR = cur.FPR,
            AUC = cur.AUC,
            F1  = cur.F1,
            N.pos = cur.n.pos,
            Gamma = as.numeric(gamma.ch),
            Gamma_actual = cur.DEA$gamma.actual,
            Seed = as.numeric(seed.i.ch),
            Method = meth
          )

          tpr.dea <- rbind(tpr.dea, cur.df)
        }
      }
    }
  }

  #####  for Gamma == 1
  for(cl in names(sc.DEA)){
    cur.DEA.cl      <- sc.DEA[[cl]]
    cur.GT.DEA.cl   <- GT.DEA[[cl]]
    cur.N.genes     <- length(cur.GT.DEA.cl)

    cur.markers.logFC_AUC     <- cur.DEA.cl$logFC_AUC
    cur.markers.logFC_AUC     <- cur.markers.logFC_AUC[names(GT.DEA.bool)]

    cur.markers.logFC_AUC[is.na(cur.markers.logFC_AUC)] <- -1000
    names(cur.markers.logFC_AUC)    <- names(GT.DEA.bool)

    cur.rocit <- prediction(as.numeric(cur.markers.logFC_AUC), as.numeric(GT.DEA.bool))
    cur.perf  <- performance(cur.rocit,"tpr","fpr")

    cur.TPR   <- cur.perf@y.values[[1]]
    cur.FPR   <- cur.perf@x.values[[1]]

    cur.AUC   <- performance(cur.rocit,"auc")
    cur.AUC   <- cur.AUC@y.values[[1]]

    cur.F1    <- performance(cur.rocit,"f")
    cur.F1    <- cur.F1@y.values[[1]]

    cur.n.pos <- cur.rocit@n.pos.pred[[1]]

    for(meth in names(DEA.list)){
      cur.df <- data.frame(
        TPR = cur.TPR,
        FPR = cur.FPR,
        AUC = cur.AUC,
        F1  = cur.F1,
        N.pos = cur.n.pos,
        Gamma = 1,
        Gamma_actual = 1,
        Seed = as.numeric(names(DEA.list[[meth]][[gamma.ch]]))[1],
        Method = meth
      )

      tpr.dea <- rbind(tpr.dea, cur.df)
    }
  }

  return(tpr.dea)
}

####### not finished !!! modoft from here

#' Plot DEA consistency
#'
#' @param clust.consistency.df output of \link{compute_consistency_of_supercell_DEA}
#' @param consistency.index.name name of the consistency index (i.e., TPR, FPR, AUC, F1)
#' @param
#' @param error_bars name of values used for errorbars (for subsampling, random grouping,
#' alternative clusteting of single cells and other methods with more than one clustering/simplification output).
#' \code{'extr'} for min/max, \code{'quartiles'} for quartiles and \code{'sd'} for meadin +- sd
#'
#' @export
#'

plot_DEA_consistency <- function(
  clust.consistency.df,
  consistency.index.name = 'TPR',
  min.value.alt.clustering = 0,
  error_bars = c('extr', 'quartiles', 'sd')[1],
  SC_meth_exclude = c(),
  to.save.plot = TRUE,
  to.save.plot.raw = FALSE,
  asp = 0.5,
  fig.folder = './plots',
  .shapes = c("Exact"=1, "Approximate"=0, "Subsampling"=2, "Random"=3,
              "Metacell_default_fp"=4, "Metacell_default_av" = 8,
              "Metacell_SC_like_fp"=4, "Metacell_SC_like_av" = 8, "Alternative" = 23),
  .colors = c("Exact"="darkred", "Approximate"="royalblue", "Subsampling"="black", "Random"="gray",
              "Metacell_default_fp"="forestgreen", "Metacell_default_av" = "forestgreen",
              "Metacell_SC_like_fp"="gold", "Metacell_SC_like_av" = "gold", "Alternative" = "darkblue"),
  verbose = FALSE,
  ...
){

  `%>%` <- dplyr::`%>%`

  if(consistency.index.name %in% colnames(clust.consistency.df)){
    clust.consistency.df[["Score"]] <- clust.consistency.df[[consistency.index.name]]
  } else {
    stop(paste("consistency.index.name:", consistency.index.name, "not found in clust.consistency.df,
               available values are:", paste(colnames(clust.consistency.df), collapse = ', ')))
  }

  clust_name <- clust.consistency.df[['Clustering_name']][1]

  ## Fiter non-realistic alternative clustering results
  clust.consistency.df <- clust.consistency.df %>%
    dplyr::filter(
      Score >= min.value.alt.clustering | Method != 'Alternative')

  ## Compute summary
  clust.consistency.df_summarized <- clust.consistency.df %>%
    dplyr::group_by(Method, Gamma) %>%
    dplyr::summarize(
      meanScore   = mean(Score),
      firstQScore = unname(summary(Score)[2]),
      thirdQScore = unname(summary(Score)[5]),
      medianScore = median(Score),
      sdScore     = sd(Score),
      minScore     = min(Score),
      maxScore     = max(Score),
      medianPsd    = min(median(Score)+sd(Score), 1),
      medianMsd    = median(Score)-sd(Score))

  clust.consistency.df_summarized[is.na(clust.consistency.df_summarized)] <- 0

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



  ## Plot across gamma
  df.to.plot <- clust.consistency.df_summarized %>%
    dplyr::filter(
      !(Method %in% SC_meth_exclude))

  df.to.plot[['min_err_bar']] <- df.to.plot[[min_err_name]]
  df.to.plot[['max_err_bar']] <- df.to.plot[[max_err_name]]


  g <- ggplot2::ggplot(df.to.plot, ggplot2::aes(x = Gamma, y = medianScore, color = Method, fill = Method,  shape = Method)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin=min_err_bar, ymax=max_err_bar), width=.0,
      position = ggplot2::position_dodge(0.02)) +
    ggplot2::scale_color_manual(values = .colors) +
    ggplot2::scale_fill_manual(values = .colors) +
    ggplot2::scale_shape_manual(values = .shapes) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = 'Graining level', y = paste0(consistency.index.name, ' (', clust_name,')'))

  plot(g)

  if(to.save.plot){
    fig.folder.save = file.path(fig.folder, 'save')
    if(!dir.exists(fig.folder.save))
      dir.create(fig.folder.save, recursive = TRUE)

    filename = paste0(consistency.index.name, '_clustering_', clust_name, '_errbar_', error_bars)
    SCBM_saveplot(p = g, folder = fig.folder.save, name = filename, save.raw.ggplot = FALSE, asp = asp, ...)
  }
  return(df.to.plot)

}


#' Compute alternative clustering at single-cell level
#'
#' @param sc.pca PCA matrix or high dimensional matrix (cells as rows and coordinates/genes as columns)
#' @param N.comp number or vector of PCs to use
#' @param N.clusters.seq vector of number of clusters to compute
#' @param hclust_methods a vector of \code{method} parameter form \link[stats]{hclust}
#' @param other_clust_funcs list of ther clustering functions (in a format of \link[stats]{kmeans})
#'
#' @return list of clusterings \code{res} in a format \code{res[["number_of_clusters"]][["clustering_method[_random_seed]"]]}
#'
#' @export
#'
compute_alternative_clustering <- function(
  sc.pca,
  N.comp = 10,
  N.clusters.seq = c(2:10),
  hclust_methods = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"),
  other_clust_funcs = list("kmeans" = kmeans),
  seed.seq = c(12346, 111, 19, 42, 7)
){
  if(length(N.comp) == 1)
    N.comp = 1:N.comp
  print(N.comp)

  d <- dist(sc.pca[,N.comp])

  sc.clust <- list()
  for(k.ch in as.character(N.clusters.seq)){
    sc.clust[[k.ch]] <- list()
  }

  # different hierarchical clusterings
  for(hclust_meth in hclust_methods){
    cur.clust <- hclust(d = d, method = hclust_meth)
    for(k in N.clusters.seq){
      k.ch <- as.character(k)
      sc.clust[[k.ch]][[hclust_meth]] <- cutree(cur.clust, k = k)
    }
  }

  # k-means with different random seeds
  for(clust_name in names(other_clust_funcs)){
    f <- other_clust_funcs[[clust_name]]
    for(k in N.clusters.seq){
      k.ch <- as.character(k)
      for(seed.i in seed.seq){
        cur.clustering.name <- paste(clust_name, seed.i, sep = "_")
        set.seed(seed.i)

        cur.clust <- kmeans(sc.pca[,N.comp], k)$cluster
        sc.clust[[k.ch]][[cur.clustering.name]] <- cur.clust
      }
    }
  }
  return(sc.clust)
}
