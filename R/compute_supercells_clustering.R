#' Computes clustering for super-cells
#'
#' @param SC.list list of super-cells and other simplifications (output of \link{compute_supercell})
#' @param N.comp number or a vector of principal components to compute clastering on
#' @param N.clusters.seq vector of number of clusters to compute
#' @param pca_name name of the PCA field, default is 'SC_PCA' as returned with \link{compute_supercells_PCA}
#'
#' @return \code{SC.list} with additional fields:
#'  \code{hclust} with hierarchical clustering results and
#'  \code{silh:hclust} with silhouette coefficients
#' @export

compute_supercells_clustering <- function(
  SC.list,
  N.comp = 10,
  N.clusters.seq = c(2:10),
  pca_name = 'SC_PCA'
){

  if(length(N.comp) == 1)
    N.comp = 1:N.comp
  print(N.comp)

  method.seq <- names(SC.list)
  for(meth in method.seq){

    cur.gamma.seq <- names(SC.list[[meth]])

    for(gamma.ch in cur.gamma.seq){
      cur.seed.seq <- names(SC.list[[meth]][[gamma.ch]])

      for(seed.i.ch in cur.seed.seq){
        print(paste("Method:", meth, "Gamma:", gamma.ch, "Seed:", seed.i.ch))
        cur.SC     <- SC.list[[meth]][[gamma.ch]][[seed.i.ch]]
        cur.dist   <- dist(cur.SC$SC_PCA$x[,N.comp])
        cur.hcl    <- SuperCell::supercell_cluster(D = cur.dist, supercell_size = cur.SC$supercell_size, return.hcl = TRUE,  k = 2)

        SC.list[[meth]][[gamma.ch]][[seed.i.ch]]$hclust <- list()
        SC.list[[meth]][[gamma.ch]][[seed.i.ch]]$hclust_silh <- c()

        for(n.cl.cur in N.clusters.seq){
          n.cl.cur.ch <- as.character(n.cl.cur)
          SC.list[[meth]][[gamma.ch]][[seed.i.ch]]$hclust[[n.cl.cur.ch]]      <- cutree(cur.hcl$hcl, k = n.cl.cur)
          SC.list[[meth]][[gamma.ch]][[seed.i.ch]][['silh:hclust']][n.cl.cur.ch]   <-
            SuperCell::supercell_silhouette(x = SC.list[[meth]][[gamma.ch]][[seed.i.ch]]$hclust[[n.cl.cur.ch]],
                                 dist = cur.dist,
                                 supercell_size = cur.SC$supercell_size)$avg.width
        }
      }
    }
  }

  return(SC.list)
}

#' Compute consistency of super-cell (meta-cell) clustering with the GT annotation (or single-cell clustering, if GT annotation is not available)
#'
#'@param SC.list list of super-cells and other simplifications (output of \link{compute_supercell})
#'@param GT_annotation ground truth clustering of single-cells
#'@param sc.clustering clustering of single-cell data, if not provided (is NULL), will be set to \code{GT_annotation}
#'@param clustering_name names of clustering (default is 'hclust')
#'@param sc.alternative.clustering output of \link{compute_alternative_clustering} or NULL
#'
#'@export

compute_consistency_of_supercell_clustering <- function(
  SC.list,
  sc.annotation,
  sc.clustering = NULL,
  clustering.name = 'hclust',
  sc.alternative.clustering = NULL
){

  clust.comp <- data.frame()

  N.GT.clusters    <- length(unique(sc.annotation)) # GT number of clusters

  for(meth in names(SC.list)){

    for(gamma.ch in names(SC.list[[meth]])){
      for(seed.i.ch in names(SC.list[[meth]][[gamma.ch]])){
        print(paste(meth, "Gamma:", gamma.ch, "Seed:", seed.i.ch))

        cur.SC           <- SC.list[[meth]][[gamma.ch]][[seed.i.ch]]



        ifelse('cells.use.idx' %in% names(cur.SC),
               cells.use.idx <- cur.SC[['cells.use.idx']],
               cells.use.idx <- 1:length(sc.annotation)) # if not all single-cells are present in sinmplified data (subsapling or metacell)

        ifelse('membership' %in% names(cur.SC),
               mmbrshp <- cur.SC[['membership']][cells.use.idx],
               mmbrshp <- 1:length(cells.use.idx)) # if not all single-cells are present in sinmplified data (subsapling or metacell)


        cur.SC.clustering <- cur.SC[[clustering.name]][[as.character(N.GT.clusters)]] # supercell clustering
        cur.SC.clustering.sc <- cur.SC.clustering[mmbrshp]
        cur.sc.annotation    <- sc.annotation[cells.use.idx] # GT single-cell clustering

        cur.clust.comp <- aricode::clustComp(cur.SC.clustering.sc, cur.sc.annotation)

        cur.clust.comp.df   <- data.frame(
          cur.clust.comp,
          Method = meth,
          Gamma = as.numeric(gamma.ch),
          Seed = as.numeric(seed.i.ch),
          Clustering_name = clustering.name,
          stringsAsFactors = FALSE
        )

        clust.comp <- rbind(clust.comp, cur.clust.comp.df)
      }
    }
  }

  ## for Gamma == 1
  if(is.null(sc.clustering)){
    sc.clustering <- sc.annotation
  }
  cur.clust.comp <- aricode::clustComp(sc.clustering, sc.annotation)

  cur.clust.comp.df   <- data.frame(
    cur.clust.comp,
    Method = names(SC.list),
    Gamma = 1,
    Seed = as.numeric(names(SC.list[[1]][[1]])[1]),
    Clustering_name = clustering.name,
    stringsAsFactors = FALSE
  )
  clust.comp <- rbind(clust.comp, cur.clust.comp.df)


  ## alternative clustering for single-cell data

  if(!is.null(sc.alternative.clustering)){

    if(!(as.character(N.GT.clusters) %in% names(sc.alternative.clustering))){
      stop(paste("Clustering for", N.GT.clusters, "is not providede (", N.GT.clusters, "is not in names(sc.alternative.clustering))."))
    }

    sc.alternative.clustering.Ncl <- sc.alternative.clustering[[as.character(N.GT.clusters)]]

    for(alt.cl.name in names(sc.alternative.clustering.Ncl)){
      cur.clust.comp <- aricode::clustComp(sc.alternative.clustering.Ncl[[alt.cl.name]], sc.annotation)

      cur.clust.comp.df   <- data.frame(
        cur.clust.comp,
        Method = 'Alternative',
        Gamma = 1,
        Seed = as.numeric(names(SC.list[[1]][[1]])[1]),
        Clustering_name = alt.cl.name,
        stringsAsFactors = FALSE
      )
      clust.comp <- rbind(clust.comp, cur.clust.comp.df)
    }
  }

  return(clust.comp)
}


#' Plot clustering consistency
#'
#' @param clust.consistency.df output of \link{compute_consistency_of_supercell_clustering}
#' @param consistency.index.name name of the consistency index (output of \link[aricode]{clustComp})
#' @param min.value.alt.clustering min index value for the alternative clustering consistency
#' @param error_bars name of values used for errorbars (for subsampling, random grouping,
#' alternative clusteting of single cells and other methods with more than one clustering/simplification output).
#' \code{'extr'} for min/max, \code{'quartiles'} for quartiles and \code{'sd'} for meadin +- sd
#'
#' @export
#'

plot_clustering_consistency <- function(
  clust.consistency.df,
  consistency.index.name = 'ARI',
  min.value.alt.clustering = 0,
  error_bars = c('extr', 'quartiles', 'sd')[1],
  SC_meth_exclude = c(),
  to.save.plot = TRUE,
  to.save.plot.raw = FALSE,
  asp = 0.5,
  fig.folder = './plots',
  .shapes = c("Exact"=1, "Approx"=0, "Subsampling"=2, "Random"=3,
              "Metacell_default_fp"=4, "Metacell_default_av" = 8,
              "Metacell_SC_like_fp"=4, "Metacell_SC_like_av" = 8, "Alternative" = 23),
  .colors = c("Exact"="darkred", "Approx"="royalblue", "Subsampling"="black", "Random"="gray",
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

  sc.clust <- list()
  for(k.ch in as.character(N.clusters.seq)){
    sc.clust[[k.ch]] <- list()
  }

  if(length(hclust_methods) > 0){
    d <- dist(sc.pca[,N.comp])
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
