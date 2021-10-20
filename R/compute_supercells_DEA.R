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
  logFC.thresh = 0,
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
          idnt.name <- paste(ident.1, collapse = "_")
          cur.res <- list()

          cur.res[[idnt.name]] <- SuperCell::supercell_FindMarkers(
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


            ## save actual gamma
            if("gamma.actual" %in% names(cur.SC)){
              gamma.actual <- cur.SC$gamma.actual
            } else {
              gamma.actual <- as.numeric(gamma.ch)
            }

            r$gamma.actual <- gamma.actual

            return(r)
          })


        res[[meth]][[gamma.ch]][[seed.i.ch]] <- cur.res
      }
    }
  }

  return(res)
}


#' Computes DEA for single-cell data and for GT cell type annotation
#' @param sc.ge single-cell gene expression
#' @param clusters single-cell clustering result or GT cell type annotation
#' @param ... rest of the parameters of \link{compute_supercells_DEA}
#' @export
compute_singglecell_DEA <- function(
  sc.ge,
  clusters,
  ident.1 = NULL,
  ident.2 = NULL,
  seed = 12345,
  pval.thresh = 0.05,
  ...
){

  if(is.null(ident.1)){ # find all genes
    cur.res <- SuperCell::supercell_FindAllMarkers(
      ge = sc.ge,
      clusters = clusters,
      ...
    )
  } else {
    idnt.name <- paste(ident.1, collapse = "_")
    cur.res <- list()

    cur.res[[idnt.name]] <- SuperCell::supercell_FindMarkers(
      ge = sc.ge,
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

      r$gamma.actual <- 1

      return(r)
    })

  return(cur.res)
}


#' Compute bool vector of cluster-specific DEG
#'
#' @param DEA.GT  GT differentially expressed genes (DEA for GT annotation)
#' @param logFC.thresh thereshold to consider DEG to be statistically significant
#' @param pval.thresh thereshold to consider DEG to be statistically significant (keep the same as in \link{compute_supercells_DEA} and \link{compute_singlecells_DEA})
#'
#' @export

compute_DEA_GT_bool <- function(
  DEA.GT,
  all.genes  = NULL,
  logFC.thresh = 1,
  pval.thresh = 0.05
){
  if(is.null(all.genes)){
    all.genes <- sort(unique(unlist(lapply(DEA.GT, rownames))))
    warning(paste("The entire set of genes was not provided, will use N =", length(all.genes), "genes to compute consistency of DEA"))
  }

  cell.types <- names(DEA.GT)

  DEA.GT.signif <- lapply(DEA.GT, function(x){
    rownames(x)[(x$adj.p.value < pval.thresh) & (x$logFC > logFC.thresh)]
  })

  signif.gene.cl <- c()
  for(cl in names(DEA.GT.signif)){
    signif.gene.cl <- c(signif.gene.cl, paste0(cl, "_", DEA.GT.signif[[cl]]))
  }

  res <- rep(FALSE, length(cell.types)*length(all.genes))
  names(res) <- paste(rep(cell.types, each = length(all.genes)), rep(all.genes, length(cell.types)),  sep = "_")

  res[signif.gene.cl] <- TRUE

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

  `%>%` <- dplyr::`%>%`

  tpr.dea <- data.frame()

  N.clusters    <- length(names(DEA.list)) # GT number of clusters (cell types)
  N.markers     <- sum(GT.DEA.bool)

  for(meth in names(DEA.list)){

    for(gamma.ch in names(DEA.list[[meth]])){
      for(seed.i.ch in names(DEA.list[[meth]][[gamma.ch]])){
        if(verbose) print(paste(meth, "Gamma:", gamma.ch, "Seed:", seed.i.ch))

        cur.DEA           <- DEA.list[[meth]][[gamma.ch]][[seed.i.ch]]

        if(is.null(names(cur.DEA)) & length(cur.DEA) > 0) names(cur.DEA) <- as.character(1:length(cur.DEA))
        markers.logFC_AUC <- c()

        for(cl in names(cur.DEA)){
          cur.DEA.cl      <- cur.DEA[[cl]]

          cur.markers.logFC_AUC         <- cur.DEA.cl$logFC_AUC
          names(cur.markers.logFC_AUC)  <- paste(cl, rownames(cur.DEA.cl), sep = "_")
          markers.logFC_AUC             <- c(markers.logFC_AUC, cur.markers.logFC_AUC)

        }

        if(!is.null(markers.logFC_AUC)){
          markers.logFC_AUC     <- markers.logFC_AUC[names(GT.DEA.bool)]
        } else {
          markers.logFC_AUC     <- rep(NA, length(GT.DEA.bool))
        }

        markers.logFC_AUC[is.na(markers.logFC_AUC)] <- -1000
        names(markers.logFC_AUC)    <- names(GT.DEA.bool)


        cur.rocit <- ROCR::prediction(as.numeric(markers.logFC_AUC), as.numeric(GT.DEA.bool))
        cur.perf  <- ROCR::performance(cur.rocit,"tpr","fpr")

        cur.TPR   <- cur.perf@y.values[[1]]
        cur.FPR   <- cur.perf@x.values[[1]]

        cur.AUC   <- ROCR::performance(cur.rocit,"auc")
        cur.AUC   <- cur.AUC@y.values[[1]]

        cur.F1    <- ROCR::performance(cur.rocit,"f")
        cur.F1    <- cur.F1@y.values[[1]]

        cur.n.pos <- cur.rocit@n.pos.pred[[1]]

        cur.df <- data.frame(
          TPR = cur.TPR,
          FPR = cur.FPR,
          AUC = cur.AUC,
          F1  = cur.F1,
          N.pos = cur.n.pos,
          Gamma = as.numeric(gamma.ch),
          Gamma_actual = unique(cur.DEA.cl$gamma.actual),
          Seed = as.numeric(seed.i.ch),
          Method = meth
        )

        cur.df <- cur.df[-nrow(cur.df),]

        cur.df <- cur.df %>%
          dplyr::filter(N.pos == min(N.markers, max(N.pos)))

        tpr.dea <- rbind(tpr.dea, cur.df)

      }
    }
  }

  #####  for Gamma == 1
  markers.logFC_AUC <- c()
  for(cl in names(cur.DEA)){
    cur.DEA.cl      <- sc.DEA[[cl]]

    cur.markers.logFC_AUC         <- cur.DEA.cl$logFC_AUC
    names(cur.markers.logFC_AUC)  <- paste(cl, rownames(cur.DEA.cl), sep = "_")
    markers.logFC_AUC             <- c(markers.logFC_AUC, cur.markers.logFC_AUC)

  }

  if(!is.null(markers.logFC_AUC)){
    markers.logFC_AUC     <- markers.logFC_AUC[names(GT.DEA.bool)]
  } else {
    markers.logFC_AUC     <- rep(NA, length(GT.DEA.bool))
  }

  markers.logFC_AUC[is.na(markers.logFC_AUC)] <- -1000
  names(markers.logFC_AUC)    <- names(GT.DEA.bool)


  cur.rocit <- ROCR::prediction(as.numeric(markers.logFC_AUC), as.numeric(GT.DEA.bool))
  cur.perf  <- ROCR::performance(cur.rocit,"tpr","fpr")

  cur.TPR   <- cur.perf@y.values[[1]]
  cur.FPR   <- cur.perf@x.values[[1]]

  cur.AUC   <- ROCR::performance(cur.rocit,"auc")
  cur.AUC   <- cur.AUC@y.values[[1]]

  cur.F1    <- ROCR::performance(cur.rocit,"f")
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

    cur.df <- cur.df[-nrow(cur.df),]

    cur.df <- cur.df %>%
      dplyr::filter(N.pos == min(N.markers, max(N.pos)))

    tpr.dea <- rbind(tpr.dea, cur.df)
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
  DEA.consistency.df,
  consistency.index.name = 'TPR',
  error_bars = c('extr', 'quartiles', 'sd')[1],
  SC_meth_exclude = c(),
  to.save.plot = TRUE,
  to.save.plot.raw = FALSE,
  asp = 0.5,
  fig.folder = './plots',
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

  if(consistency.index.name %in% colnames(DEA.consistency.df)){
    DEA.consistency.df[["Score"]] <- DEA.consistency.df[[consistency.index.name]]
  } else {
    stop(paste("consistency.index.name:", consistency.index.name, "not found in clust.consistency.df,
               available values are:", paste(colnames(DEA.consistency.df), collapse = ', ')))
  }

  DEA.consistency.df_summarized <- DEA.consistency.df %>%
    dplyr::group_by(Method, Gamma_actual) %>%
    dplyr::summarise(
      meanScore   = mean(Score),
      firstQScore = unname(summary(Score)[2]),
      thirdQScore = unname(summary(Score)[5]),
      medianScore = median(Score),
      sdScore     = ifelse(!is.na(sd(Score)), sd(Score), 0),
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



  ## Plot across gamma
  df.to.plot <- DEA.consistency.df_summarized %>%
    dplyr::filter(
      !(Method %in% SC_meth_exclude))

  df.to.plot[['min_err_bar']] <- df.to.plot[[min_err_name]]
  df.to.plot[['max_err_bar']] <- df.to.plot[[max_err_name]]


  g <- ggplot2::ggplot(df.to.plot, ggplot2::aes(x = Gamma_actual, y = medianScore, color = Method, fill = Method,  shape = Method)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin=min_err_bar, ymax=max_err_bar), width=.0,
      position = ggplot2::position_dodge(0.02)) +
    ggplot2::scale_color_manual(values = .colors) +
    ggplot2::scale_fill_manual(values = .colors) +
    ggplot2::scale_shape_manual(values = .shapes) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = 'Graining level', y = paste0(consistency.index.name))

  plot(g)

  if(to.save.plot){
    fig.folder.save = file.path(fig.folder, 'save')
    if(!dir.exists(fig.folder.save))
      dir.create(fig.folder.save, recursive = TRUE)

    filename = paste0(consistency.index.name, '_errbar_', error_bars)
    SCBM_saveplot(p = g, folder = fig.folder.save, name = filename, save.raw.ggplot = FALSE, asp = asp, ...)
  }
  return(df.to.plot)

}


