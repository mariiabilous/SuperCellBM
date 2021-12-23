#' Plots optimal number of clusters (as argmax of silhouette coefficient)
#'
#' @param SC.list list of super-cell-like structures
#' @param silh.name name of silhouette field
#' @param sc.silh a vector of single-cell silhouettes with names being number of clusters and values being silhouettes coefficients
#'
#' @return barplot of (median) argmax silhouettes for a certain gamma and SC approach and corresponding data frame
#' @export


plot_argmax_silh <- function(
  SC.list,
  silh.name = "silh:hclust",
  sc.silh = NULL,
  error_bars = c('extr', 'quartiles', 'sd')[1],
  to.save.plot = TRUE,
  to.save.plot.raw = FALSE,
  asp = 0.5,
  fig.folder = './plots',
  fig.name = "",
  ignore.gammas = c(),
  ignore.methods = c(),
  N.clusters.seq.breaks = NULL,
  max.N.cl = NULL,
  bar_width = 0.1,
  dodge_width = 0.12,
  err_bar_line_size = 0.2,
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

  if(fig.name != "") fig.name <- paste0("_", fig.name)

  silh.df <- data.frame()

  for(meth in names(SC.list)){
    for(gamma.ch in names(SC.list[[meth]])){
      gamma = as.numeric(gamma.ch)
      for(seed.i.ch in names(SC.list[[meth]][[gamma.ch]])){
        seed.i <- as.numeric(seed.i.ch)

        cur.SC <- SC.list[[meth]][[gamma.ch]][[seed.i.ch]]

        if(!(silh.name %in% names(cur.SC))){
          stop(paste("Silhouette name:", silh.name, "is not found in super-cell-like structure"))
        }

        if(is.null(max.N.cl)) max.N.cl <- max(as.numeric(names(cur.SC[[silh.name]])))

        cur.silh      <- cur.SC[[silh.name]]
        cur.silh      <- cur.silh[as.numeric(names(cur.silh)) <= max.N.cl]

        cur.silh.max  <- as.numeric(names(cur.silh)[which.max(cur.silh)])

        cur.silh.df <- data.frame(
          Method = meth,
          Gamma = gamma,
          Seed = seed.i,
          argmax.silh = cur.silh.max
        )
        silh.df <- rbind(silh.df, cur.silh.df)
      }
    }
  }

  ## For Gamma == 1
  if(!is.null(sc.silh)){

    if(is.null(max.N.cl)) max.N.cl <- max(as.numeric(names(sc.silh)))

    sc.silh        <- sc.silh[as.numeric(names(sc.silh)) <= max.N.cl]
    sc.silh.argmax <- as.numeric(names(sc.silh)[which.max(sc.silh)])
    cur.silh.df <- data.frame(
      Method = "Exact",
      Gamma = 1,
      Seed = 12345,
      argmax.silh = sc.silh.argmax
    )
    silh.df <- rbind(silh.df, cur.silh.df)

  }

  silh.df.summary <- silh.df %>%
    dplyr::group_by(Method, Gamma) %>%
    dplyr::summarize(
      meanScore   = mean(argmax.silh),
      firstQScore = ifelse(!is.na(mean(argmax.silh)), unname(summary(argmax.silh)[2]), NA),
      thirdQScore = ifelse(!is.na(mean(argmax.silh)), unname(summary(argmax.silh)[5]), NA),
      medianScore = median(argmax.silh),
      sdScore     = ifelse(!is.na(sd(argmax.silh)), sd(argmax.silh), 0),
      minScore     = min(argmax.silh),
      maxScore     = max(argmax.silh),
      medianPsd    = min(median(argmax.silh)+ifelse(!is.na(sd(argmax.silh)), sd(argmax.silh), 0), 1),
      medianMsd    = max(median(argmax.silh)-ifelse(!is.na(sd(argmax.silh)), sd(argmax.silh), 0), 0)
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
  df.to.plot <- silh.df.summary %>%
    dplyr::filter(
      !(Method %in% ignore.methods) & !(Gamma %in% ignore.gammas))

  df.to.plot[['min_err_bar']] <- df.to.plot[[min_err_name]]
  df.to.plot[['max_err_bar']] <- df.to.plot[[max_err_name]]

  .colors <- .colors[unique(df.to.plot$Method)]

  g <- ggplot2::ggplot(df.to.plot,
                       ggplot2::aes(x = Gamma, y = medianScore, fill = Method)) +
    ggplot2::geom_bar(color = NA, alpha = 0.7, orientation = "x", stat = "identity", position = ggplot2::position_dodge(width = dodge_width, preserve = "single"), width = bar_width) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = firstQScore,
                   ymax = thirdQScore,
                   color = Method), size = err_bar_line_size,  position = ggplot2::position_dodge(width = dodge_width, preserve = "single"), width = 0.) +
    #ggplot2::scale_x_log10() +
    ggplot2::labs(x = "Graining level", y = "Median Argmax Silhouette")

  if(!is.null(.colors)){
    g <- g + ggplot2::scale_color_manual(values  = .colors) +
      ggplot2::scale_fill_manual(values = .colors)
  }

  if(!is.null(N.clusters.seq.breaks)){
    g <- g +
      ggplot2::scale_y_continuous(breaks = N.clusters.seq.breaks)
  }

  plot(g)

  if(to.save.plot){
    fig.folder.save = file.path(fig.folder, 'save')
    if(!dir.exists(fig.folder.save))
      dir.create(fig.folder.save, recursive = TRUE)

    filename = paste0('argmax_silhouette_', silh.name, '_errbar_', error_bars, fig.name)
    SCBM_saveplot(p = g, folder = fig.folder.save, name = filename, save.raw.ggplot = FALSE, asp = asp, ...)
  }
  return(invisible(df.to.plot))

}



#' Plots optimal number of clusters (as argmax of silhouette coefficient) for multple clustering algorithms
#'
#' @param SC.list list of super-cell-like structures
#' @param silh.name name of silhouette field
#' @param sc.silh a vector of single-cell silhouettes with names being number of clusters and values being silhouettes coefficients
#'
#' @return barplot of (median) argmax silhouettes for a certain gamma and SC approach and corresponding data frame
#' @export


plot_argmax_silh_alt_clusterings <- function(
  SC.list,
  silh.names = c("silh:hclust", "silh:unw_hclust", "silh:unw_kmeans", "silh:unw_seurat"),
  sc.silh = NULL,
  error_bars = c('extr', 'quartiles', 'sd')[1],
  to.save.plot = TRUE,
  to.save.plot.raw = FALSE,
  asp = 0.5,
  fig.folder = './plots',
  fig.name = "",
  ignore.gammas = c(),
  ignore.methods = c(),
  N.clusters.seq.breaks = NULL,
  bar_width = 0.1,
  dodge_width = 0.12,
  err_bar_line_size = 0.2,
  .shapes = c("hclust" = 1, "unw_hclust" = 2, "unw_kmeans" = 2, "unw_seurat" = 2),
  .colors = c("hclust" = "darkred", "unw_hclust" = "#882255", "unw_kmeans" = "#117733", "unw_seurat" = "#8CC6EC"),
  verbose = FALSE,
  ...
){

  `%>%` <- dplyr::`%>%`

  if(fig.name != "") fig.name <- paste0("_", fig.name)

  silh.df <- data.frame()

  for(meth in names(SC.list)){
    for(gamma.ch in names(SC.list[[meth]])){
      gamma = as.numeric(gamma.ch)
      for(seed.i.ch in names(SC.list[[meth]][[gamma.ch]])){
        seed.i <- as.numeric(seed.i.ch)

        cur.SC <- SC.list[[meth]][[gamma.ch]][[seed.i.ch]]

        for(silh.name in silh.names){
          if(silh.name %in% names(cur.SC)){


          if(silh.name != "silh:unw_seurat"){
            cur.silh      <- cur.SC[[silh.name]]
            cur.silh.max  <- as.numeric(names(cur.silh)[which.max(cur.silh)])
          } else {
            cur.silh.max  <- length(unique(cur.SC[["default_unw_seurat"]]))
          }

          cur.silh.df <- data.frame(
            Method = meth,
            Gamma = gamma,
            Seed = seed.i,
            argmax.silh = cur.silh.max,
            Clustering_name = gsub("silh:", "", silh.name),
            silh.name = silh.name
          )
          silh.df <- rbind(silh.df, cur.silh.df)

          } else {
            warning(paste("Silhouette name:", silh.name, "is not found in super-cell-like structure at", meth, ", Gamma:", gamma.ch, ", Seed:", seed.i.ch))
          }
        }
      }
    }
  }

  ## For Gamma == 1
  if(!is.null(sc.silh)){

    cur.SC <- sc.silh

    for(silh.name in silh.names){
      if(!(silh.name %in% names(cur.SC))){
        stop(paste("Silhouette name:", silh.name, "is not found in super-cell-like structure"))
      }

      if(silh.name != "silh:unw_seurat"){
        cur.silh      <- cur.SC[[silh.name]]
        cur.silh.max  <- as.numeric(names(cur.silh)[which.max(cur.silh)])
      } else {
        cur.silh.max  <- length(unique(cur.SC[["default_unw_seurat"]]))
      }


      cur.silh.df <- data.frame(
        Method = c(unique(silh.df$Method)),
        Gamma = 1,
        Seed = 12345,
        argmax.silh = cur.silh.max,
        Clustering_name = gsub("silh:", "", silh.name),
        silh.name = silh.name
      )
      silh.df <- rbind(silh.df, cur.silh.df)
    }
  }

  silh.df.summary <- silh.df %>%
    dplyr::group_by(Method, Gamma, Clustering_name, silh.name) %>%
    dplyr::summarize(
      meanScore   = mean(argmax.silh),
      firstQScore = ifelse(!is.na(mean(argmax.silh)), unname(summary(argmax.silh)[2]), NA),
      thirdQScore = ifelse(!is.na(mean(argmax.silh)), unname(summary(argmax.silh)[5]), NA),
      medianScore = median(argmax.silh),
      sdScore     = ifelse(!is.na(sd(argmax.silh)), sd(argmax.silh), 0),
      minScore     = min(argmax.silh),
      maxScore     = max(argmax.silh),
      medianPsd    = min(median(argmax.silh)+ifelse(!is.na(sd(argmax.silh)), sd(argmax.silh), 0), 1),
      medianMsd    = max(median(argmax.silh)-ifelse(!is.na(sd(argmax.silh)), sd(argmax.silh), 0), 0)
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
  df.to.plot <- silh.df.summary %>%
    dplyr::filter(
      !(Method %in% ignore.methods) & !(Gamma %in% ignore.gammas))

  df.to.plot[['min_err_bar']] <- df.to.plot[[min_err_name]]
  df.to.plot[['max_err_bar']] <- df.to.plot[[max_err_name]]

  .colors <- .colors[unique(df.to.plot$Clustering_name)]

  for(meth in unique(df.to.plot$Method)){
    g <- ggplot2::ggplot(df.to.plot[df.to.plot$Method == meth, ],
                         ggplot2::aes(x = Gamma, y = medianScore, fill = Clustering_name)) +
      ggplot2::geom_bar(color = NA, alpha = 0.7, orientation = "x", stat = "identity", position = ggplot2::position_dodge(width = dodge_width, preserve = "single"), width = bar_width) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = firstQScore,
                     ymax = thirdQScore,
                     color = Clustering_name), size = err_bar_line_size,  position = ggplot2::position_dodge(width = dodge_width, preserve = "single"), width = 0.) +
      #ggplot2::scale_x_log10() +
      ggplot2::labs(x = "Graining level", y = "Optimal numner of clusters")


    if(!is.null(.colors)){
      g <- g + ggplot2::scale_color_manual(values  = .colors) +
        ggplot2::scale_fill_manual(values = .colors)
    }

    if(!is.null(N.clusters.seq.breaks)){
      g <- g +
        ggplot2::scale_y_continuous(breaks = N.clusters.seq.breaks)
    }

    plot(g)

    if(to.save.plot){
      fig.folder.save = file.path(fig.folder, 'save')
      if(!dir.exists(fig.folder.save))
        dir.create(fig.folder.save, recursive = TRUE)

      filename = paste0('argmax_silhouette_', silh.name, '_errbar_', error_bars, fig.name, "_", meth)
      SCBM_saveplot(p = g, folder = fig.folder.save, name = filename, save.raw.ggplot = FALSE, asp = asp, ...)
    }
  }
  return(invisible(df.to.plot))

}


