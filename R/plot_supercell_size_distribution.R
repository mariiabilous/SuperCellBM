 #' Plot super-cell size distribution
 #'
 #' @param SC.list list of super-cell objects (result of \list[SuperCell]{SCimplify})
 #' @param methods_to_plot methods to plot
 #' @param ignore.gammas graining levels not to plot
 #' @param fig.name fig name to add to the standard fig name
 #'
 #' plots and saves
 #'
 #' @return data frame with super-cell size distributions
 #'
 #' @export

plot_supercell_size_distribution <- function(
  SC.list,
  methods_to_plot = NULL,
  ignore.gammas = NULL,
  fig.folder = './plots',
  fig.name = "",
  boxplot.width = 1,
  outlier.alpha = 0,
  do.log.y = FALSE,
  scale.vln = "width",
  trim.vln = TRUE,
  median.p.size = 1,
  ...
){

  `%>%` <- dplyr::`%>%`

  res.df <- data.frame()

  for(meth in names(SC.list)){
    for(gamma.ch in names(SC.list[[meth]])){
      for(seed.i.ch in names(SC.list[[meth]][[gamma.ch]])){
        cur.SC <- SC.list[[meth]][[gamma.ch]][[seed.i.ch]]

        cur.df <- data.frame(
          Method = meth,
          Gamma = as.numeric(gamma.ch),
          Seed  = as.numeric(seed.i.ch),
          supercell_size = cur.SC$supercell_size
        )

        res.df <- rbind(res.df, cur.df)
      }
    }

    cur.df <- data.frame(
      Method = meth,
      Gamma = 1,
      Seed  = as.numeric(names(SC.list[[meth]][[gamma.ch]])[1]),
      supercell_size = 1
    )

    res.df <- rbind(res.df, cur.df)
  }

  if(is.null(methods_to_plot)){
    methods_to_plot <- unique(res.df[["Method"]])
  }


  fig.folder.save = file.path(fig.folder, 'save')
  if(!dir.exists(fig.folder.save))
    dir.create(fig.folder.save, recursive = TRUE)

  for(meth in methods_to_plot){
    cur.df <- res.df %>%
      dplyr::filter(Method == meth & !((Gamma) %in% ignore.gammas))

    cur.df$Gamma_f <- as.factor(cur.df$Gamma)

    g <- ggplot2::ggplot(
      data = cur.df,
      ggplot2::aes(x = Gamma_f, y = supercell_size)) +
      ggplot2::geom_violin(
        scale = scale.vln,
        trim = trim.vln,
        ...) +
      ggplot2::geom_boxplot(
        width = boxplot.width,
        outlier.alpha = outlier.alpha,
        fill = "black",
        outlier.colour = NA,
        ...) +
      ggplot2::stat_summary(fun = median, geom = "point", fill = "white", shape = 21, size = median.p.size) +
      ggplot2::labs(title = paste0(meth, ", super-cell size"),
                    x = "Graining level",
                    y = "Super-cell size")

      if(do.log.y){
        g <- g + ggplot2::scale_y_log10()
      }

    plot(g)

    filename = paste0('supercell_size_distribution_', fig.name, "_", meth)
    SCBM_saveplot(p = g, folder = fig.folder.save, name = filename, save.raw.ggplot = FALSE, do.log = FALSE, ...)

  }

  return(invisible(res.df))
}

