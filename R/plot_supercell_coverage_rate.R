### Plot dropout rate distribution
#'
#' @param SC.GE.list list of super-cell gene expression
#' @param methods_to_plot methods to plot
#' @param ignore.gammas graining levels not to plot
#' @param fig.name fig name to add to the standard fig name
#'
#' plots and saves
#'
#' @return data frame with super-cell size distributions
#'
#' @export

plot_supercell_coverage_distribution <- function(
  SC.GE.list,
  sc.GE = NULL,
  methods_to_plot = NULL,
  ignore.gammas = NULL,
  fig.folder = './plots',
  fig.name = "",
  boxplot.width = 1,
  outlier.alpha = 0,
  scale.vln = "width",
  trim.vln = TRUE,
  median.p.size = 1,
  ...
){
  `%>%` <- dplyr::`%>%`

  if(is.null(sc.GE)){
    warning("Single-cell gene expression was not provided (`sc.GE == NULL`), coverage at single-cell level will not be computed")
    abs.sc.coverage <- NA
    sc.coverage     <- NA
  } else {
    abs.sc.coverage <- Matrix::colSums(sc.GE>0)
    sc.coverage     <- abs.sc.coverage/nrow(sc.GE)
  }

  res.df <- data.frame()

  for(meth in names(SC.GE.list)){
    for(gamma.ch in names(SC.GE.list[[meth]])){
      for(seed.i.ch in names(SC.GE.list[[meth]][[gamma.ch]])){
        cur.SC.GE <- SC.GE.list[[meth]][[gamma.ch]][[seed.i.ch]]
        abs_coverage <- Matrix::colSums(cur.SC.GE>0)

        N.g <- nrow(cur.SC.GE)

        cur.df <- data.frame(
          Method = meth,
          Gamma = as.numeric(gamma.ch),
          Seed  = as.numeric(seed.i.ch),
          abs_coverage = abs_coverage,
          coverage = abs_coverage/N.g
        )

        res.df <- rbind(res.df, cur.df)
      }
    }

    ## Gamma == 1
    cur.df <- data.frame(
      Method = meth,
      Gamma = 1,
      Seed  = as.numeric(names(SC.GE.list[[meth]][[gamma.ch]])[1]),
      abs_coverage = abs.sc.coverage,
      coverage = sc.coverage
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
      dplyr::filter(Method == meth & !(Gamma %in% ignore.gammas))

    cur.df$Gamma_f <- as.factor(cur.df$Gamma)

    g <- ggplot2::ggplot(data = cur.df,
                         ggplot2::aes(x = Gamma_f, y = coverage)) +
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
      ggplot2::labs(title = paste0(meth, ", 1-dropout"),
                  x = "Graining level",
                  y = "Coverage")

    plot(g)

    filename = paste0('supercell_coverage_distribution_', fig.name, "_", meth)
    SCBM_saveplot(p = g, folder = fig.folder.save, name = filename, save.raw.ggplot = FALSE, do.log = FALSE, ...)

  }

  return(invisible(res.df))
}
