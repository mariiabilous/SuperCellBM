#' Computes annotation of super-cell types based on the annotation of single cells and purity of super-cell
#'
#' Runs \code{\link[SuperCell]{supercell_assign}} and \code{\link[SuperCell]{supercell_purity}} for all super-cell elements in SC.list (output of \code{compute_supercells})
#'
#'
#' @param SC.list list of super-cell-like structures (output of \code{compute_supercells})
#' @param sc.annotation vector of cell type / cluster / etc annotation at the single-cell level (length on the vector == number of single cells)
#' @param annotation.name a name of the annotation, will be added as additional SC.list field
#' @param annotation.meth \code{method} parameter of \link[SuperCell]{supercell_assign} function. Default is 'jaccard'
#' @param verbose flag whether to pring progress steps
#'
#' @return SC.list with each super-cell element containing additional field \code{annotation.name} that contains super-cell annotation
#' @export
#'
#'
annotate_supercells_to_cluster <-function(
  SC.list,
  sc.annotation,
  annotation.name,
  annotation.meth = "jaccard",
  verbose = FALSE
){

  SC_methods <- names(SC.list)

  for(meth in SC_methods){
    if(verbose) print(meth)
    gamma.seq.ch <- names(SC.list[[meth]])

    for(gamma.ch in gamma.seq.ch){
      if(verbose) print(paste("GAMMA:", gamma.ch))
      seed.seq.ch <- names(SC.list[[meth]][[gamma.ch]])

      for(seed.i.ch in seed.seq.ch){
        if(verbose) print(paste("Seed:", seed.i.ch))

        cur.SC <- SC.list[[meth]][[gamma.ch]][[seed.i.ch]]

        if('cells.use.idx' %in% names(cur.SC)){ ## Metacell or Subsampling
          cur_cell_idx <- cur.SC$cells.use.idx
        } else { # Exact, Approx, Random
          cur_cell_idx <- 1:length(cur.SC$membership)
        }

        if(meth != 'Subsampling'){

          annotation <- SuperCell::supercell_assign(
            clusters = sc.annotation[cur_cell_idx], # single-cell assigment to cell lines (clusters)
            supercell_membership = cur.SC$membership[cur_cell_idx], # single-cell assignment to super-cells
            method = annotation.meth)

          purity <- SuperCell::supercell_purity(
            sc.annotation[cur_cell_idx],
            cur.SC$membership[cur_cell_idx])

        } else { # Subsampling: 'super-cell' == single-cell -> super-cell annotation == single-cell annotation
          annotation <- sc.annotation[cur_cell_idx]
          purity <- rep(1, length(annotation))
        }

        SC.list[[meth]][[gamma.ch]][[seed.i.ch]][[annotation.name]] <- annotation
        SC.list[[meth]][[gamma.ch]][[seed.i.ch]][[paste0('purity:', annotation.name)]] <- purity

      }
    }
  }

  return(SC.list)
}



#' @import ggplot2
NULL

#' Plots super-cell purity
#' @param SC.list list of super-cell-like structures (output of \code{compute_supercells})
#' @param annotation.name name of the reference annotation
#'
#' @export

plot_annotation_purity <- function(
  SC.list,
  annotation.name,
  SC_meth_exclude = c('Metacell_default_av', 'Metacell_SC_like_av', 'Subsampling'),
  seed = 12345,
  to.save.plot = TRUE,
  to.save.plot.raw = FALSE,
  asp = 0.5,
  fig.folder = './plots',
  ignore.gammas = c(),
  ignore.methods = c(),
  .shapes = c("Exact"=1, "Approx"=0, "Subsampling"=2, "Random"=3,
             "Metacell_default_fp"=4, "Metacell_default_av" = 8, "Metacell_SC_like_fp"=4, "Metacell_SC_like_av" = 8),
  .colors = c("Exact"="darkred", "Approx"="royalblue", "Subsampling"="black", "Random"="gray",
             "Metacell_default_fp"="forestgreen", "Metacell_default_av" = "forestgreen",
             "Metacell_SC_like_fp"="gold", "Metacell_SC_like_av" = "gold"),
  verbose = FALSE,
  ...
){


  SC_methods <- names(SC.list)

  purity_df <- data.frame()
  purity_key <- paste0('purity:', annotation.name)

  for(meth in SC_methods){
    if(verbose) print(meth)
    gamma.seq.ch <- names(SC.list[[meth]])

    for(gamma.ch in gamma.seq.ch){
      if(verbose) print(paste("GAMMA:", gamma.ch))
      seed.seq.ch <- names(SC.list[[meth]][[gamma.ch]])

      for(seed.i.ch in seed.seq.ch){
        if(verbose) print(paste("Seed:", seed.i.ch))

        cur.SC <- SC.list[[meth]][[gamma.ch]][[seed.i.ch]]
        gamma  <- as.numeric(gamma.ch)
        if ('gamma.actual' %in% names(cur.SC)){ # Metacell
          gamma.actual <- cur.SC[['gamma.actual']]
        } else {
          gamma.actual <- cur.SC[['gamma']]
        }
        seed.i <- as.numeric(seed.i.ch)

        cur_purity           <- cur.SC[[purity_key]]

        cur_df <- data.frame(
          Method = meth,
          Gamma = gamma,
          Gamma_actual = gamma.actual,
          Seed = seed.i,
          Annotation = annotation.name,
          Purity = cur_purity,
          stringsAsFactors = FALSE
        )

        purity_df <- rbind(purity_df, cur_df)
      }
    }
  }


  `%>%` <- dplyr::`%>%`
  # Add fake perfect purity at gamma = 1
  gammas <- sort(unique(purity_df$Gamma))
  gamma1 <- purity_df %>%
    dplyr::filter(Gamma == gammas[1])
  gamma1$Gamma <- 1
  gamma1$Gamma_actual <- 1
  gamma1$Purity <- 1

  purity_df <- rbind(purity_df, gamma1)

  ## Compute summary of the purity
  purity_df_summarized <- purity_df %>%
    dplyr::group_by(Method, Gamma, Gamma_actual, Seed) %>%
    dplyr::summarize(
      meanPurity   = mean(Purity),
      firstQPurity = unname(summary(Purity)[2]),
      thirdQPurity = unname(summary(Purity)[5]),
      medianPurity = median(Purity),
      sdPurity     = sd(Purity))


  ## Plot purity across gamma
  df.to.plot <- purity_df_summarized %>%
    dplyr::filter(
      !(Method %in% SC_meth_exclude) &
        Seed == seed)


  df.to.plot <- df.to.plot %>%
    dplyr::filter(
      !(Method %in% ignore.methods) & !(Gamma %in% ignore.gammas))

  u.Methods <- unique(df.to.plot$Method)

  g <- ggplot2::ggplot(df.to.plot, ggplot2::aes(x = Gamma_actual, y = medianPurity, color = Method, shape = Method)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin=firstQPurity, ymax=thirdQPurity), width=.0,
      position = ggplot2::position_dodge(0.02)) +
    ggplot2::scale_color_manual(values = .colors[u.Methods]) +
    ggplot2::scale_shape_manual(values = .shapes[u.Methods]) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(x = 'Graining level', y = 'Annotation purity')

  plot(g)

  if(to.save.plot){
    fig.folder.save = file.path(fig.folder, 'save')
    if(!dir.exists(fig.folder.save))
      dir.create(fig.folder.save, recursive = TRUE)

    filename = paste0('annotation_purity_', annotation.name)
    SCBM_saveplot(p = g, folder = fig.folder.save, name = filename, save.raw.ggplot = FALSE, asp = asp, ...)
  }
  return(df.to.plot)
}
