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




#'
#' Plots super-cell purity
#' @param SC.list list of super-cell-like structures (output of \code{compute_supercells})
#' @param annotation.name name of the reference annotation
#'
#' @export

plot_purity <- function(
  SC.list,
  annotation.name,
  verbose = FALSE
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
  return(purity_df)
}
