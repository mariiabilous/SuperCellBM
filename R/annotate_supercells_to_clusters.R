#' Computes annotation of super-cell types based on the annotation of single cells
#'
#' Runs \code{\link[SuperCell]{supercell_assign}} for all super-cell elements in SC.list (output of \code{compute_supercells})
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
        annotation  <- SuperCell::supercell_assign(
          clusters = sc.annotation, # single-cell assigment to cell lines (clusters)
          supercell_membership = SC.list$membership, # single-cell assignment to super-cells
          method = annotation.meth)

        SC.list[[meth]][[gamma.ch]][[seed.i.ch]][[annotation.name]] <- annotation
      }
    }
  }

  return(SC.list)
}
