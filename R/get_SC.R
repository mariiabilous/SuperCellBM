#' Returns super-cell-like structure at specific graining levels
#' @param SC.list super-cell--like structure
#' @param gamma.seq vector of gammas to exptact
#'
#' @return `SC.list` at specific graining levels
#' @export
#'
get_SC_by_gamma <- function(
  SC.list,
  gamma.seq
){
  gamma.seq.ch <- as.character(gamma.seq)

  SC.res <- list()
  for(meth in names(SC.list)){
    SC.res[[meth]] <- list()
    for(gamma.ch in gamma.seq.ch){
      if(gamma.ch %in% names(SC.list[[meth]])){
        SC.res[[meth]][[gamma.ch]] <- SC.list[[meth]][[gamma.ch]]
      }
    }
  }
  return(SC.res)
}
