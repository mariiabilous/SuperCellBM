#' Remove MC field from SC.list object when not needed anymore (i.e., when super-cell GE has been computed)
#'
#' @param SC.list list of super-cell like strucctures for which mc fields have to be removed
#'
#' @export

remove_mc_fields <- function(
  SC.list
){

  mc.methods <- grep("Metacell", names(SC.list), value = T)

  for(meth in mc.methods){
    for(gamma.ch in names(SC.list[[meth]])){
      for(seed.i.ch in names(SC.list[[meth]][[gamma.ch]])){

        SC.list[[meth]][[gamma.ch]][[seed.i.ch]][["mc2d_name"]] <- NULL
        SC.list[[meth]][[gamma.ch]][[seed.i.ch]][["mc2d"]] <- NULL
        SC.list[[meth]][[gamma.ch]][[seed.i.ch]][["mc_info"]] <- NULL
      }
    }
  }
  return(SC.list)
}
