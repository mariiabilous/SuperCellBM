#' Concatenates super-cell-like structures by method
#'
#' @param SC.list1
#' @param SC.list2
#' @return concatenated `SC.list`
#' @export

bind_by_method <- function(
  SC.list1,
  SC.list2
){
  meth1 <- names(SC.list1)
  meth2 <- names(SC.list2)

  if(length(intersect(meth1, meth2))>0){
    stop(paste("Don't know how to concatenate when `SC.list1` and `SC.list2` have common method(s)", paste(intersect(meth1, meth2), collapse = ", ")))
  }

  SC.list <- c(SC.list1, SC.list2)
  return(SC.list)
}



#' Concatenates super-cell-like structures by graining level
#' @param SC.list1
#' @param SC.list2
#'
#' @return concatenated `SC.list`
#' @export

bind_by_graining_level <- function(
  SC.list1,
  SC.list2
){
  SC.list <- SC.list1

  for(meth in names(SC.list1)){
      if(meth %in% names(SC.list2)){
        gammas1 <- names(SC.list1[[meth]])
        gammas2 <- names(SC.list2[[meth]])

        if(length(intersect(gammas1, gammas2))){
          stop(paste("Don't know how to concatenate when `SC.list1` and `SC.list2` have common gamma(s)", paste(intersect(gammas1, gammas2), collapse = ", ")))
        }

        for(gamma.ch in gammas2){
          SC.list[[meth]][[gamma.ch]] <- SC.list2[[meth]][[gamma.ch]]
        }
      }
  }
return(SC.list)
}


#' Concatenates super-cell-like structures by seeds
#' @param SC.list1
#' @param SC.list2
#'
#' @return concatenated `SC.list`
#' @export

bind_by_seed <- function(
  SC.list1,
  SC.list2
){
  SC.list <- SC.list1

  for(meth in names(SC.list1)){
    if(meth %in% names(SC.list2)){

      for(gamma.ch in names(SC.list1[[meth]])){
        if(gamma.ch %in% names(SC.list2[[meth]])){

          seeds1 <- names(SC.list1[[meth]][[gamma.ch]])
          seeds2 <- names(SC.list2[[meth]][[gamma.ch]])
          if(length(intersect(seeds1, seeds2))){
            stop(paste("Don't know how to concatenate when `SC.list1` and `SC.list2` have common seeds(s)", paste(intersect(seeds1, seeds2), collapse = ", ")))
          }

          for(seed.i.ch in seeds2){
            SC.list[[meth]][[gamma.ch]][[seed.i.ch]] <- SC.list2[[meth]][[gamma.ch]][[seed.i.ch]]
          }
        }
      }
    }
  }
  return(SC.list)
}
