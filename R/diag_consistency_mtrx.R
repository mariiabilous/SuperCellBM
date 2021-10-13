
#' Diagonalize two clustering consistency matricies to map clustering results
#' @param m matrix of two clustering consistency (\code{as.matrix(table(clustering, GT_annotation))})
#'
#' @export

diag_consistency_mtrx <- function(m){
  if(!is.matrix(m)){
    stop("m is not a Matrix")
  }
  rnames       <- if(!is.null(rownames(m))){rownames(m)} else{1:nrow(m)}
  cnames       <- if(!is.null(colnames(m))){colnames(m)} else{1:ncol(m)}
  dnames       <- list("1" = rnames, "2" = cnames)
  min.margin   <- which.min(c(nrow(m), ncol(m)))
  max.margin   <- setdiff(1:2, min.margin)
  mnames       <- dnames[[max.margin]]

  ord          <- apply(m, min.margin, which.max)
  ord.val      <- apply(m, min.margin, max)
  if(length(unique(ord)) != length(ord)){
    used.margs   <- unique(ord)
    unused.margs <- setdiff(1:length(mnames), used.margs)

    dupl.num <- which(table(ord) > 1)
    dupl     <- as.numeric(names(dupl.num))
    for(d in dupl){
      #print(d)
      idx        <- which(ord == d)
      #print(idx)
      keep.idx   <- idx[which.max(ord.val[idx])]
      #print(keep.idx)
      rep.idx    <- setdiff(idx, keep.idx)
      #print(rep.idx)
      ord[rep.idx] <- unused.margs[1:length(rep.idx)]
      #print(ord)
      used.margs   <- unique(ord)
      #print(used.margs)
      unused.margs <-  setdiff(1:length(mnames), used.margs)
    }
  }
  #### chech uniquness of $ord$ ####
  res          <- c(mnames[ord], mnames[setdiff(1:length(mnames), ord)])
  new.order    <- c(ord, setdiff(1:length(mnames), ord))
  if(max.margin == 2){ # col
    res <- m[, new.order]
  } else {
    res <- m[new.order,]
  }
  ord.col     <- if(max.margin == 2){new.order} else {1:ncol(m)}
  ord.row     <- if(max.margin == 1){new.order} else {1:nrow(m)}
  res <- list(m = res, vrow = ord.row, vcol = ord.col)
  return(res)
}


#' Map super-cell clustering results to GT cell type annotation
#'
#' @param SC.list list of super-cell like structures
#' @param GT.cell.type vector with GT cell type annotation to which clusering has to be mapped
#' @param clust.name name of clustering that has to be mapped to \code{clust.name} (name of the SC.list field that corresponds to clustering result)
#'
#' @return SC.list with a new field named "\code{paste0(clust.name, "_mapped_to_GT")}" that contained clustering result of \code{clust.name} mapped to \code{GT.cell.type}
#'
#' @export

map_clustering_to_cell_type <- function(
  SC.list,
  GT.cell.type,
  clust.name
){
  GT.cell.type.names <- sort(unique(GT.cell.type))

  for(meth in names(SC.list)){
    for(gamma.ch in names(SC.list[[meth]])){
      for(seed.i.ch in names(SC.list[[meth]][[gamma.ch]])){
        cur.SC <- SC.list[[meth]][[gamma.ch]][[seed.i.ch]]

        if("membership" %in% names(cur.SC)){
          mmbrshp <- cur.SC[["membership"]]
        } else {
          mmbrshp <- 1:length(cur.SC$supercell_size)
        }

        if("cells.use.idx" %in% names(cur.SC)){
          cells.use.idx <- cur.SC[["cells.use.idx"]]
        } else {
          cells.use.idx <- 1:length(mmbrshp)
        }

        if(clust.name %in% names(cur.SC)){
          cur.cl    <- cur.SC[[clust.name]]
        } else {
          stop(paste0("Clustering result '", clust.name, "' is not found in SC.list, please provide a correct clsuterin result name"))
        }

        ## extrapolate SC clusterimg to single cells
        cur.cl.sc <- cur.cl[mmbrshp]
        cur.SC[[paste0(clust.name, "_mapped_to_GT")]] <-
          GT.cell.type.names[diag_consistency_mtrx(m = as.matrix(table(cur.cl.sc, GT.cell.type[cells.use.idx])))$vcol[as.character(cur.cl)]]

        SC.list[[meth]][[gamma.ch]][[seed.i.ch]] <- cur.SC

      }
    }
  }
  return(SC.list)
}
