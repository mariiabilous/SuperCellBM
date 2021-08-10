
#' Computes metacell from super-cells
#'
#'
#'
#' @export

SCimple2Metacell <- function(X, gamma, min_mc_size = NULL,  p_resamp = 0.75, n_resamp = 500, mc.alpha = 2, mc.K=30){

  N.c <- length(X$cgraph@cell_names)
  N.SC <- round(N.c/gamma)
  # compute min_mc_size such that it approximates N.SC (number of super-cells)
  if(is.null(min_mc_size)){
    min_mc_size <- gamma + 1
  }#
  mc <- mc_build_mc(X, min_mc_size = min_mc_size, p_resamp = p_resamp, n_resamp = n_resamp, mc.K = mc.K, mc.alpha = mc.alpha)
  print(paste("min_mc_size", min_mc_size))

  all.cell.ids      <- X$cgraph@cell_names
  present.cell.ids  <- names(mc$mc@mc)


  membership        <- rep(NA, length = length(all.cell.ids))
  names(membership) <- all.cell.ids
  membership[present.cell.ids] <- mc$mc@mc


  mc2d_name <- paste(mc$mc_name, "_2d")


  metacell::mcell_mc2d_force_knn(mc2d_id = mc2d_name,
                       mc_id = mc$mc_name,
                       graph_id = X$cgraph_name)
  mc2d <- metacell::scdb_mc2d(mc2d_name)


  #### from metacell graph to igraph
  N.mc <- ncol(mc$mc@mc_fp)
  W <- matrix(0, ncol = N.mc, nrow = N.mc)
  rownames(W) <- names(mc2d@mc_x)
  colnames(W) <- names(mc2d@mc_y)
  for(ii in 1:nrow(mc2d@graph)){
    i <- mc2d@graph[ii,1]
    j <- mc2d@graph[ii,2]
    W[i, j] <- 1
  }
  mc.igraph <- igraph::graph_from_adjacency_matrix(W, weighted = T, mode = "undirected", diag = FALSE)

  supercell_size                     <- as.numeric(table(mc$mc@mc))
  igraph::E(mc.igraph)$width         <- sqrt(igraph::E(mc.igraph)$weight/10)
  igraph::V(mc.igraph)$size          <- supercell_size
  igraph::V(mc.igraph)$sizesqrt      <- sqrt(igraph::V(mc.igraph)$size)


  lay              <- cbind(mc2d@mc_x, mc2d@mc_y)
  rownames(lay)    <- NULL
  mc.igraph$layout <- lay


  res <- list(membership            = membership,
              graph.supercells      = mc.igraph,
              supercell_size        = supercell_size,
              simplification.algo   = "Metacell",
              gamma                 = gamma,

              gamma.actual          = N.c/length(table(mc$mc@mc)),
              cells.use.ids         = names(mc$mc@mc),
              cells.use.idx         = pmatch(names(mc$mc@mc), mc$mc@cell_names),
              mc2d_name             = mc2d_name,
              mc2d                  = mc2d,
              mc_info               = mc)

  return(res)
}
