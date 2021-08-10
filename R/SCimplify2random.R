require(SuperCell)
require(igraph)
#' from Scimplify result to random simplification
#'
#' @export

SCimple2Random <- function(SC, gamma, seed = 12345){
  if(is.null(SC$graph.singlecell)){
    stop("SCimply2random is only available for SCiplify result whith return.singlecell.NW == TRUE")
  }

  N.c          <- length(SC$membership)
  N.SC         <- round(N.c/gamma)

  set.seed(seed)
  #r.membership <- sample(round(N.c/gamma), N.c, replace = T)
  r.membership <- c(1:N.SC, sample(N.SC, N.c-N.SC, replace = T))

  SC.random                      <- SC
  SC.random$membership           <- r.membership
  SC.random$supercell_size       <- as.vector(table(r.membership))
  SC.random$simplification.algo  <- "Random_grouping"

  SC.NW                          <- igraph::contract(SC.random$graph.singlecell, r.membership)
  SC.NW                          <- igraph::simplify(SC.NW, remove.loops = T, edge.attr.comb="sum")

  igraph::E(SC.NW)$width         <- sqrt(igraph::E(SC.NW)$weight/10)
  igraph::V(SC.NW)$size          <- SC.random$supercell_size
  igraph::V(SC.NW)$sizesqrt      <- sqrt(igraph::V(SC.NW)$size)

  SC.random$graph.supercells     <- SC.NW
  SC.random$seed                 <- seed

  return(SC.random)
}
