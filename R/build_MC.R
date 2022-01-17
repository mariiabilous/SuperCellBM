#' Build metacell oblect
#' Wrapper that inits and builds metacell object from count matrix
#'
#' @param sc.counts single-cell count matrix (genes x cells)
#' @param proj.name project_name parameter in \link[metacell]
#' @param genes.to.use genes to use ti build metacell
#' @param mc.k.knn k in kNN single-cell network in metacell
#' @param MC.T_vm cutoff to select HVG in metacell
#' @param MC.folder folder to store metacell files
#' @param min_mc_size minimum size of metacell (used to vary graining level)
#' @param p_resamp
#'
#' @return metacell::scdb_mc - like output with some additional information stored as a list
#' @export
#'
build_MC <- function(
  sc.counts,
  proj.name,
  genes.to.use = NULL, # NULL == all genes
  mc.k.knn = 100,
  MC.T_vm = 0.08,
  MC.folder = "MC",
  min_mc_size = 20,
  p_resamp = 0.75, n_resamp = 500, mc.alpha = 2, mc.K = 30,
  verbose = FALSE
){
  SuperCellBM::init_mc(mc.dir = MC.folder)

  mc.bknn  <- SuperCellBM::mc_build_cgraph(
    ge.mtrx = as.matrix(sc.counts),
    genes.to.use = genes.to.use,
    k.knn = mc.k.knn,
    T_vm = MC.T_vm,
    project_name = proj.name,
    balanced.Knn = T,
    return.igraph = F,
    verbose = verbose)

  mc <- mc_build_mc(mc.bknn, min_mc_size = min_mc_size, p_resamp = p_resamp,
                    n_resamp = n_resamp, mc.K = mc.K, mc.alpha = mc.alpha)
}
