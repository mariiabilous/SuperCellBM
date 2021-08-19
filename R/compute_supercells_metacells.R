#' compute super-cell like structure of metacells
#'
#' @param ToComputeSC if TRUE computes metacells, if FALSE, loads from corresponding folder
#' @export
compute_supercells_metacells <- function(
  sc.counts,
  gamma.seq,
  SC.list,
  proj.name,
  ToComputeSC,
  mc.k.knn = 100,
  T_vm_def = 0.08,
  MC.folder = "MC",
  MC_gene_settings = c('Metacell_default', 'Metacell_SC_like'),
  verbose = FALSE
){

  mc.bknn.list <- list()
  SC.mc <- list()

  SC.gene.used = SC.list[[1]][[1]][[1]]$genes.use
  seed.ch = names(SC.list[[1]][[1]])

  for(MC_gene_setting in MC_gene_settings){
    mc.dir <- file.path(data.folder, MC.folder, MC_gene_setting)

    if(ToComputeSC){

      init_mc(mc.dir = mc.dir)

      if(MC_gene_setting == 'Metacell_default'){
        MC.genes.to.use <- NULL
        MC.T_vm <- T_vm_def
      } else {
        MC.genes.to.use <- SC.gene.used
        MC.T_vm <- -Inf
      }

      mc.bknn.list[[MC_gene_setting]]  <- mc_build_cgraph(
        ge.mtrx = as.matrix(sc.counts),
        genes.to.use = MC.genes.to.use,
        k.knn = mc.k.knn,
        T_vm = MC.T_vm,
        project_name = proj.name,
        balanced.Knn = T,
        return.igraph = T,
        verbose = verbose)

      SC.mc[[MC_gene_setting]] <- list()

      ## compute metacell partition
      for(gamma in gamma.seq){
        gamma.ch <- as.character(gamma)
        SC.mc[[MC_gene_setting]][[gamma.ch]] <- list()

        SC.mc[[MC_gene_setting]][[gamma.ch]][[seed.ch]] <- SCimple2Metacell(
          X = mc.bknn.list[[MC_gene_setting]],
          gamma = gamma,
          min_mc_size = min(SC.list[['Exact']][[gamma.ch]][[1]]$supercell_size))
      }

      SC.mc_save <- SC.mc[[MC_gene_setting]]
      saveRDS(SC.mc_save, file.path(mc.dir, 'SC.Rds'))
    } else {
      SC.mc[[MC_gene_setting]] <- readRDS(file.path(mc.dir, 'SC.Rds'))
    }
  }
  return(SC.mc)
}

SC.list$Exact$`5`$`12345`
#' returns list of actual gammas of metacells
#'
#' @param SC.mc the output of \code{\link{compute_supercells_metacells}}
#'
#' @return sorted list of rounded to int actuall gammas
#' @export
#'
get_actual_gammas_metacell <- function(
  SC.mc
){
  gamma_actual <- list()
  for(mc_type in names(SC.mc)){
    gamma_actual[[mc_type]] <- c()
    for (gamma.ch in names(SC.mc[[mc_type]])){
      for(seed.i.ch in names(SC.mc[[mc_type]][[gamma.ch]])){
        gamma_actual[[mc_type]] <- c(gamma_actual[[mc_type]], SC.mc[[mc_type]][[gamma.ch]][[seed.i.ch]]$gamma.actual)
      }
    }
  }

  return(round(sort(unique(unlist(gamma_actual[[mc_type]])))))
}
