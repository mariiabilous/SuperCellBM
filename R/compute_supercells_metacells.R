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

  gamma.seq.SC <- as.numeric(names(SC.list[[1]]))

  if(setdiff(gamma.seq, gamma.seq.SC)){
    abscent_gammas <- setdiff(gamma.seq, gamma.seq.SC)
    warning(paste("Some gammas in gamma.seq:", paste(abscent_gammas, collapse = ", "),
                  "not found in SC.list, Metacell for these gammas will not be computed!"))
    gamma.seq <- intersect(gamma.seq, gamma.seq.SC)
  }

  SC.gene.used <- SC.list[[1]][[1]][[1]]$genes.use
  seed.ch <- names(SC.list[[1]][[1]])

  for(MC_gene_setting in MC_gene_settings){
    mc.dir <- file.path(data.folder, MC.folder, MC_gene_setting)

    if(ToComputeSC){

      init_mc(mc.dir = mc.dir)

      if(MC_gene_setting == 'Metacell_default'){
        MC.genes.to.use <- NULL
        MC.T_vm <- T_vm_def
      } else {
        MC.genes.to.use <- SC.gene.used
        MC.T_vm <- +Inf #filter out all the genes in order to manually add SC.gene.used
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


#' compute super-cell like structure of metacells (based on provided `min_mc_size` instead of `min_mc_size` as min supercell size)
#'
#' @param sc.counts single-cell count matrix (genes by cells)
#' @param min_mc_size_seq vector of min metacell sizes (to varry graining levels)
#' @param proj.name project name to initiame MC and store MC outputs
#' @param ToComputeSC if TRUE computes metacells, if FALSE, loads from corresponding folder
#' @param genes.use genes used to build super-cells (for MC_gene_settings == 'Metacell_SC_like'), can be ommited if \code{SC.list} is provided
#' @param mc.k.knn metacell k for knn
#' @param T_vm_def metacell cutoff for variable genes
#' @param MC.folder MC folder to store output
#' @param MC_gene_settings do not touch :)
#'
#' @export
compute_supercells_metacells_with_min_mc_size <- function(
  sc.counts,
  min_mc_size_seq,
  proj.name,
  ToComputeSC,
  genes.use = NULL,
  SC.list = NULL,
  mc.k.knn = 100,
  T_vm_def = 0.08,
  MC.folder = "MC",
  MC_gene_settings = c('Metacell_default', 'Metacell_SC_like'),
  verbose = FALSE
){

  mc.bknn.list <- list()
  SC.mc <- list()

  if(!is.null(SC.list)){
    seed.ch <- names(SC.list[[1]][[1]])
  } else {
    seed.ch <- "12345"
  }

  if(is.null(genes.use) & is.null(SC.list)){
    stop("Please provide either genes.use or SC.list parameter!")
  }

  if(is.null(genes.use)){
    SC.gene.used <- SC.list[[1]][[1]][[1]]$genes.use
  } else {
    SC.gene.used <- genes.use
  }


  for(MC_gene_setting in MC_gene_settings){
    mc.dir <- file.path(data.folder, MC.folder, MC_gene_setting)

    if(ToComputeSC){

      init_mc(mc.dir = mc.dir)

      if(MC_gene_setting == 'Metacell_default'){
        MC.genes.to.use <- NULL
        MC.T_vm <- T_vm_def
      } else {
        MC.genes.to.use <- SC.gene.used
        MC.T_vm <- +Inf #filter out all the genes in order to manually add SC.gene.used
      }

      mc.bknn.list[[MC_gene_setting]]  <- mc_build_cgraph(
        ge.mtrx = as.matrix(sc.counts),
        genes.to.use = MC.genes.to.use,
        k.knn = mc.k.knn,
        T_vm = MC.T_vm,
        project_name = proj.name,
        balanced.Knn = T,
        return.igraph = F,
        verbose = verbose)

      SC.mc[[MC_gene_setting]] <- list()

      ## compute metacell partition

      for(min_mc_size in min_mc_size_seq){
        cur.MC <- SCimple2Metacell(
          X = mc.bknn.list[[MC_gene_setting]],
          min_mc_size = min_mc_size)

        gamma.ch <- as.character(cur.MC$gamma.actual)
        SC.mc[[MC_gene_setting]][[gamma.ch]] <- list()

        SC.mc[[MC_gene_setting]][[gamma.ch]][[seed.ch]] <- cur.MC
      }

      SC.mc_save <- SC.mc[[MC_gene_setting]]
      saveRDS(SC.mc_save, file.path(mc.dir, 'SC.Rds'))
    } else {
      SC.mc[[MC_gene_setting]] <- readRDS(file.path(mc.dir, 'SC.Rds'))
    }
  }
  return(SC.mc)
}


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

  return(round(sort(unique(unlist(gamma_actual)))))
}
