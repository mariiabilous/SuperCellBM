#' Computes super-cells for the Exact, Approximate super-cells and for Subsampling and Random grouping
#'
#' Computes super-cells for the Exact, Approximate super-cells and for Subsampling and Random grouping
#' for gammas in \code{"gamma.seq"} and for different Random seeds \code{"seed.seq"}
#'
#'
#' @param sc.GE single-cell gene expression matrix with genes as rows and cells as columns
#' @param gamma.seq a vector of graining levels (i.e., \code{c(1, 5, 10, 50, 100, 200)})
#' @param ToComputeSC a bool flag to indicate whether to compute super-cells (if TRUE) or load (if FALSE) (if was previously computed and saved)
#' @param filename filename to save/load the resulting list of super-cells
#' @param data.folder path to the data folder of a certain analysis, usualy in a format './data/project_name/'
#' @param n.var.genes \code{"n.var.genes"} parameter of SuperCell function \code{\link[SuperCell]{SCimplify()}}
#' @param k.knn \code{"k.knn"} parameter of SuperCell function \code{\link[SuperCell]{SCimplify()}}
#' @param n.pc \code{"n.pc"} parameter of SuperCell function \code{\link[SuperCell]{SCimplify()}}, can be a number, then 1:n.pc PC will be used to compute super-cells or a specific vector of PC to use (e.g., 2:15, c(2, 5, 6, 19))
#' @param approx.N \code{"approx.N"} parameter of SuperCell function \code{\link[SuperCell]{SCimplify()}}
#' @param fast.pca whethre to use fast PCA implementation from ilbra package (default TRUE)
#' @param genes.use \code{"genes.use"} parameter of SuperCell function \code{\link[SuperCell]{SCimplify()}}. A specific genes set to use for computing super-cells, if provided, then parameter \code{"n.var.genes"} will be omitted, if not provided, top \code{"n.var.genes"} most variable genes will be used
#' @param genes.exclude genes exclude from the gene set used for the computation of super-cells (e.g., a vector containing mito genes or ribo genes)
#' @param seed.seq a vector of Random seeds (used for the Approximate coarse-graining, Subsampling and Random grouping)
#'
#' @return list of super-cell in a format SC.list[[SC_method]][[gamma_i]][[seed_i]]
#'
#' @examples
#' \dontrun{
#'
#' }
#' @export
#'
compute_supercells <- function(
  sc.GE,
  gamma.seq,
  ToComputeSC,
  filename = 'initial',
  data.folder = './data',
  n.var.genes = 1000,
  k.knn = 5,
  n.pc = 10,
  approx.N = 1000,
  fast.pca = TRUE,
  genes.use = NULL,
  genes.exclude = NULL,
  seed.seq = c(12345, 111, 19, 42, 7),
  verbose = FALSE
){

  filepath = file.path(data.folder, 'SC')
  filename = paste0(filename, '.Rds')

  if(ToComputeSC){

    SC.list <- list()
    SC.list[['Exact']] = list()
    SC.list[['Approx']] = list()
    SC.list[['Random']] = list()
    SC.list[['Subsampling']] = list()

    for(gamma in gamma.seq){

      if(verbose) print(paste('GAMMMA:', gamma))
      gamma.ch <- as.character(gamma)

      SC.list[['Exact']][[gamma.ch]] = list()
      SC.list[['Approx']][[gamma.ch]] = list()
      SC.list[['Random']][[gamma.ch]] = list()
      SC.list[['Subsampling']][[gamma.ch]] = list()

      ### Exact -- thi one can be replaced with rescaling for gamma > gamma.seq[1]
      if(verbose) print('Exact')
      SC.list[['Exact']][[gamma.ch]][[as.character(seed.seq[1])]] = SuperCell::SCimplify(
        X = sc.GE,
        genes.use = genes.use,
        genes.exclude = genes.exclude,
        n.var.genes = n.var.genes,
        k.knn = k.knn,
        gamma = gamma,
        n.pc = n.pc,
        fast.pca = fast.pca)

      for(seed.i in seed.seq){
        if(verbose) print(paste('SEED:', seed.i))
        seed.i.ch = as.character(seed.i)

        ### Approx
        if(verbose) print('Approx')
        SC.list[['Approx']][[gamma.ch]][[seed.i.ch]] =  SuperCell::SCimplify(
          X = sc.GE,
          genes.use = genes.use,
          genes.exclude = genes.exclude,
          n.var.genes = n.var.genes,
          k.knn = k.knn,
          gamma = gamma,
          n.pc = n.pc,
          fast.pca = fast.pca,
          do.approx = TRUE,
          seed = seed.i,
          approx.N = approx.N)

        ### Random
        if(verbose) print("Random")
        SC.list[['Random']][[gamma.ch]][[seed.i.ch]] <- SCimple2Random(SC = SC.list[['Exact']][[gamma.ch]][[1]],
                                                                       gamma = gamma,
                                                                       seed = seed.i)
        if(verbose) print("Subsampling")
        SC.list[['Subsampling']][[gamma.ch]][[seed.i.ch]] <- SCimple2Subsampling(X = sc.GE,
                                                                                 SC.list[['Exact']][[gamma.ch]][[1]],
                                                                                 gamma = gamma,
                                                                                 seed = seed.i)

      }

    }


    if(!dir.exists(filepath)){
      dir.create(filepath)
    }

    saveRDS(SC.list, file = file.path(filepath, filename))

  } else {
    SC.list <- readRDS(file = file.path(filepath, filename))
  }

  return(SC.list)
}

#' Computed super-cells for the additional gamma as concatenates the results to the existing ones
#'
#' The same parameters as in \code{\link{compute_supercells()}}, with 2 exceptions
#'
#' @param SC.list output of \code{\link{compute_supercells()}} (i.e., list of super-cells in a format \code{"SC.list[['SC_method']][[gamma_i]][[seed_i]]"})
#' @param additional_gamma_seq vector of aggitional gammas for which SC has to be computed
#'
#' @return apdated list of super-cells with additional gammas

compute_supercells_additional_gammas <- function(
  SC.list,
  additional_gamma_seq,
  ToComputeSC,
  data.folder = './data',
  filename = 'additional_gammas',
  approx.N = 1000,
  fast.pca = TRUE,
  verbose = FALSE
){

  filepath = file.path(data.folder, 'SC')
  filename = paste0(filename, '.Rds')

  genes.use <- SC.list[[1]][[1]][[1]]$genes.use
  n.var.genes <- length(genes.use) # redundant when `genes.use` provided
  n.pc      <- SC.list[[1]][[1]][[1]]$n.pc
  k.knn     <- SC.list[[1]][[1]][[1]]$k.knn

  seed.seq  <- as.numeric(names(SC.list[['Approx']][[1]]))
  old_gammas <- as.numeric(names(SC.list[[1]]))

  if(length(intersect(old_gammas, additional_gamma_seq)) > 0){
    common_gammas <- intersect(old_gammas, additional_gamma_seq)
    warning(paste("Gamma(s):", paste(common_gammas, collapse = ","), "already computed. No recalculation will be done!"))
    additional_gamma_seq <- setdiff(additional_gamma_seq, common_gammas)
  }

  if(ToComputeSC){

    for(gamma in additional_gamma_seq){

      print(paste('GAMMMA:', gamma))
      gamma.ch <- as.character(gamma)

      SC.list[['Exact']][[gamma.ch]] = list()
      SC.list[['Approx']][[gamma.ch]] = list()
      SC.list[['Random']][[gamma.ch]] = list()
      SC.list[['Subsampling']][[gamma.ch]] = list()

      ### Exact -- thi one can be replaced with rescaling for gamma > gamma.seq[1]
      if(verbose) print('Exact')
      SC.list[['Exact']][[gamma.ch]][[as.character(seed.seq[1])]] = SuperCell::SCimplify(
        X = sc.GE,
        genes.use = genes.use,
        n.var.genes = n.var.genes,
        k.knn = k.knn,
        gamma = gamma,
        n.pc = n.pc,
        fast.pca = fast.pca)

      for(seed.i in seed.seq){
        if(verbose) print(paste('SEED:', seed.i))
        seed.i.ch = as.character(seed.i)

        ### Approx
        print('Approx')
        SC.list[['Approx']][[gamma.ch]][[seed.i.ch]] =  SuperCell::SCimplify(
          X = sc.GE,
          genes.use = genes.use,
          n.var.genes = n.var.genes,
          k.knn = k.knn,
          gamma = gamma,
          n.pc = n.pc,
          fast.pca = fast.pca,
          do.approx = TRUE,
          seed = seed.i,
          approx.N = approx.N)

        ### Random
        if(verbose) print("Random")
        SC.list[['Random']][[gamma.ch]][[seed.i.ch]] <- SCimple2Random(SC = SC.list[['Exact']][[gamma.ch]][[1]],
                                                                       gamma = gamma,
                                                                       seed = seed.i)
        print("Subsampling")
        SC.list[['Subsampling']][[gamma.ch]][[seed.i.ch]] <- SCimple2Subsampling(X = sc.GE,
                                                                                 SC.list[['Exact']][[gamma.ch]][[1]],
                                                                                 gamma = gamma,
                                                                                 seed = seed.i)

      }

    }


    if(!dir.exists(filepath)){
      dir.create(filepath)
    }

    saveRDS(SC.list, file = file.path(filepath, filename))

  } else {
    SC.list <- readRDS(file = file.path(filepath, filename))
  }

  return(SC.list)

}

