#' Computes PCA for supee-cells over all simplification methods, graining levels (gammas) and random seeds
#'
#' @param SC.list super-cell list
#' @param SC.GE.list super-cell gene expression list
#' @param N.comp numbe of PC to compute
#' @param genes.omit list of genes to exclude from hvg
#' @param pca_name name of the fiels to add to SC.list
#'
#' @return \code{SC.list} with an additional fiels named \code{pca_name} containing PCA embedding
#' @export

compute_supercells_PCA <- function(
  SC.list,
  SC.GE.list,
  N.comp,
  genes.omit = NULL,
  pca_name = 'SC_PCA'
){

  method.seq   <- names(SC.list)
  gamma.seq    <- as.numeric(names(SC.list[[1]]))
  seed.seq     <- as.numeric(names(SC.list[[2]]))

  hvg.for.pca  <- SC.list[[1]][[1]][[1]]$genes.use

  for(meth in method.seq){
    print(meth)
    cur.gamma.seq <- names(SC.list[[meth]])

    for(gamma.ch in cur.gamma.seq){
      cur.seed.seq <- names(SC.list[[meth]][[gamma.ch]])

      for(seed.i.ch in cur.seed.seq){
        print(paste("Method:", meth, "Gamma:", gamma.ch, "Seed:", seed.i.ch))
        cur.SC     <- SC.list[[meth]][[gamma.ch]][[seed.i.ch]]
        cur.GE     <- SC.GE.list[[meth]][[gamma.ch]][[seed.i.ch]]

        cur.PCA    <-  SuperCell::supercell_prcomp(
          X = Matrix::t(cur.GE[hvg.for.pca,]),
          genes.exclude = genes.omit,
          supercell_size = cur.SC$supercell_size,
          fast.pca = TRUE,
          k = max(N.comp)
        )

        SC.list[[meth]][[gamma.ch]][[seed.i.ch]][[pca_name]] <- cur.PCA
      }
    }
  }
  return(SC.list)
}
