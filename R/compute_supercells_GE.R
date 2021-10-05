#' Computes gene expression of metacells
#'
#' @param sc.GE single-cell gene expression
#' @param SC.list list of super-cells (an output of \code{\link{compute_supercells}} or \code{\link{compute_supercells_metacells}})
#' @param ToComputeSC_GE whether to compute or to load already computed GE
#' @param filename file name where to save or load GE
#' @param data.folder path to the data folder, usually has a fprmat './data/project_name'
#'
#' @return a list of GE in the same format as \code{"SC.list"}, i.e., \code{"SC.GE.list[[MS_method]][[gamma_i]][[seed_i]]"}
#'
#' @export
compute_supercells_GE <- function(
  sc.GE,
  SC.list,
  ToComputeSC_GE = TRUE,
  filename = 'additional_gammas',
  data.folder = './data',
  verbose = FALSE
){

  filepath <- file.path(data.folder, 'SC_GE')
  filename <- paste0(filename, '.Rds')

  SC.GE.list <- list()

  SC.methods <- names(SC.list)
  gamma.seq <- as.numeric(names(SC.list[[1]]))
  seed.seq <- as.numeric(names(SC.list[['approx']][[1]]))

  metacell_fp_names = grep('^Metacell_(.*)_fp', names(SC.list), value = TRUE)

  if(ToComputeSC_GE){
    for(SC.meth in SC.methods){

      if(verbose) print(SC.meth)
      SC.GE.list[[SC.meth]] <- list()

      for(gamma.ch  in names(SC.list[[SC.meth]])){
        if(verbose) print(paste("GAMMA:", gamma.ch))
        SC.GE.list[[SC.meth]][[gamma.ch]] <- list()

        for(seed.i.ch in names(SC.list[[SC.meth]][[gamma.ch]])){
          if(verbose) print(paste("Seed:", seed.i.ch))

          # set switcher: 'subsampling' - to subsample GE, 'metacell_fp' - to get metacell footprint as a SC profile, and 'normal' - to average GE within super-cells
          if (SC.meth %in% c('Subsampling', metacell_fp_names)) switcher <- SC.meth else switcher <- 'normal'
          if (switcher %in% metacell_fp_names) switcher <- 'Metacell_fp'

          if(verbose) print(paste("Switcher:", switcher))

          switch (switcher,

                  normal = {
                    if("cells.use.idx" %in% names(SC.list[[SC.meth]][[gamma.ch]][[seed.i.ch]])){
                      cur.idx = SC.list[[SC.meth]][[gamma.ch]][[seed.i.ch]]$cells.use.idx
                    } else {
                      cur.idx <- 1:length(SC.list[[SC.meth]][[gamma.ch]][[seed.i.ch]]$membership)
                    }

                    SC.GE.list[[SC.meth]][[gamma.ch]][[seed.i.ch]] <-
                      supercell_GE(
                        ge = sc.GE[,cur.idx],
                        groups = SC.list[[SC.meth]][[gamma.ch]][[seed.i.ch]]$membership[cur.idx]
                      )
                  },

                  Subsampling = {
                    cur.idx = SC.list[[SC.meth]][[gamma.ch]][[seed.i.ch]]$cells.use.idx
                    SC.GE.list[[SC.meth]][[gamma.ch]][[seed.i.ch]] <- sc.GE[,cur.idx]
                  },

                  Metacell_fp = {
                    fp <- SC.list[[SC.meth]][[gamma.ch]][[seed.i.ch]]$mc_info$mc@mc_fp
                    SC.GE.list[[SC.meth]][[gamma.ch]][[seed.i.ch]] <- fp
                  }
          )
        }
      }
    }

    if(!dir.exists(filepath)){
      dir.create(filepath)
    }

    saveRDS(SC.GE.list, file = file.path(filepath, filename))

  } else {
    SC.GE.list <- readRDS(file = file.path(filepath, filename))
  }
  return(SC.GE.list)
}
