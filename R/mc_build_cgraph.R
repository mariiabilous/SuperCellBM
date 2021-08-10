###### write wrappers for metacell package: build balanced or raw knn graph (mc_build_cgraph), comute metacells (mc_build_mc) using output form mc_build_cgraph
require(igraph)
require(metacell)

############### --------------- mc_build_mc -------------- ################

#### Function implements several steps of metacell analysis:
#### 1) creates sc matrix of metacell type
#### 2) builds balanced (bknn) or raw (raw_knn) graph

#'
#'  Builds cell graph object using Metacell package
#'
#'  Builds cell object of Metacell package using UMI data. No gene filtering performed inside,
#'  so either provide \code{ge.mtrx} as a UMI matrix with row containing only genes of interest or specify which genes to use via \code{genes.to.use}
#'
#'  @param ge.mtrx UMI matrix with genes stored as rows, cells stored as cols
#'  @param k.knn k-kearest neighbors parameter  (should  be large (~100) if cgraph is balanced)
#'  @param genes.to.use vector of genes to use for building cgraph. If NULL (default), all the genes from \code{ge.mtrx} will be used
#'  @param T_vm min variance of geness to keep for MC construction (default is -Inf - no MC gene filtering)
#'  @param project_name used for storing and accessing .Rda files of Metacell output
#'  @param balanced.Knn whether cgraph have to be computed using balanced kNN from Metacell package (keep k.knn large then (~100))
#'  @param return.igraph whether to return igraph object of cgraph
#'
#'
#'  @return \code{mat} -- mat object of \code{\link[metacell]{metacell}
#'  @return \code{mat_name} -- name of mat object
#'  @return \code{gstat} -- gstat oject of \code{\link[metacell]{metacell}
#'  @return \code{gstat_name} -- name of gstat object
#'  @return \code{gset} -- gset oject of \code{\link[metacell]{metacell}
#'  @return \code{gset_name} -- name of gset object
#'  @return \code{cgraph} -- cgraph oject of \code{\link[metacell]{metacell}
#'  @return \code{cgraph_name} -- name of cgraph object '
#'  @return \code{ograph} -- \code{igraph} object of cgraph (if return.igraph == TRUE)
#'  @return \code{k.knn} -- \code{k.knn} input parameter
#'  @return \code{balanced.Knn} -- \code{balanced.Knn} input parameter
#'  @return \code{project.name} -- \code{project.name} parameter
#'
#' @export
mc_build_cgraph <- function(ge.mtrx, k.knn, genes.to.use = NULL, T_vm = -Inf, project_name = "metacell", balanced.Knn = TRUE, return.igraph = TRUE){

  if(is.null(genes.to.use)){
      genes    <- rownames(ge.mtrx)
  } else {
      genes    <- genes.to.use[genes.to.use %in% rownames(ge.mtrx)]
  }

  ge.mtrx = ge.mtrx[genes,]
  # save GE data into .tsv file in order to provide yts path into load (mcell_import_scmat_tsv) function
  temp_filename      <- "ge_mtrx.tsv"
  temp_filename_path <- paste0(tempdir(), "/", temp_filename)
  write.table(ge.mtrx, file = temp_filename_path)

  cell.ids <- colnames(ge.mtrx)


  #import GE data
  mat_name <- project_name
  metacell::mcell_import_scmat_tsv(mat_nm = mat_name, fn = temp_filename_path, dset_nm = cell.ids)
  mat <- metacell::scdb_mat(mat_name)
  print("import done")

  # create gene set (for this I have to run gene statistics and "filter" bad genes... although I keep all of them since they have been already filtered)

  gstat_name <- project_name
  metacell::mcell_add_gene_stat(gstat_id = gstat_name, mat_id = mat_name, force=T)
  gstat <- metacell::scdb_gstat(gstat_name)
  print("gene stat done")

  # feature (important genes) set (must be the same as )
  gset_name <- paste0(gstat_name, "_feats")
  metacell::mcell_gset_filter_varmean(gset_id = gset_name, gstat_id = gstat_name, T_vm=T_vm, force_new=T) # no filtering, but this function needed to initialize feats
  gset <- metacell::scdb_gset(gset_name)
  print("gene set done")

  # build cell graph (cgraph) from raw or balanced kNN
  cgraph_name <- paste0(project_name, "_cgraph_from_", if(balanced.Knn){"bknn"} else {"rknn"})
  if(balanced.Knn){
    if(k.knn < 100){
      warning("Knn is too low, unexpected results")
    }
    metacell::mcell_add_cgraph_from_mat_bknn(mat_id = mat_name,
                                   graph_id = cgraph_name,
                                   gset_id = gset_name,
                                   K = k.knn,
                                   dsamp = F)
  } else {
    metacell::mcell_add_cgraph_from_mat_raw_knn(mat_id = mat_name,
                                      graph_id = cgraph_name,
                                      gset_id = gset_name,
                                      K = k.knn,
                                      dsamp = F)
  }

  cgraph <- metacell::scdb_cgraph(cgraph_name)
  print("cgraph done")

  res <- list(mat = mat, mat_name = mat_name,
              gstat = gstat, gstat_name = gstat_name,
              gset = gset, gset_name = gset_name,
              cgraph = cgraph, cgraph_name = cgraph_name)

  if(return.igraph){
    i.ep <- match(cgraph@edges$mc1, cgraph@cell_names)
    j.ep <- match(cgraph@edges$mc2, cgraph@cell_names)
    x.ep <- cgraph@edges$w
    print("done i,j,x")
    W <- Matrix::sparseMatrix(i = i.ep, j = j.ep, x = x.ep)
    print("done W as sparse matrix")
    colnames(W) <- rownames(W) <- cgraph@cell_names
    print("done colnames, rownames of W")

    W <- W[cgraph@nodes, cgraph@nodes] # remove outliers
    print("done remove outliers")

    graph <- igraph::graph_from_adjacency_matrix(W, weighted = T)  #, mode = "undirected")
    print("done graph form adj matrix (sparse)")

    graph <- igraph::as.undirected(graph)
    print("done as.undirected graph")
    graph$node.ids <- cgraph@nodes
    graph$node.idx <- pmatch(graph$node.ids, cgraph@cell_names)
  }

  res$igraph       <- graph
  res$balanced.Knn <- balanced.Knn
  res$k.knn        <- k.knn
  res$project_name <- project_name
  return(res)
}


#'
#'  Builds metacell partition using output of \code{mc_build_cgraph}
#'
#' @param X -- an output of \code{mc_build_cgraph}
#' @param min_mc_size -- \code{min_mc_size} input parameter of \code{mcell_coclust_from_graph_resamp} and
#'  \code{mcell_mc_from_coclust_balanced} from \code{metacell}
#'  @param p_resamp -- \code{mcell_coclust_from_graph_resamp} input parameter of \code{mcell_coclust_from_graph_resamp} from \code{metacell}
#'  @param n_resamp -- \code{n_resamp} input param of \code{mcell_coclust_from_graph_resamp} from \code{metacell}
#'  @param mc.K -- \code{mc.K} input parameter of \code{mcell_mc_from_coclust_balanced} from \code{metacell}
#'  @param mc.alpha -- \code{mc.alpha} input parameter of \code{mcell_mc_from_coclust_balanced} from \code{metacell}
#'
#'  @return \code{coclust} --  coclust oject of \code{\link[metacell]{metacell}
#'  @return \code{coclust_name} -- name of coclust
#'  @return \code{mc} --  mc oject of \code{\link[metacell]{metacell}
#'  @return \code{min_mc_size} -- \code{min_mc_size} input parameter
#'  @return \code{p_resamp} -- \code{p_resamp} input parameter
#'  @return \code{n_resamp} -- \code{n_resamp} input parameter
#'  @return \code{mc.K} -- \code{mc.K} input parameter
#'
#' @export
#'
#'


############### --------------- mc_build_mc -------------- ################
mc_build_mc <- function(X, min_mc_size = 15, p_resamp = 0.75, n_resamp = 500, mc.K = 30, mc.alpha = 2){ # X - output of mc_build_cgraph
  if(is.null(X$project_name)){
    stop("X must be an output of mc_build_cgraph function! X$project_name is missing")
  }
  ######## coclust
  coclust_name <- paste0(X$project_name, "_", n_resamp)
  metacell::mcell_coclust_from_graph_resamp(coc_id = coclust_name,
                                  graph_id = X$cgraph_name,
                                  min_mc_size = min_mc_size,
                                  p_resamp = p_resamp, n_resamp = n_resamp)
  coclust   <- metacell::scdb_coclust(coclust_name)


  ###### metacell
  mc_name <- paste0(X$project_name, "_", n_resamp, "_", min_mc_size)

  metacell::mcell_mc_from_coclust_balanced(
    coc_id = coclust_name,
    mat_id = X$mat_name,
    mc_id  = mc_name,
    K = mc.K, min_mc_size = min_mc_size, alpha = mc.alpha)

  mc <-  metacell::scdb_mc(mc_name)

  res <- list(coclust = coclust,
              coclust_name = coclust_name,
              min_mc_size = min_mc_size,
              p_resamp = p_resamp,
              n_resamp = n_resamp,
              mc.K = mc.K,
              mc = mc,
              mc_name = mc_name)
  return(res)
}



