#' Convert \list[SeuratObject]{Graph} to \list[igraph]{igraph}
#'
#' @param graph \list[SeuratObject]{Graph}
#'
#' @return \list[igraph]{igraph}
#' @export

Graph_to_igraph <- function(
 graph,
 ...
){

  adj.matrix <- Matrix::sparseMatrix(i = graph@i+1, p = graph@p, x = graph@x, dims = graph@Dim, dimnames = graph@Dimnames)
  igr        <- igraph::graph_from_adjacency_matrix(
    adjmatrix = adj.matrix,
    weighted = TRUE,
    diag = FALSE,
    ...
  )

  igr        <- igraph::simplify(igr, remove.multiple = T)

  return(igr)
}

#' Builds super-cell object from a single-cell NW and a membership vector
#'
#' @param sc.nw igraph single-cell network
#' @param membership super-cell membership vector
#' @param gname graph name
#' @param resolution clustering resolution (for louvain clustering)
#' @param N.comp number (or vector) of principal components used to build super-cells
#' @param genes.use genes used to build super-cells
#' @param return.singlecell.NW whether to return single-cell network
#'
#' @return the same as \link[SuperCell]{SCimplify} with an additional field - `res` for the resolution usign in \link[Seurat]{FindClusters} to obtain super-cell membership
#'
#' @export

SC_from_membership <- function(
  sc.nw,
  membership,
  gname,
  resolution = NULL,
  N.comp = NULL,
  genes.use = NULL,
  return.singlecell.NW = FALSE
){
  supercell_size   <- as.vector(table(membership))
  N.SC             <- length(unique(membership))

  actual.gamma     <- round(length(membership)/N.SC)
  if(actual.gamma == 1 & length(membership)/N.SC != 1)
    actual.gamma <- actual.gamma + 1
  actual.gamma.ch  <- as.character(actual.gamma)

  # build super-cell network by merging single-cell NW nodes into super-cells
  SC.NW            <- igraph::contract(sc.nw, membership)
  SC.NW            <- igraph::simplify(SC.NW, remove.loops = T, edge.attr.comb="sum")

  igraph::E(SC.NW)$width         <- sqrt(igraph::E(SC.NW)$weight/10)
  igraph::V(SC.NW)$size          <- supercell_size
  igraph::V(SC.NW)$sizesqrt      <- sqrt(igraph::V(SC.NW)$size)


  res <- list(
    graph.supercells = SC.NW,
    gamma = actual.gamma,
    resolution = resolution,
    N.SC = N.SC,
    membership = membership,
    supercell_size = supercell_size,
    genes.use = genes.use,
    simplification.algo = gname,
    do.approx = FALSE,
    n.pc = N.comp,
    k.knn = tidyr::extract_numeric(gname)
  )

  if(return.singlecell.NW){res$graph.singlecell <- sc.nw}
  return(res)
}


#' Compute super-cells using alternative algorithms behind.
#' Compute super-cells trying different number of neighbors for KNN, SNN implements in Seurat and also different clustering algorithms (i.e., current walktrap and Seurat-based Louvain)
#'
#' @param counts count matrix (genes x cells)
#' @param genes.use a set of genes used in the original super-cell construction (output of \link[SuperCell]{SCimplify}, field `genes.use`)
#' @param meta.data meta data as an input to \link[SeuratObject]{CreateSeuratObject}
#' @param ge log-normalized gene expression matrix (if computation is different from whar Seurat outputs)
#' @param k.seq set of k to compute knn network
#' @param res.seq set of resolutions to obtain different graining levels using louvain clustering in \link[Seurat]{FindClusters}
#' @param gamma.seq graining levels to compute for walktrap clustering (is NULL, graining levels are retreived from the graining levels obtained in louvain clustering with \code{res.seq} resolutions)
#' @param return.singlecell.NW whether to return single-cell network
#' @param common.gammas whether to compute walktrap super-cells for both provided graining leveles (if providede) and obtained graning levels with louvain


#' @return list of SuperCell -like structures with the first result layer corresponding to a NW type and the second layer corresponding to the clustering type (walktrap and louvain for the moment), the third layer is a graining level
#'
#' @export

compute_alternative_SC <- function(
  counts,
  genes.use,
  meta.data = NULL,
  ge = NULL,
  N.comp = 10,
  k.seq = sort(c(5, 10, 50, 100, round(0.01*ncol(counts)), round(0.05*ncol(counts)))),
  res.seq = c(1, 2, 5, 10, 20, 30, 40, 50, 60, 70, 80),
  gamma.seq = NULL,
  return.singlecell.NW = FALSE,
  verbose = FALSE,
  common.gammas = TRUE,
  do.directed = c(T, F),
  group.singletons = FALSE,
  max.gamma = 150,
  ...
){

  N.c <- ncol(counts)
  cell.ids  <- colnames(counts)
  if(length(N.comp) == 1) N.comp <- 1:N.comp

  m.seurat <- Seurat::CreateSeuratObject(counts, meta.data = meta.data, ...)
  m.seurat <- Seurat::NormalizeData(m.seurat, ...)
  m.seurat <- Seurat::FindVariableFeatures(m.seurat, ...)
  m.seurat@assays$RNA@var.features <- genes.use

  m.seurat <- Seurat::ScaleData(m.seurat, features =  genes.use, ...)
  m.seurat <- Seurat::RunPCA(m.seurat, npcs = max(N.comp), ...)
  m.seurat <- Seurat::FindNeighbors(m.seurat, dims = N.comp, ...)

  # retreive seurat networks
  nw.seurat.list  <- list(
    'knn_seurat' = m.seurat@graphs$RNA_nn,
    'snn_seurat' = m.seurat@graphs$RNA_snn
  )

  nw.igraph.list <- list(
    'knn_seurat' = Graph_to_igraph(m.seurat@graphs$RNA_nn),
    'snn_seurat' = Graph_to_igraph(m.seurat@graphs$RNA_snn)
  )

  if(!is.null(ge)){
    X.for.pca            <- Matrix::t(ge[genes.use, ])
    X.for.pca            <- scale(X.for.pca)
    X.for.pca[is.na(X.for.pca)] <- 0
    sc.PCA               <- irlba::irlba(X.for.pca, nv = max(N.comp, 25))
    sc.PCA$x             <- sc.PCA$u %*% diag(sc.PCA$d)
    sc.PCA               <- sc.PCA$x[, N.comp]
  } else {
    sc.PCA <- Seurat::Embeddings(m.seurat)[,N.comp]
  }

  for(k in k.seq){
    for(directed in do.directed){
      graph.name <- paste0("knn_k_", k, "_")
      directed.ch <- ifelse(directed, "dir", "undir")
      graph.name <- paste0(graph.name, directed.ch)

      cur.nw <- SuperCell::build_knn_graph(
        X = sc.PCA,
        from = "coordinates",
        directed = directed,
        k = k
      )
      if(verbose) {
        print(paste("Graph supposed to be", directed.ch, ", and is.directed():", igraph::is.directed(cur.nw$graph.knn)))
      }

      nw.igraph.list[[graph.name]]    <- cur.nw$graph.knn
      adj.mtx                         <- igraph::get.adjacency(cur.nw$graph.knn)
      rownames(adj.mtx)               <- cell.ids
      colnames(adj.mtx)               <- cell.ids
      nw.seurat.list[[graph.name]]    <- Seurat::as.Graph(adj.mtx)
    }
  }

  # now graphs exist in 2 classes: seurat graphs (`nw.seurat.list`) and igraph graphs (`nw.igraph.list`)
  for(gname in names(nw.seurat.list)){
    m.seurat@graphs[[gname]] <- nw.seurat.list[[gname]]
  }

  res.list <- list()
  # from computed networks, compute super-cells with either walktrap or louvain
  for(gname in names(nw.seurat.list)){
    res.list[[gname]] <- list()

    # igraph graph
    sc.nw              <- nw.igraph.list[[gname]]

    # scale resolution: for larger k, number of clusters correponding to a given resolution increases
    k <- min(igraph::degree(sc.nw))

    ## louvain
    res.list[[gname]][["louvain"]] <- list()
    for(resolution in sort(res.seq, decreasing = TRUE)){
      resolution.normalized <- round(resolution/(sqrt(sqrt(k))), digits = 3)

      if(verbose) print(paste(gname, "louvain, res =",  resolution.normalized))

      m.seurat         <- FindClusters(m.seurat, resolution = resolution.normalized, graph.name = gname, group.singletons = group.singletons)

      cur.membership   <- as.numeric(as.vector(m.seurat$seurat_clusters)) + 1

      if(!group.singletons){
        singleton.idx  <- which(is.na(cur.membership))

        N.first        <- max(cur.membership, na.rm = T) + 1
        cur.membership[singleton.idx] <- N.first:(N.first+length(singleton.idx))
      }

      names(cur.membership) <- cell.ids

      res <- SC_from_membership(
        sc.nw = sc.nw,
        membership = cur.membership,
        gname = paste0(gname, "_louvain"),
        res = resolution.normalized,
        N.comp = N.comp,
        genes.use = genes.use,
        return.singlecell.NW = return.singlecell.NW
      )

      actual.gamma.ch <- as.character(res$gamma)
      if(verbose) print(paste(gname, "louvain, res =",  resolution.normalized, ", graining level =", actual.gamma.ch))
      res.list[[gname]][["louvain"]][[actual.gamma.ch]] <- res

      if(as.numeric(actual.gamma.ch) > max.gamma){
        warning("Gamma:", actual.gamma.ch , ">", max.gamma, ", thus break resolution loop")
        break()

      }
    }


    ## walktrap
    #get louvain gammas
    if(is.null(gamma.seq)){
      cur.gamma.seq <- as.numeric(names(res.list[[gname]][["louvain"]]))
      warning(paste("Since `gamma.seq` is NULL, it was petreived from the louvain partition!
                    \nNew gamma.seq = ", paste(gamma.seq, collapse = ",")))

    } else {
      if(common.gammas){
        cur.gamma.seq <- sort(unique(c(as.numeric(names(res.list[[gname]][["louvain"]])), gamma.seq)))
      } else {
        cur.gamma.seq <- gamma.seq
      }
    }

    for(gamma in cur.gamma.seq){
      k   <- round(N.c/gamma)

      g.s              <- igraph::cluster_walktrap(sc.nw)
      cur.membership   <- igraph::cut_at(g.s, k)
      names(cur.membership) <- cell.ids

      res <- SC_from_membership(
        sc.nw = sc.nw,
        membership = cur.membership,
        gname = paste0(gname, "_walktrap"),
        res = NULL,
        N.comp = N.comp,
        genes.use = genes.use,
        return.singlecell.NW = return.singlecell.NW
      )

      actual.gamma.ch <- as.character(res$gamma)
      res.list[[gname]][["walktrap"]][[actual.gamma.ch]] <- res
    }

  }
  return(res.list)
}
