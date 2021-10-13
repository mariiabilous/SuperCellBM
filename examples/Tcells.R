##### "Super-cell analysis of Tcells (from Zheng et al., 2018)"
#
#
#
#
##### Copy Tcell.Rmd


if (!requireNamespace("igraph"))           install.packages("igraph")
if (!requireNamespace("RANN"))             install.packages("RANN")
if (!requireNamespace("WeightedCluster"))  install.packages("WeightedCluster")
if (!requireNamespace("corpcor"))          install.packages("corpcor")
if (!requireNamespace("weights"))          install.packages("weights")
if (!requireNamespace("Hmisc"))            install.packages("Hmisc")
if (!requireNamespace("Matrix"))           install.packages("Matrix")
if (!requireNamespace("patchwork"))        install.packages("patchwork")
if (!requireNamespace("plyr"))             install.packages("plyr")
if (!requireNamespace("irlba"))            install.packages("irlba")
if (!requireNamespace("scater"))           install.packages("scater")

if (!requireNamespace("remotes")) install.packages("remotes")
remotes::install_github("GfellerLab/SuperCell")
#remotes::install_github("tanaylab/metacell")

library(SuperCell)
#library(metacell)

devtools::install()
devtools::document()
library(SuperCellBM)

source("./config/Tcells_config.R")


###### ------------------------------ Flasgs ------------------------------ ######
ToLoadPreprocData <- T
ToComputeSC       <- F
ToComputeSC_GE    <- F
ToClusterSinglecell <- F
ToComputeAltClustering <- T
ToComputePCACLustering <- T

ToTestPackage     <- F
N.clusters        <- 3

filename_suf <- "" # variable to add a suffix to the saved files in case of testing of the package

if(ToTestPackage){
  testing_gamma_seq <- c(1, 10, 100)
  testing_seed_seq <- .seed.seq[1:3]

  warning(paste("The reduced set of graining leveles and seeds will be used, to get real output, turn it ti FALSE"))
  warning(paste("Original set of graining levels is:", paste(.gamma.seq, collapse = ", "),
                "but used testing set is:", paste(testing_gamma_seq, collapse = ", ")))
  warning(paste("Original set of seeds  is:", paste(.seed.seq, collapse = ", "),
                "but used testing set is:", paste(testing_seed_seq, collapse = ", ")))

  .gamma.seq <- testing_gamma_seq
  .seed.seq <- testing_seed_seq

  filename_suf = "_testing_package"
} else {
  .seed.seq = .seed.seq[1:5]
}

###### ------------------------------ Load data ------------------------------ ######


if(!ToLoadPreprocData){
  source("./config/Tcells_load_data.R")
} else {
  sc.GE      <- readRDS(file = file.path(data.folder, "sc_ge.Rds"))
  sc.counts  <- readRDS(file = file.path(data.folder, "sc_counts.Rds"))
  cell.meta  <- readRDS(file = file.path(data.folder, "cell_meta.Rds"))
}


###### ------------------------------ Get main values from the dataset ------------------------------ ######

GT.cell.type                <- cell.meta$GT.cell.type
names(GT.cell.type)         <- rownames(cell.meta)

GT.cell.type.names          <- names(table(cell.meta$GT.cell.type))
GT.cell.type.2.num          <- 1:length(unique(GT.cell.type))
names(GT.cell.type.2.num)   <- GT.cell.type.names
GT.cell.type.num            <- GT.cell.type.2.num[GT.cell.type]
names(GT.cell.type.num)     <- names(GT.cell.type)

cell.ids                    <- colnames(sc.GE)
N.c                         <- length(cell.ids)
gene.names                  <- rownames(sc.GE)
N.g                         <- length(gene.names)

N.clusters                  <- length(unique(GT.cell.type))

mito.genes                  <- grep(pattern = "^MT",     x = gene.names, value = TRUE)
ribo.genes                  <- grep(pattern = "^RP[LS]", x = gene.names, value = TRUE)
mito.ribo.genes             <- c(mito.genes, ribo.genes)
length(mito.ribo.genes)
.genes.omit                  <- mito.ribo.genes


## uncomment this when needed

# pal.GT                      <- color.tsne.
# scales::show_col(pal.GT)
# scales::show_col(color.tsne)

###### ------------------------------ Compute SuperCells ------------------------------ ######

filename <- paste0('initial', filename_suf)

SC.list <- compute_supercells(
  sc.GE,
  ToComputeSC = ToComputeSC,
  data.folder = data.folder,
  filename = filename,
  gamma.seq = .gamma.seq,
  n.var.genes = .N.var.genes,
  k.knn = .k.knn,
  n.pc = .N.comp,
  approx.N = .approx.N,
  fast.pca = TRUE,
  genes.use = .genes.use,
  genes.exclude = .genes.omit,
  seed.seq = .seed.seq,
  verbose = TRUE
)

cat(paste("Super-cell computed for:", paste(names(SC.list), collapse = ", "),
          "\nat graining levels:", paste(names(SC.list[['Approx']]), collapse = ", "),
          "\nfor seeds:", paste(names(SC.list[['Approx']][[1]]), collapse = ", "), "\n",
          "\nand saved to / loaded from", paste0(filename, ".Rds")))


###### ------------------------------ Compute MetaCells ------------------------------ ######

SC.mc <- compute_supercells_metacells(
  sc.counts = sc.counts,
  gamma.seq = .gamma.seq,
  SC.list = SC.list,
  proj.name = proj.name,
  ToComputeSC = ToComputeSC,
  mc.k.knn = 100,
  T_vm_def = 0.08,
  MC.folder = "MC",
  MC_gene_settings = c('Metacell_default', 'Metacell_SC_like') # do not change
)

###### ------------------------------ Get actual graining levels obtained with Metacell ------------------------------ ######

additional_gamma_seq <- get_actual_gammas_metacell(SC.mc)

cat(paste("Metacells were computed in", length(names(SC.mc)), "settings:", paste(names(SC.mc), collapse = ", "),
          "\nfor Gammas:", paste(names(SC.mc[[1]]), collapse = ", "),
          "\nbut actual gammas are:", paste(additional_gamma_seq, collapse = ", ")
))


# manually expand MC because later we will have 2 different setting for MC profile: fp - footpring of MC, av - averaged
SC.mc.fp <- SC.mc
names(SC.mc.fp) <- sapply(names(SC.mc), FUN = function(x){paste0(x, '_fp')})

SC.mc.av <- SC.mc
names(SC.mc.av) <- sapply(names(SC.mc), FUN = function(x){paste0(x, '_av')})

SC.mc.expanded <- c(SC.mc.fp, SC.mc.av)

names(SC.mc.expanded)
rm(SC.mc.fp, SC.mc.av, SC.mc)


###### ------------------------------ Compute SC at additional gammas ------------------------------ ######
filename <- paste0('additional_gammas', filename_suf)

SC.list <- compute_supercells_additional_gammas(
  SC.list,
  additional_gamma_seq = additional_gamma_seq,
  ToComputeSC = ToComputeSC,
  data.folder = data.folder,
  filename = filename,
  approx.N = .approx.N,
  fast.pca = TRUE
)

cat(paste("Super-cells of methods:", paste(names(SC.list), collapse = ", "),
          "\nwere computed at aggitional graining levels:", paste(additional_gamma_seq, collapse = ", "),
          "\nand added to SC.list"
))

###### ------------------------------ Concat Supercell and Metacell ------------------------------ ######
SC.list <- c(SC.list, SC.mc.expanded)
rm(SC.mc.expanded)

filename <- paste0("all", filename_suf)
saveRDS(SC.list, file = file.path(data.folder, "SC", paste0(filename, ".Rds")))

cat(paste(
  "Metacell data added to SC.list \nand now it contains:",
  paste(names(SC.list), collapse = ", "),
  "\nSC.list was saved to", file.path(data.folder, "SC", paste0(filename, ".Rds"))
))

###### ------------------------------ Cluster single-cell data ------------------------------ ######
genes.use <- SC.list$Exact$`10`$`12345`$genes.use
if(length(.N.comp) == 1) {N.comp <- 1:.N.comp} else {N.comp <- .N.comp}

X.for.pca       <- scale(Matrix::t(sc.GE[genes.use,]))
sc.pca          <- irlba::irlba(X.for.pca, nv = max(N.comp, 25))
sc.pca$x        <- sc.pca$u %*% diag(sc.pca$d)

if(ToClusterSinglecell){
  sc.dist <- dist(sc.pca$x[,N.comp])
  sc.clustering.hcl <- hclust(sc.dist, method = "ward.D2")

  saveRDS(sc.clustering.hcl, file = file.path(data.folder, "sc_clustering_hcl.Rds"))
} else {
  sc.clustering.hcl <- readRDS(file.path(data.folder, "sc_clustering_hcl.Rds"))
}

sc.clustering <- cutree(sc.clustering.hcl, k = N.clusters)

###### ------------------------------ Annotate super-cells based on single-cell clustering ------------------------------ ######

SC.list <- annotate_supercells_to_cluster(
  SC.list = SC.list,
  sc.annotation = sc.clustering,
  annotation.name = 'sc_clustering'
)

purity.df <- plot_annotation_purity(
  SC.list,
  annotation.name = 'sc_clustering', w = 2.5, h = 2.5,
)


###### ------------------------------ Compute GE ------------------------------ ######

filename <- paste0("all", filename_suf)

SC.GE.list <- compute_supercells_GE(
  sc.GE = sc.GE,
  SC.list = SC.list,
  ToComputeSC_GE = ToComputeSC_GE,
  data.folder = data.folder,
  filename = filename,
  verbose = FALSE
)

cat(paste("Gene expression profile computed for:", paste(names(SC.GE.list), collapse = ", "),
          "\nat graining levels:", paste(sort(as.numeric(names(SC.GE.list[['Approx']]))), collapse = ", "),
          "\nfor seeds:", paste(names(SC.GE.list[['Approx']][[1]]), collapse = ", "),
          "\nand saved to / loaded from", paste0(filename, ".Rds")
))

###### ------------------------------ PCA + clustering of super-cells ------------------------------ ######
filename <- paste0("all_PCA_clust", filename_suf)

if(ToComputePCACLustering){
  SC.list <- compute_supercells_PCA(
    SC.list = SC.list,
    SC.GE.list = SC.GE.list,
    N.comp = .N.comp
  )

  SC.list <- compute_supercells_clustering(
    SC.list,
    N.comp = 10,
    N.clusters.seq = c(2:10),
    pca_name = 'SC_PCA'
  )

  saveRDS(SC.list, file = file.path(data.folder, 'SC', paste0(filename, '.Rds')))
} else {
  SC.list <- readRDS(file = file.path(data.folder, 'SC', paste0(filename, '.Rds')))
}

###### ------------------------------ Alternative clustering of single-cell data ------------------------------ ######

if(ToComputeAltClustering){
  alt.clusterings <- compute_alternative_clustering(
    sc.pca = sc.pca$x,
    N.comp = .N.comp,
    N.clusters.seq = .N.clusters.seq,
    hclust_methods = c("average", "mcquitty", "median"),
    seed.seq = .seed.seq,
    verbose = T
  )
  saveRDS(sc.clustering, file = file.path(data.folder, "sc_alternative_clustering.Rds"))
} else {
  alt.clusterings <- readRDS(file.path(data.folder, "sc_alternative_clustering.Rds"))
}

###### ------------------------------ Consistency between sc and SC clustering ------------------------------ ######

clust.consistency <- compute_consistency_of_supercell_clustering(
  SC.list = SC.list,
  sc.annotation = sc.clustering,
  sc.clustering = sc.clustering,
  sc.alternative.clustering = alt.clusterings
)

clust.consistency_summary <- plot_clustering_consistency(
  clust.consistency,
  min.value.alt.clustering = 0.,
  error_bars = "extr"
)

