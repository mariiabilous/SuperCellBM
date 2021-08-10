An example how to build multiple super-cell-like objects, including
‘exact’ (Super-cells obtained with the exact coarse-gaining), ‘approx’
(Super-cells obtained with the aaproximate coarse-graining),
’metacell(.\*)’ (Metacell build on the same genes as super-cells –
‘metacell\_SC\_like’; and Metacell build in a default set of genes –
‘metacell\_default’) for a set of graining levels and random seeds.

``` r
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("SingleCellExperiment")
# 
# if (!requireNamespace("remotes")) install.packages("remotes")
# remotes::install_github("GfellerLab/SuperCell")
# remotes::install_github("mariiabilous/SuperCellBM")

library(SingleCellExperiment)
library(SuperCell)
library(SuperCellBM)
```

Load some default parameters
----------------------------

Such as `.gamma.seq`for the set of fraining levels, `.seed.seq` for the
set of random seeds, `adata.folder` and `fig.folder` for the folders
where to write data and plots

``` r
source("./examples/config/Tian_config.R")
```

### Flags

whether to compute super-cell (`ToComputeSC`) or whether to compute
super-cell gene expression (`ToComputeSC_GE`) or load saved files. Make
sure, this file exists :)

``` r
ToComputeSC <- T
ToComputeSC_GE <- T
```

\#\#Load `cell_lines` data from
(<a href="https://pubmed.ncbi.nlm.nih.gov/31133762/" class="uri">https://pubmed.ncbi.nlm.nih.gov/31133762/</a>).

``` r
RData.file.path <- file.path(data.folder, 'cell_lines_git.RData')

if(!file.exists(RData.file.path)){
  if(!dir.exists(data.folder)) dir.create(data.folder, recursive = T)
  download.file('https://github.com/LuyiTian/sc_mixology/blob/master/data/sincell_with_class_5cl.RData?raw=true', 
                RData.file.path)
}

load(RData.file.path)

# keep used dataset 
cell_lines_SCE <- sce_sc_10x_5cl_qc

#remove not-used datasets
rm(sc_Celseq2_5cl_p1, sc_Celseq2_5cl_p2, sc_Celseq2_5cl_p3, sce_sc_10x_5cl_qc)
```

Get and set the main variables, such as single-cell gene expression
(`sc.GE`), single-cell counts (`sc.counts`), number of single cells
(`N.c`) and total number of genes (`N.g`). Set matrix column names to
cellIDs (`cell.ids`) and row names to gene names (`gene.names`)

``` r
cell.ids     <- cell_lines_SCE@colData@rownames
N.c          <- cell_lines_SCE@colData@nrows 


gene.names   <- cell_lines_SCE@int_elementMetadata$external_gene_name
N.g          <- length(gene.names)

sc.GE           <- cell_lines_SCE@assays$data$logcounts
colnames(sc.GE) <- cell.ids
rownames(sc.GE) <- gene.names

sc.counts   <- cell_lines_SCE@assays$data$counts
colnames(sc.counts) <- cell.ids
rownames(sc.counts) <- gene.names

## this is not needed at this point, but will be used later
GT.cell.type <- cell_lines_SCE@colData$cell_line_demuxlet
names(GT.cell.type) <- cell.ids
N.clusters         <- length(unique(GT.cell.type))

GT.cell.type.names          <- names(table(GT.cell.type))
GT.cell.type.2.num          <- 1:length(unique(GT.cell.type))
names(GT.cell.type.2.num)   <- GT.cell.type.names
GT.cell.type.num            <- GT.cell.type.2.num[GT.cell.type]
names(GT.cell.type.num)     <- names(GT.cell.type)

## uncomment this when needed 
#.pal.GT <- .color.tsne.Tian  ## to Global Config
#scales::show_col(.pal.GT)


#mito.genes <-  grep(pattern = "^MT", x = gene.names, value = TRUE)
#ribo.genes <-  grep(pattern = "^RP[LS]", x = gene.names, value = TRUE)
#mito.ribo.genes <- c(mito.genes, ribo.genes)
#length(mito.ribo.genes)

#gene.meta <- data.frame(name = gene.names, inNcells = rowSums(sc.GE>0), mean.expr = rowMeans(sc.GE), sd = rowSds(sc.GE))
#head(gene.meta)
```

Compute Super-cell structure
----------------------------

for the Exact, Aprox (Super-cells obtained with the exact or approximate
coarse-graining), Subsampling or Random (random grouping of cells into
super-cells)

``` r
SC.list <- compute_supercells(
  sc.GE,
  ToComputeSC = ToComputeSC,
  data.folder = data.folder,
  filename = 'initial',
  gamma.seq = .gamma.seq,
  n.var.genes = .N.var.genes,
  k.knn = .k.knn,
  n.pc = .N.comp,
  approx.N = 1000,
  fast.pca = TRUE,
  genes.use = .genes.use, 
  genes.exclude = .genes.omit,
  seed.seq = .seed.seq
  )
```

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is set to N.SC 3918

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is set to N.SC 3918

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is set to N.SC 3918

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is set to N.SC 3918

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is set to N.SC 3918

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is set to N.SC 3918

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is set to N.SC 3918

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is set to N.SC 3918

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is set to N.SC 3918

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is set to N.SC 3918

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is set to N.SC 3918

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is set to N.SC 3918

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is set to N.SC 3918

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is set to N.SC 1959

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is set to N.SC 1959

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is set to N.SC 1959

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is set to N.SC 1959

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is set to N.SC 1959

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is set to N.SC 1959

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is set to N.SC 1959

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is set to N.SC 1959

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is set to N.SC 1959

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is set to N.SC 1959

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is set to N.SC 1959

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is set to N.SC 1959

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is set to N.SC 1959

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

    ## Warning in SuperCell::SCimplify(X = sc.GE, genes.use = genes.use, genes.exclude
    ## = genes.exclude, : N.approx is not much larger than desired number of super-
    ## cells, so an approximate simplification may take londer than an exact one!

``` r
cat(paste("Super-cell computed for:", paste(names(SC.list), collapse = ", "), 
          "\nat graining levels:", paste(names(SC.list[['approx']]), collapse = ", "),
          "\nfor seeds:", paste(names(SC.list[['approx']][[1]]), collapse = ", "), "\n"))
```

    ## Super-cell computed for: Exact, Approx, Random, Subsampling 
    ## at graining levels:  
    ## for seeds:

Compute metacells in two settings:
----------------------------------

-   ‘metacell\_default’ - when metacell is computed with the default
    parameters (from the tutorial), using gene set filtered by MC
-   ‘metacell\_SC\_like’ - when metacell is computed with the same
    default parameters, but at the same set of genes as Super-cells

``` r
SC.mc <- compute_supercells_metacells(
  sc.counts = sc.counts, 
  gamma.seq = .gamma.seq,
  SC.list = SC.list,
  proj.name = proj.name,
  ToComputeSC = ToComputeSC, 
  mc.k.knn = 100,
  T_vm_def = 0.08,
  MC.folder = "MC", 
  MC_gene_settings = c('metacell_default', 'metacell_SC_like') # do not change
)
```

    ## initializing scdb to examples/data/Tian/MC/metacell_default

    ## Calculating gene statistics...

    ## will downsamp

    ## done downsamp

    ## will gen mat_n

    ## done gen mat_n

    ## done computing basic gstat, will compute trends

    ## ..done

    ## will build balanced knn graph on 3918 cells and 2940 genes, this can be a bit heavy for >20,000 cells

    ## [1] "done i,j,x"
    ## [1] "done W as sparse matrix"
    ## [1] "done colnames, rownames of W"
    ## [1] "done remove outliers"
    ## [1] "done graph form adj matrix (sparse)"
    ## [1] "done as.undirected graph"

    ## running bootstrap to generate cocluster

    ## done resampling

    ## filtered 379554 left with 153208 based on co-cluster imbalance

    ## building metacell object, #mc 87

    ## add batch counts

    ## compute footprints

    ## compute absolute ps

    ## compute coverage ps

    ## reordering metacells by hclust and most variable two markers

    ## reorder on IFI27 vs KRT81

    ## [1] "min_mc_size 1"

    ## comp mc graph using the graph Tian_cgraph_from_bknn and K 20

    ## running bootstrap to generate cocluster

    ## done resampling

    ## filtered 379554 left with 153208 based on co-cluster imbalance

    ## building metacell object, #mc 87

    ## add batch counts

    ## compute footprints

    ## compute absolute ps

    ## compute coverage ps

    ## reordering metacells by hclust and most variable two markers

    ## reorder on IFI27 vs KRT81

    ## [1] "min_mc_size 1"

    ## comp mc graph using the graph Tian_cgraph_from_bknn and K 20

    ## running bootstrap to generate cocluster

    ## done resampling

    ## filtered 379554 left with 153208 based on co-cluster imbalance

    ## building metacell object, #mc 87

    ## add batch counts

    ## compute footprints

    ## compute absolute ps

    ## compute coverage ps

    ## reordering metacells by hclust and most variable two markers

    ## reorder on IFI27 vs KRT81

    ## [1] "min_mc_size 1"

    ## comp mc graph using the graph Tian_cgraph_from_bknn and K 20

    ## running bootstrap to generate cocluster

    ## done resampling

    ## filtered 420498 left with 164711 based on co-cluster imbalance

    ## building metacell object, #mc 74

    ## add batch counts

    ## compute footprints

    ## compute absolute ps

    ## compute coverage ps

    ## reordering metacells by hclust and most variable two markers

    ## reorder on IFI27 vs AKR1B10

    ## [1] "min_mc_size 3"

    ## comp mc graph using the graph Tian_cgraph_from_bknn and K 20

    ## running bootstrap to generate cocluster

    ## done resampling

    ## filtered 432803 left with 170010 based on co-cluster imbalance

    ## building metacell object, #mc 66

    ## add batch counts

    ## compute footprints

    ## compute absolute ps

    ## compute coverage ps

    ## reordering metacells by hclust and most variable two markers

    ## reorder on IFI27 vs AKR1B10

    ## [1] "min_mc_size 4"

    ## comp mc graph using the graph Tian_cgraph_from_bknn and K 20

    ## running bootstrap to generate cocluster

    ## done resampling

    ## filtered 432803 left with 170010 based on co-cluster imbalance

    ## building metacell object, #mc 66

    ## add batch counts

    ## compute footprints

    ## compute absolute ps

    ## compute coverage ps

    ## reordering metacells by hclust and most variable two markers

    ## reorder on IFI27 vs AKR1B10

    ## [1] "min_mc_size 4"

    ## comp mc graph using the graph Tian_cgraph_from_bknn and K 20

    ## running bootstrap to generate cocluster

    ## done resampling

    ## filtered 473321 left with 187723 based on co-cluster imbalance

    ## building metacell object, #mc 58

    ## add batch counts

    ## compute footprints

    ## compute absolute ps

    ## compute coverage ps

    ## reordering metacells by hclust and most variable two markers

    ## reorder on IFI27 vs KRT81

    ## [1] "min_mc_size 8"

    ## comp mc graph using the graph Tian_cgraph_from_bknn and K 20

    ## running bootstrap to generate cocluster

    ## done resampling

    ## filtered 473321 left with 187723 based on co-cluster imbalance

    ## building metacell object, #mc 58

    ## add batch counts

    ## compute footprints

    ## compute absolute ps

    ## compute coverage ps

    ## reordering metacells by hclust and most variable two markers

    ## reorder on IFI27 vs KRT81

    ## [1] "min_mc_size 8"

    ## comp mc graph using the graph Tian_cgraph_from_bknn and K 20

    ## initializing scdb to examples/data/Tian/MC/metacell_SC_like

    ## Calculating gene statistics...

    ## will downsamp

    ## done downsamp

    ## will gen mat_n

    ## done gen mat_n

    ## done computing basic gstat, will compute trends

    ## ..done

    ## will build balanced knn graph on 3918 cells and 1000 genes, this can be a bit heavy for >20,000 cells

    ## [1] "done i,j,x"
    ## [1] "done W as sparse matrix"
    ## [1] "done colnames, rownames of W"
    ## [1] "done remove outliers"
    ## [1] "done graph form adj matrix (sparse)"
    ## [1] "done as.undirected graph"

    ## running bootstrap to generate cocluster

    ## done resampling

    ## filtered 422867 left with 160291 based on co-cluster imbalance

    ## building metacell object, #mc 85

    ## add batch counts

    ## compute footprints

    ## compute absolute ps

    ## compute coverage ps

    ## reordering metacells by hclust and most variable two markers

    ## reorder on GPX2 vs IFITM3

    ## [1] "min_mc_size 1"

    ## comp mc graph using the graph Tian_cgraph_from_bknn and K 20

    ## running bootstrap to generate cocluster

    ## done resampling

    ## filtered 422867 left with 160291 based on co-cluster imbalance

    ## building metacell object, #mc 85

    ## add batch counts

    ## compute footprints

    ## compute absolute ps

    ## compute coverage ps

    ## reordering metacells by hclust and most variable two markers

    ## reorder on GPX2 vs IFITM3

    ## [1] "min_mc_size 1"

    ## comp mc graph using the graph Tian_cgraph_from_bknn and K 20

    ## running bootstrap to generate cocluster

    ## done resampling

    ## filtered 422867 left with 160291 based on co-cluster imbalance

    ## building metacell object, #mc 85

    ## add batch counts

    ## compute footprints

    ## compute absolute ps

    ## compute coverage ps

    ## reordering metacells by hclust and most variable two markers

    ## reorder on GPX2 vs IFITM3

    ## [1] "min_mc_size 1"

    ## comp mc graph using the graph Tian_cgraph_from_bknn and K 20

    ## running bootstrap to generate cocluster

    ## done resampling

    ## filtered 469880 left with 174233 based on co-cluster imbalance

    ## building metacell object, #mc 72

    ## add batch counts

    ## compute footprints

    ## compute absolute ps

    ## compute coverage ps

    ## reordering metacells by hclust and most variable two markers

    ## reorder on IFI27 vs CYP24A1

    ## [1] "min_mc_size 3"

    ## comp mc graph using the graph Tian_cgraph_from_bknn and K 20

    ## running bootstrap to generate cocluster

    ## done resampling

    ## filtered 488185 left with 180369 based on co-cluster imbalance

    ## building metacell object, #mc 64

    ## add batch counts

    ## compute footprints

    ## compute absolute ps

    ## compute coverage ps

    ## reordering metacells by hclust and most variable two markers

    ## reorder on IFI27 vs ALDH1A1

    ## [1] "min_mc_size 4"

    ## comp mc graph using the graph Tian_cgraph_from_bknn and K 20

    ## running bootstrap to generate cocluster

    ## done resampling

    ## filtered 488185 left with 180369 based on co-cluster imbalance

    ## building metacell object, #mc 64

    ## add batch counts

    ## compute footprints

    ## compute absolute ps

    ## compute coverage ps

    ## reordering metacells by hclust and most variable two markers

    ## reorder on IFI27 vs ALDH1A1

    ## [1] "min_mc_size 4"

    ## comp mc graph using the graph Tian_cgraph_from_bknn and K 20

    ## running bootstrap to generate cocluster

    ## done resampling

    ## filtered 540717 left with 201103 based on co-cluster imbalance

    ## building metacell object, #mc 57

    ## add batch counts

    ## compute footprints

    ## compute absolute ps

    ## compute coverage ps

    ## reordering metacells by hclust and most variable two markers

    ## reorder on IFI27 vs KRT81

    ## [1] "min_mc_size 8"

    ## comp mc graph using the graph Tian_cgraph_from_bknn and K 20

    ## running bootstrap to generate cocluster

    ## done resampling

    ## filtered 540717 left with 201103 based on co-cluster imbalance

    ## building metacell object, #mc 57

    ## add batch counts

    ## compute footprints

    ## compute absolute ps

    ## compute coverage ps

    ## reordering metacells by hclust and most variable two markers

    ## reorder on IFI27 vs KRT81

    ## [1] "min_mc_size 8"

    ## comp mc graph using the graph Tian_cgraph_from_bknn and K 20

``` r
additional_gamma_seq <- get_actual_gammas_metacell(SC.mc)

cat(paste("Metacells were computed in", length(names(SC.mc)), "settings:", paste(names(SC.mc), collapse = ", "), "\n",
          "For Gammas:", paste(names(SC.mc[[1]]), collapse = ", "), "\n",
          "But actual gammas are:", paste(additional_gamma_seq, collapse = ", "), "\n"
))
```

    ## Metacells were computed in 2 settings: metacell_default, metacell_SC_like 
    ##  For Gammas: 1, 2, 5, 10, 20, 50, 100, 200 
    ##  But actual gammas are: 46, 54, 61, 69

``` r
# manually expand MC because later we will have 2 diferetn setting for MC profile: fp - footpring of MC, av - averaged 
SC.mc.fp <- SC.mc
names(SC.mc.fp) <- sapply(names(SC.mc), FUN = function(x){paste0(x, '_fp')})

SC.mc.av <- SC.mc
names(SC.mc.av) <- sapply(names(SC.mc), FUN = function(x){paste0(x, '_av')})

SC.mc.expanded <- c(SC.mc.fp, SC.mc.av)

names(SC.mc.expanded)
```

    ## [1] "metacell_default_fp" "metacell_SC_like_fp" "metacell_default_av"
    ## [4] "metacell_SC_like_av"

``` r
rm(SC.mc.fp, SC.mc.av, SC.mc)
```

Compute super-cells (Exact, Approx, Subsampling and Random) at the addidional grainig levels obtained with Metacell.
--------------------------------------------------------------------------------------------------------------------

So that we can dirrectly compare the results of Super-cells and
Metacells at the same graining levels.

``` r
SC.list <- compute_supercells_additional_gammas(
  SC.list,
  additional_gamma_seq = additional_gamma_seq,
  ToComputeSC = TRUE,
  data.folder = data.folder,
  filename = 'additional_gammas',
  approx.N = 1000,
  fast.pca = TRUE
)
```

    ## [1] "GAMMMA: 46"
    ## [1] "Exact"
    ## [1] "SEED: 12345"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 111"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 19"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 42"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 7"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 559241"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 123"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 987"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 234"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 91"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 877"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 451"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 817"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "GAMMMA: 54"
    ## [1] "Exact"
    ## [1] "SEED: 12345"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 111"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 19"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 42"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 7"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 559241"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 123"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 987"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 234"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 91"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 877"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 451"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 817"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "GAMMMA: 61"
    ## [1] "Exact"
    ## [1] "SEED: 12345"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 111"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 19"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 42"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 7"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 559241"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 123"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 987"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 234"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 91"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 877"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 451"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 817"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "GAMMMA: 69"
    ## [1] "Exact"
    ## [1] "SEED: 12345"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 111"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 19"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 42"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 7"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 559241"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 123"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 987"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 234"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 91"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 877"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 451"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"
    ## [1] "SEED: 817"
    ## [1] "Approx"
    ## [1] "Random"
    ## [1] "Subsampling"

``` r
print(paste("Super-cells of methods:", paste(names(SC.list), collapse = ", "), 
      "were computed at aggitional graining levels:", paste(additional_gamma_seq, collapse = ", "), "
      and added to SC.list"
      ))
```

    ## [1] "Super-cells of methods: Exact, Approx, Random, Subsampling were computed at aggitional graining levels: 46, 54, 61, 69 \n      and added to SC.list"

### Concatenate Metacells to the list of Super-cells

``` r
SC.list <- c(SC.list, SC.mc.expanded)
names(SC.list)
```

    ## [1] "Exact"               "Approx"              "Random"             
    ## [4] "Subsampling"         "metacell_default_fp" "metacell_SC_like_fp"
    ## [7] "metacell_default_av" "metacell_SC_like_av"

``` r
rm(SC.mc.expanded)
```

Compute GE for Super-cell data
------------------------------

GE profile for the super-cell data is computede: - for super-cells
(Exact, Approx) by averaging gene expression within super-cells, - for
Random, also averaging gene expression within super-cells - for the
Subsampling, sc.GE matrix is just subsampled, - for Metacells, gene
expression is computed in 2 ways: 1) the same as super-cells (averaging
gene expression within Metacells) -&gt; `metacell_(.*)_av` 2) using the
default output of Metacell maned footprint -&gt; `metacell_(.*)_fp`

``` r
SC.GE.list <- compute_supercells_GE(
  sc.GE = sc.GE, 
  SC.list = SC.list,
  ToComputeSC_GE = ToComputeSC_GE, 
  data.folder = data.folder,
  filename = 'all'
)

cat(paste("Gene expression profile computed for:", paste(names(SC.GE.list), collapse = ", "), 
    "\nat graining levels:", paste(sort(as.numeric(names(SC.GE.list[['Approx']]))), collapse = ", "),
    "\nfor seeds:", paste(names(SC.GE.list[['Approx']][[1]]), collapse = ", ")))
```

    ## Gene expression profile computed for: Exact, Approx, Random, Subsampling, metacell_default_fp, metacell_SC_like_fp, metacell_default_av, metacell_SC_like_av 
    ## at graining levels: 1, 2, 5, 10, 20, 46, 50, 54, 61, 69, 100, 200 
    ## for seeds: 12345, 111, 19, 42, 7, 559241, 123, 987, 234, 91, 877, 451, 817
