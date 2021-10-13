cell.type.name <- c("CD4_naive", "CD8_naive", "CD8_cytotoxic", "CD4_helper")

print(cell.type.name)
GE.list <- list()

for(ct.name in cell.type.name){
  cur.data.path <- file.path(data.folder, ct.name)

  cur.ge        <- Matrix::readMM(file.path(cur.data.path, "matrix.mtx"))

  cur.barcodes  <- as.vector(read.table(file = file.path(cur.data.path, "barcodes.tsv"), sep = '\t', header = FALSE)$V1)
  cur.barcodes  <- paste0(ct.name, "_", cur.barcodes)

  cur.genes     <- as.vector(read.table(file = file.path(cur.data.path, "genes.tsv"), sep = '\t', header = FALSE)$V2)

  cur.ge@Dimnames[[1]] <- cur.genes
  cur.ge@Dimnames[[2]] <- cur.barcodes

  GE.list[[ct.name]] <- cur.ge

  print(ct.name)
}

mutual.genes <- GE.list[[1]]@Dimnames[[1]]
for(ct.name in cell.type.name){
  mutual.genes <- intersect(mutual.genes, GE.list[[ct.name]]@Dimnames[[1]])
}

sc.counts              <- c()

GT.cell.type.fine <- c()
GT.cell.type      <- c()
for(ct.name in cell.type.name){
  cur.counts <- GE.list[[ct.name]]
  sc.counts  <- cbind(sc.counts, cur.counts[mutual.genes,])

  GT.cell.type.fine <- c(GT.cell.type.fine, rep(ct.name, GE.list[[ct.name]]@Dim[[2]]))
  GT.cell.type      <- c(GT.cell.type, rep(gsub("_.*", "", ct.name), GE.list[[ct.name]]@Dim[[2]]))
}


names(GT.cell.type.fine) <- colnames(sc.counts)
names(GT.cell.type)      <- colnames(sc.counts)

n.UMI <- colSums(sc.counts)
hist(n.UMI)
summary(n.UMI)
remove.cells.UMI <- names(n.UMI)[n.UMI > 3000  | n.UMI < 800]


gene.detection              <- rowSums(sc.counts > 0)/ncol(sc.counts)
genes.to.remove             <- names(gene.detection)[gene.detection < 0.01 | gene.detection > 0.9]
genes.to.keep               <- setdiff(rownames(sc.counts), genes.to.remove)

gene.names                  <- genes.to.keep
sc.counts                   <- sc.counts[genes.to.keep,]


sc.GE    <- scater::normalizeCounts(sc.counts)

remove.cells <- colnames(sc.GE)[(sc.GE["CD4",] > 0.01 & GT.cell.type == "CD8") |
                                  ((sc.GE["CD8A",] > 0.01 | sc.GE["CD8B",] > 0.01) & GT.cell.type == "CD4")]

remove.cells <- unique(c(remove.cells, remove.cells.UMI))


keep.cells         <- setdiff(colnames(sc.counts), remove.cells)
print(paste("Keep:", length(keep.cells), "cells"))

sc.counts          <- sc.counts[, keep.cells]
sc.GE              <- sc.GE[, keep.cells]
GT.cell.type.fine  <- GT.cell.type.fine[keep.cells]
GT.cell.type       <- GT.cell.type[keep.cells]

cell.meta <- data.frame(GT.cell.type.fine = GT.cell.type.fine, GT.cell.type = GT.cell.type, cell.ids = keep.cells)

saveRDS(sc.GE, file = file.path("data", proj.name, "sc_ge.Rds"))
saveRDS(sc.counts, file = file.path("data", proj.name, "sc_counts.Rds"))
saveRDS(cell.meta, file = file.path("data", proj.name, "cell_meta.Rds"))


rm(GE.list, cur.ge, cur.counts, cur.barcodes, cur.data.path, cur.genes,
   gene.detection, genes.to.keep, genes.to.remove, keep.cells,
   GT.cell.type, GT.cell.type.fine)
