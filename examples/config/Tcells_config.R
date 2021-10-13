### config file for Zheng dataset 

proj.name   <- "Tcells"
data.folder <- file.path("data", proj.name)
fig.folder  <- file.path("plots", proj.name)


.seed.seq           <- c(12345, 111, 19, 42, 7)
.seed               <-.seed.seq[1]

.N.clusters.seq     <- 2:10
.N.clusters.seq.breaks <- unique(c(0, seq(2, max(.N.clusters.seq), 2), max(.N.clusters.seq)))
.gamma.seq          <- c(1, 2, 5, 10, 50, 100, 200, 500, 1000) #1,2
.N.var.genes        <- 500
.genes.use          <- NULL
.N.comp             <- 10

.approx.N           <- 1000

.k.knn              <- 5

.SELECTED_METHODS   <- c("Exact", "Random", "Subsampling", "Metacell")


.logFC.thresh   <- 0.25
.pval.thresh    <- 0.05

.logFC.thresh.SC        <- 0.00

.do.extra.log.rescale  <- TRUE
.do.extra.sqrt.rescale <- FALSE

.do.frames <- TRUE

.min.SC.size.to.plot         <- c(0, 0, 0, 0, 1, 5, 10, 100, 200)+1 #c(0, 0, 0, 0, 1, 2, 3, 10, 20)+1
names(.min.SC.size.to.plot)  <- as.character(.gamma.seq)

.log.base                    <- c(2,2, 2,2, 1.9, 1.9, 1.6, 1.5, 1.5)
names(.log.base)  <- as.character(.gamma.seq)
