### config file for TIAN dataset

proj.name   <- "Tian"
data.folder <- file.path("examples", "data", proj.name)
fig.folder  <- file.path("examples", "plots", proj.name)


.seed.seq           <- c(12345, 111, 19, 42, 7, 559241, 123, 987, 234, 91, 877, 451, 817)
.seed               <- .seed.seq[1]

.N.clusters.seq     <- 2:10
.N.clusters.seq.breaks <- unique(c(0, seq(2, max(.N.clusters.seq), 2), max(.N.clusters.seq)))
.gamma.seq          <- c(1, 2, 5, 10, 20, 50, 100, 200)
.N.var.genes        <- 1000
.genes.use          <- NULL
.genes.omit         <- NULL
.N.comp             <- 10

.approx.N           <- 1000

.k.knn              <- 5

.SELECTED_METHODS   <- c("Exact", "Random", "Subsampling", "Metacell")


.logFC.thresh   <- 0.5
.pval.thresh    <- 0.05

.logFC.thresh.SC        <- 0.25

.do.extra.log.rescale  <- FALSE
.do.extra.sqrt.rescale <- FALSE

.do.frames <- TRUE


.min.SC.size.to.plot         <- c(0, 0, 0, 0, 0, 0, 0, 0)+1
names(.min.SC.size.to.plot)  <- as.character(.gamma.seq)
