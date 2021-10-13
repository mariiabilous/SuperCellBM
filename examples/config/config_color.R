##### Colors for plots in the manuscript

.color.simplification          <- "darkred"
.color.approx.simplification   <- "royalblue"
.color.subsampling             <- "black"
.color.random.grouping         <- "gray"
.color.metacell                <- "salmon"

.pal.simplification.3          <- c(.color.simplification, .color.subsampling, .color.random.grouping)
names(.pal.simplification.3)   <- c("Exact", "Subsampling", "Random")
.shape.simplification.3        <- c(1, 2, 3)
names(.shape.simplification.3) <- names(.pal.simplification.3)

.pal.simplification.4          <- c(.color.simplification, .color.approx.simplification, .color.subsampling, .color.random.grouping)
names(.pal.simplification.4)   <- c("Exact", "Approximate", "Subsampling", "Random")
.shape.simplification.4        <- c(1, 0, 2, 3)
names(.shape.simplification.4) <- names(.pal.simplification.4)

.pal.simplification.5          <- c(.color.simplification, .color.approx.simplification, .color.subsampling, .color.random.grouping, .color.metacell)
names(.pal.simplification.5)   <- c("Exact", "Approximate", "Subsampling", "Random", "Metacell")
.shape.simplification.5        <- c(1, 0, 2, 3, 4)
names(.shape.simplification.5) <- names(.pal.simplification.5)

.pal.alt.clusterings           <- c(inlmisc::GetColors(9,  scheme = "muted")[c(6, 9, 4, 8)], "#8CC6EC")
names(.pal.alt.clusterings)    <- c("hcl", "hcl_weighted", "kmeans", "PAM_weighted", "Seurat")

##### tSNE colors #####
#install.packages("inlmisc")
.color.tsne                    <- c(.cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                                             "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) 
scales::show_col(.color.tsne)
#inlmisc::GetColors(9, scheme = "muted")

.color.tsne.Tian               <- .color.tsne[c(2, 3, 4, 5, 8)] #[c(1:5)]
names(.color.tsne.Tian)        <- c("A549", "H838", "H1975", "H2228", "HCC827")

.color.tsne.Carmona            <- .color.tsne[c(2:4)] #[c(4,2,6)]
names(.color.tsne.Carmona)     <- c("Exhausted/ML", "EM-like", "Naive")

.color.tsne.Zheng              <- .color.tsne[c(3, 2)]
names(.color.tsne.Zheng)       <- c("CD4", "CD8")


.color.tsne.Zheng.PBMCs       <- .color.tsne[4:8]
names(.color.tsne.Zheng.PBMCs)<- c("T",  "CD14",  "B",  "NK",  "CD34")

.color.tsne.Paul               <- c(inlmisc::GetColors(9, scheme = "muted"), "black")
names(.color.tsne.Paul)             <- c("Erythrocyte", "MP/EP", "Megakaryocyte", "GMP", "DC", "Basophil", "Monocyte", "Neutrophil", "Eosinophil","NK")

.color.tsne.Zilionis           <- c("#4867B1", "#4A2C4B", "#F8991C", "#33C3EC", "#1E6936", "#9F93FD", "#6ABD46")
names(.color.tsne.Zilionis)    <- c("B cells", "Basophils", "MoMacDC", "Neutrophils",
                                   "NK cells", "pDC", "T cells")

.gg.h <- 3
.gg.w  <- 3

.gg.h.genegene <- 3.5
.gg.w.genegene <- 3.5

.Y_lims <- c(-0.1, 1.01)

.my_ggsave <- function(p, folder, name, w = .gg.w, h = .gg.h, ext = "pdf", save.raw.ggplot = TRUE){
  ggsave(p, filename = paste0(folder, name, "_legend.", ext), width = 2*w, height = 2*h, useDingbats=FALSE)
  ggsave(p + theme(legend.position = "none"), filename = paste0(folder, name, ".", ext), width = w, height = h, useDingbats=FALSE)
  if(save.raw.ggplot)
    saveRDS(p, file = paste0(folder, name, ".Rds"))
}

.my_saveplot <- function(p, folder, name, w = .gg.w, h = .gg.h, ext = "pdf", save.raw.ggplot = TRUE, do.log = TRUE, asp = 1){
  fname         <- paste0(folder, name, ".", ext)
  fname.rds     <- paste0(folder, name, ".Rds")
  
  legend               <- cowplot::get_legend(p)
  title                <- cowplot::get_title(p)
  p.wo.legend          <- p + theme(legend.position = "none", asp = asp)
  p.wo.legend.void     <- p + theme_void() + theme(legend.position = "none", asp = asp)
  p.wo.legend.classic  <- p + theme_classic() + theme(legend.position = "none", asp = asp)
  p.wo.legend.bw       <- p + theme_bw() + theme(legend.position = "none", asp = asp, 
                                                 panel.grid.minor  = element_blank(),
                                                 panel.grid.major  = element_line(size = 0.1))
  if(do.log){
  p.list        <- list(p, 
                        ggpubr::as_ggplot(legend),
                        ggpubr::as_ggplot(title),
                        p.wo.legend,
                        p.wo.legend + scale_x_log10(),
                        p.wo.legend + scale_y_log10(),
                        p.wo.legend + scale_x_log10() + scale_y_log10(),
                        p.wo.legend.void,
                        p.wo.legend.void + scale_x_log10(),
                        p.wo.legend.void + scale_y_log10(),
                        p.wo.legend.void + scale_x_log10() + scale_y_log10(),
                        p.wo.legend.classic,
                        p.wo.legend.classic + scale_x_log10(),
                        p.wo.legend.classic + scale_y_log10(),
                        p.wo.legend.classic + scale_x_log10() + scale_y_log10(),
                        p.wo.legend.bw,
                        p.wo.legend.bw + scale_x_log10(),
                        p.wo.legend.bw + scale_y_log10(),
                        p.wo.legend.bw + scale_x_log10() + scale_y_log10())
  } else {
    p.list        <- list(p, 
                          ggpubr::as_ggplot(legend),
                          ggpubr::as_ggplot(title),
                          p.wo.legend,
                          p.wo.legend.void,
                          p.wo.legend.classic,
                          p.wo.legend.bw)
  }
    
  

  pdf(fname, width = w, height = h)
  invisible(lapply(p.list, print))
  dev.off()
  
  if(save.raw.ggplot) saveRDS(p, file = fname.rds)
}

.my_savepdf <- function(p, folder, name, h = .h.gg, w = .w.gg){
  pdf(paste0(folder, name, ".pdf"))
  for(p.i in p){
    print(p.i)
  }
  dev.off()
}

.ltype.gamma <- c("solid", "F1", "longdash", "dotted")
names(.ltype.gamma) <- as.character(c(1, 10, 100, 1000))
