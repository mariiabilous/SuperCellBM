#' Plots and saves a legend mapping cell types and colors
#'
#' @param pal.GT a palette as a named vector with names being cell types and values being corresponding colors
#' @param fig.folder folder for figures, usually has a format './plots/project_name'
#'
#' @return void, saves a legend plot to file.path(fig.folder, "Legend_cell_type.pdf")
#'
#' @export
plot_legend_celltype <- function(pal.GT, fig.folder){
  N.clusters <- length(pal.GT)
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(aes(x = 1:N.clusters, y = 1:N.clusters, color = names(pal.GT))) +
    ggplot2::scale_color_manual(values = pal.GT, name = "Cell type") +
    ggplot2::theme_bw()

  legend <- ggpubr::as_ggplot(cowplot::get_legend(p))
  if(!dir.exists(fig.folder)){
    dir.create(fig.folder)
  }
  ggplot2::ggsave(legend, filename = file.path(fig.folder, "Legend_cell_type.pdf"))
}
