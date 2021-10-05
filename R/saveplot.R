#' Saves plot for manuscript
#'
#' @param p ggplot object
#' @param folder figure folder
#' @param name figure name
#' @param w widht
#' @param h height
#' @param ext file expension
#' @param save.raw.ggplot flag whether to save raw ggplot
#' @param do.log do log normalization
#' @param asp plot aspect ratio
#'
#' @export

SCBM_saveplot <- function(p,
                          folder,
                          name,
                          w = 3,
                          h = 3,
                          ext = "pdf",
                          save.raw.ggplot = TRUE,
                          do.log = TRUE,
                          asp = 1,
                          gamma_field = 'Gamma',
                          do_filter_x_breaks = TRUE){
  fname         <- file.path(folder, paste0(name, ".", ext))
  fname.rds     <- file.path(folder, paste0(name, ".Rds"))

  legend               <- cowplot::get_legend(p)
  title                <- cowplot::get_title(p)
  p.wo.legend          <- p + ggplot2::theme(legend.position = "none", asp = asp)
  p.wo.legend.void     <- p + ggplot2::theme_void() + ggplot2::theme(legend.position = "none", asp = asp)
  p.wo.legend.classic  <- p + ggplot2::theme_classic() + ggplot2::theme(legend.position = "none", asp = asp)
  p.wo.legend.bw       <- p + ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none", asp = asp,
                   panel.grid.minor  = ggplot2::element_blank(),
                   panel.grid.major  = ggplot2::element_line(size = 0.1))
  if(do.log){
    if(!(gamma_field %in% names(p$data))){
      gamma_field <- rlang::quo_get_expr(p$mapping$x)
    }
    x_breaks <- sort(unique(p$data[[gamma_field]]))

    ## filter breaks
    if(do_filter_x_breaks){
      x_breaks <- x_breaks[!(x_breaks > 10 & x_breaks%%10 != 0)]
    }


    p.list        <- list(p,
                          ggpubr::as_ggplot(legend),
                          ggpubr::as_ggplot(title),
                          p.wo.legend,
                          p.wo.legend + ggplot2::scale_x_log10(breaks = x_breaks),
                          p.wo.legend + ggplot2::scale_y_log10(),
                          p.wo.legend + ggplot2::scale_x_log10(breaks = x_breaks) + ggplot2::scale_y_log10(),
                          p.wo.legend.void,
                          p.wo.legend.void + ggplot2::scale_x_log10(breaks = x_breaks),
                          p.wo.legend.void + ggplot2::scale_y_log10(),
                          p.wo.legend.void + ggplot2::scale_x_log10(breaks = x_breaks) + ggplot2::scale_y_log10(),
                          p.wo.legend.classic,
                          p.wo.legend.classic + ggplot2::scale_x_log10(breaks = x_breaks),
                          p.wo.legend.classic + ggplot2::scale_y_log10(),
                          p.wo.legend.classic + ggplot2::scale_x_log10(breaks = x_breaks) + ggplot2::scale_y_log10(),
                          p.wo.legend.bw,
                          p.wo.legend.bw + ggplot2::scale_x_log10(breaks = x_breaks),
                          p.wo.legend.bw + ggplot2::scale_y_log10(),
                          p.wo.legend.bw + ggplot2::scale_x_log10(breaks = x_breaks) + ggplot2::scale_y_log10())
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
