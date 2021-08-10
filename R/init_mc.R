#' Initiate MC folders
#'
#' This function initiate MetaCell folder by running \code{scdb_init()}, \code{scfigs_init()}, \code{set_param()}
#'
#' @param mc.dir folder for MC data and outputs
#'
#' @return void
#'
#' @examples
#' \dontrun{
#'
#' }
#' @export
#'
init_mc <- function(mc.dir = NULL){
  if(is.null(mc.dir)){
    stop("specify path to metacell output folder")
  }

  if(!dir.exists(mc.dir)) dir.create(mc.dir, recursive = TRUE)
  metacell::scdb_init(mc.dir, force_reinit=T)

  if(!dir.exists(paste0(mc.dir, "/figs"))) dir.create(paste0(mc.dir, "/figs"))
  metacell::scfigs_init(paste0(mc.dir, "/figs"))
  tgconfig::set_param("mcell_mc2d_height", 1000, "metacell")
  tgconfig::set_param("mcell_mc2d_width", 1000, "metacell")

  return()
}
