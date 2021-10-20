#' Retreives and saves GO data
#' 
#' @param genes vector of genes for which get and save GO ids
#' @param martDS mart data set 
#' @param data.folder folder to save GO dataset 
#' @param DO.save whether to save GO dataset
#' 
#' @return \link[biomaRt]{getBM} result for genes
#' 
#' @export

get_GO <- function(
  genes, 
  martDS = "hsapiens_gene_ensembl",
  data.folder = "/data",
  go_domain = "all",
  DO.save = TRUE
  
){
  
  if(go_domain == "all" | is.null(go_domain)){
    go_domain <- c('biological_process','molecular_function','cellular_component')
  }
  
  ensembl <- biomaRt::useMart("ensembl", dataset = martDS)
  bm      <- biomaRt::getBM(attributes=c('hgnc_symbol', 'go_id', 'name_1006', 'namespace_1003'), filters = 'hgnc_symbol', values= genes, mart = ensembl)
  bm      <- bm[bm$namespace_1003 %in% go_domain,]
  
  if(DO.save){
    if(!dir.exists(data.folder)) dir.create(data.folder, recursive = TRUE)
    saveRDS(bm, file = file.path(data.folder, "BM.Rds"))
  }
  return(bm)
}
