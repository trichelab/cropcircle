#' accessor for assays(x, 'MzR')
#'
#' @param   x     the object
#' 
#' @return        assay(sce, "MzR") 
#'
#' @export
#'
MzR <- function(x) {
  assay(x, "MzR")
}
