#' normalize metabolomics data with 'Exp', 'TC', and perhaps 'CAL' columns
#' 
#' assay(sce, "normMzR") <- log1p((Exp|CAL) / TC) for all features in MzR
#'
#' @param sce     a SingleCellExperiment with !is.null(assay(sce, "MzR"))
#' 
#' @import Matrix
#' @import SingleCellExperiment
#' 
#' @export
normalizeMzR <- function(sce) {


}
