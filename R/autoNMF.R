#' take an SCE of metabolomics data with 'MzR' in the assays, norm, and NMF it
#' 
#' This function assumes that TC and CAL samples have been used to normalize, 
#' but if they have not (if there is no normMzR in assayNames(sce)), it will
#' run normalizeMzR() and the latter will complain if no TCs/CALs are present 
#' 
#' @param sce     a SingleCellExperiment of mass spec data
#' @param k       target rank for NMF decomposition 
#' @param ...     additional arguments passed to RcppML::nmf
#' 
#' @import RcppML
#' @import Matrix
#' @import SingleCellExperiment
#' 
#' @export 
autoNMF <- function(sce, k, ...) { 

  

}
