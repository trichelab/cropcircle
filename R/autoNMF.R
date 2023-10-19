#' take an SCE of metabolomics data with 'MzR' in the assays, norm, and NMF it
#' 
#' This function assumes that TC and CAL samples have been used to normalize, 
#' but if they have not (if there is no normMzR in assayNames(sce)), it will
#' run normalizeMzR() and the latter will complain if no TCs/CALs are present 
#' 
#' FIXME: drop the dependency on singlet!
#'
#' @param x       a SingleCellExperiment (or similar) of mass spec data
#' @param k       target rank for NMF decomposition (default: use ard)
#' @param asy     name of the assay to use (default: normMzR) 
#' @param ...     additional arguments passed to RcppML::nmf or singlet::ard_nmf
#' 
#' @return        an RcppML::nmf object
#'
#' @import singlet
#' @import RcppML
#' @import Matrix
#' @import SingleCellExperiment
#' 
#' @export 
#'
autoNMF <- function(x, k=NULL, asy="normMzR", ...) { 

  if (is.null(k)) {
    as(ard_nmf(assay(x, asy), ...), "nmf")
  } else { 
    nmf(assay(x, asy), k=k, ...)
  }

}
