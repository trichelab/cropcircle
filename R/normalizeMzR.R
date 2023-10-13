#' normalize metabolomics data with 'Exp', 'TC', and perhaps 'CAL' columns
#' 
#' assay(sce, "normMzR") <- log1p((Exp|CAL) / TC) for all features in MzR
#'
#' @param sce     a SingleCellExperiment with !is.null(assay(sce, "MzR"))
#' @param norm    name of the column indicating normalization controls ("norm")
#' @param batch   name of the column indicating batch (for matching; "batch")
#' @param how     how to normalize (default: "logFC", see Details)
#' 
#' @return        a SingleCellExperiment normalized per batch against controls
#' 
#' @details 
#' logFC normalization is log((1 + (X/sigma_norm))/(1 + (mu_norm/sigma_norm))).
#' sigma_norm is padded by 1 to avoid division by zero (it shrinks, a bit). 
#' Be careful not to use loess controls as if they were blanks!
#'
#' @import Matrix
#' @import SingleCellExperiment
#' 
#' @export
normalizeMzR <- function(sce, norm="norm", batch="batch", how=c("logFC")) {

  how <- match.arg(how) 
  if (is.null(assay(sce, "MzR"))) stop("MzR assay not found. Cannot normalize.")
  batches <- unique(colData(sce)[[batch]])
  batches <- sort(batches[!is.na(batches)])
  names(batches) <- batches

  batchnormed <- lapply(batches, function(b) {
    batchsce <- sce[, which(sce$batch == b)]
    mu <- rowMeans(assay(batchsce[,which(colData(batchsce)[[norm]])],"MzR"))
    sigma <- rowSds(assay(batchsce[,which(colData(batchsce)[[norm]])],"MzR"))+1
    muZ <- mu/sigma
    xZ <- assay(batchsce, "MzR")/sigma
    assay(batchsce, "normMzR") <- log( (1+xZ)/(1+muZ))
    return(batchsce)
  })
  do.call(cbind, batchnormed)

}
