#' read untargeted metabolomics data into R and churn out an SCE
#' 
#' Not to be confused with readMzR
#' 
#' To a first approximation, the function just reads in and parses colnames
#' 
#' @param x       the CSV file of the output 
#' @param coldat  extract column data from the column names? (TRUE)
#' @param out     format of the output object (default is "SCE"; alt is "raw")
#' @param verbose be verbose? (FALSE)
#'
#' @details 
#'   TC (technical controls) have various subtypes:
#'     TC_Blank  (solvent blank)
#'     TC_PB     (process blank)
#'     TC_QC     (QC within batch -- pooled sample -- linearity/shift)
#'   
#'   TC_*_conditioning can be dropped altogether
#'   TC_*_restart can be dropped altogether 
#' 
#'   The remaining TCs are useful for normalization (Blank) and linearization 
#'   (QC) -- linearization via GAMM is sensible -- ZachM problem -- mgcv 
#' 
#'   It is important to block on Batch (B[1234]) and Experiment (MC000XYZ) 
#'   since the QC controls are within batch linearizers for run order, while
#'   TC_PB is a true blank (noise), and TC_Blank is a solvent blank (carrier)
#'   In terms of TC_QC, it is critical to model the dropoff per feature well.
#' 
#' @examples 
#' 
#' if (FALSE) { 
#'  
#'   library(Matrix)
#'   library(RcppML)
#'   library(singlet)
#'
#'   fl <- "MC00295_HILIC_A_All_PoolSize.csv"
#'   HILAC <- readUntargeted(fl)
#'   HILAC <- HILAC[-which(rownames(HILAC) == "group"), ] 
#'   HILAC <- Matrix(log1p(HILAC)) # unstandardized; fixme
#'   nmf_fit <- RcppML::nmf(HILAC, k=5) # need ard_nmf
#'
#' }
#'
#' @import Matrix
#' @import SingleCellExperiment
#' 
#' @export
readUntargeted <- function(x, coldat=TRUE, out=c("SCE","raw"), verbose=FALSE) {
  
  out <- match.arg(out)
  res <- read.csv(x, head=FALSE)
  colnames(res) <- res[1, ]
  res <- res[-1, ] 
  rownames(res) <- make.unique(res[, 1])
  res <- res[, -1]
  firstHit <- head(grep("^(Exp|TC|CAL)_", ignore=TRUE, colnames(res)), 1)
  lastHit <- tail(grep("^(Exp|TC|CAL)_", ignore=TRUE, colnames(res)), 1)
  rowDatCols <- colnames(res)[seq(1, (firstHit - 1))]
  if (ncol(res) > lastHit) {
    addCols <- colnames(res)[seq((lastHit + 1), ncol(res))]
    rowDatCols <- c(rowDatCols, addCols)
  }
  if (length(rowDatCols) > 0) {
    if (verbose) message("Separating rowData from assay data...", append=FALSE)
    rowdat <- res[, rowDatCols]
    assayDatCols <- setdiff(colnames(res), rowDatCols)
    res <- data.matrix(res[, assayDatCols])
    if (verbose) message("...done.")
  } 
  if (coldat) {
    if (verbose) message("Extracting colData from column names...",append=FALSE)
    cdat <- parseColnames(colnames(res)) 
    if (verbose) message("...done.")
  }

  if (out == "SCE") {
    if (verbose) message("Constructing SingleCellExperiment...", append=FALSE)
    sce <- SingleCellExperiment(assays=list(MzR=as(res, "dgCMatrix")))
    if (coldat) colData(sce) <- DataFrame(cdat)
    if (exists("rowdat")) rowData(sce) <- rowdat
    mainExpName(sce) <- "MzR"
    return(sce)
    if (verbose) message("...done.")
  } else { 
    if (coldat) attr(res, "colData") <- cdat
    return(res)
  }

}
