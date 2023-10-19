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
  if (verbose) message("Reading ", x, "...", appendLF=FALSE)
  res <- read.csv(x, head=FALSE)
  if (verbose) message("...done.")
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
    rowdat <- res[, rowDatCols]
    nrdc <- length(rowDatCols)
    if (verbose) message("Extracting ", nrdc, " rowData cols...",appendLF=FALSE)
    assayDatCols <- setdiff(colnames(res), rowDatCols)
    res <- data.matrix(res[, assayDatCols])
    if (verbose) message("...done.")
  } 
  if (coldat) {
    if (verbose) message("Extracting colData...", appendLF=FALSE)
    cdat <- parseColnames(colnames(res)) 
    if (verbose) message("...done.")
  }

  if (out == "SCE") {
    
    if (verbose) message("Creating SingleCellExperiment...")
    # sce <- SingleCellExperiment(assays=list(MzR=as(res, "dgCMatrix")))
    sce <- SingleCellExperiment(assays=list(MzR=res))

    # this should probably move into parseColnames
    if (coldat) {
      cdat$condition <- factor(cdat$condition) # TC or Exp
      ncontrols <- table(cdat$condition)["TC"]
      if (verbose) message("Found ", ncontrols, " controls.")
      nexperimental <- table(cdat$condition)["Exp"]
      if (verbose) message("Found ", nexperimental, " experimental samples.")
      others <- setdiff(levels(cdat$condition), c("TC", "Exp"))
      nothers <- length(others)
      if (verbose & nothers > 0) {
        message("Found ",nothers," other conditions:")
        for (i in others) message("  ", i)
      }
      cdat$experiment <- factor(cdat$experiment)
      if (verbose & nlevels(cdat$experiment) > 1) {
        message("Found ", nlevels(cdat$experiment), " experiments:")
        for (i in levels(cdat$experiment)) message("  ", i)
      }
      cdat$method <- factor(cdat$method)
      if (verbose & nlevels(cdat$method) > 1) {
        message("Found ", nlevels(cdat$method), " methods:")
        for (i in levels(cdat$method)) message("  ", i)
      }
      cdat$control <- cdat$condition == "TC"
      cdat$batch <- NA
      cdat[, "batch"] <- paste0(as.character(cdat$experiment), ":", 
                                as.character(cdat$method),
                                ifelse(grepl("^B", cdat$sample), 
                                       paste0(":", cdat$sample), ""))
      cdat$batch <- factor(cdat$batch)
      if (verbose & nlevels(cdat$batch) > 0) {
        message("Found ", nlevels(cdat$batch), " batches:")
        for (i in levels(cdat$batch)) message("  ", i)
      }
      cdat$replicate <- as.integer(cdat$replicate)
      colData(sce) <- DataFrame(cdat)
    }

    if (exists("rowdat")) {
      names(rowdat) <- sub("RT_min", "rtime", names(rowdat)) # for xcms
      # otherwise have to specify rtime="RT_min" in groupFeatures(sce, ...) 
      # e.g. groupFeatures(merged, param=SimilarRtimeParam(10), rtime="RT_min")
      rowData(sce) <- rowdat
    }
    mainExpName(sce) <- "MzR"
    if (verbose) message("...done.")
    return(sce)

  } else { 
    
    if (verbose) message("Returning the raw data matrix.")
    if (coldat) attr(res, "rowData") <- rowdat
    if (coldat) attr(res, "colData") <- cdat
    return(res)

  }

}
