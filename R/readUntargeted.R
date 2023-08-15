#' read untargeted metabolomics data into R and churn out an SCE
#' 
#' Not to be confused with readMzR
#' 
#' To a first approximation, the function just reads in and parses colnames
#' 
#' @param x       the CSV file of the output 
#' @param coldat  extract column data from the column names? (FALSE)
#' @param out     format of the output object (default is "raw"; alt "SCE")
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
#'   HILAC <- Matrix(log1p(HILAC)) # unstandardized 
#'   nmf_fit <- singlet::ard_nmf(HILAC)
#'
#' }
#'
#' @import Matrix
#' @import SingleCellExperiment
#' 
#' @export
readUntargeted <- function(x, coldat=TRUE, out=c("SCE","raw")) {
  
  res <- read.csv(x, head=FALSE)
  colnames(res) <- res[1, ]
  res <- res[-1, ] 
  rownames(res) <- make.unique(res[, 1])
  firsthit <- grep("Exp", colnames(res))[1]
  rdn <- colnames(res)[seq(1, (firsthit - 1))]
  if (length(rdn) > 0) {
    rowdat <- res[, rdn]
    adn <- setdiff(colnames(res), rdn)
    res <- as.matrix(res[, adn])
  } 
  out <- match.arg(out)
  
  if (coldat) {
    cdat <- parseColnames(colnames(res)) 
    attr(res, "colData") <- cdat
  }

  if (out == "SCE") {
    sce <- SingleCellExperiment(assays=list(MzR=as.matrix(res)))
    if (coldat) colData(sce) <- DataFrame(cdat)
    if (exists("rowdat")) rowData(sce) <- rowdat
    mainExpName(sce) <- "MzR"
    return(sce)
  } else { 
    return(res)
  }

}
