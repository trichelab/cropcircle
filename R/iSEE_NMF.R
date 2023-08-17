#' iSEE_NMF
#'
#' Simple wrapper function to launch an iSEE web application for an NMF'ed SCE.
#' iSEE is a Bioconductor package that you must install to run this function.
#'
#' @param     sce     A SingleCellExperiment with NMF and UMAP reducedDims.
#'
#' @return            nothing; an event loop will be initiated and iSEE() run
#'
#' @seealso   SingleCellExperiment::SingleCellExperiment
#' @seealso   RcppML::nmf
#' @seealso   uwot::umap
#' @seealso   iSEE::iSEE
#'
#' @examples 
#' if (require("iSEE")) {
#'   data(C18)
#'   iSEE_NMF(C18)
#' } 
#'
#' @export
iSEE_NMF <- function(sce) { 

  if (!require("iSEE")) {
    message("Please install necessary prerequisites before running iSEEapp:")
    message("install.packages('BiocManager')")
    message("BiocManager::install('iSEE')")
    stop("Cannot launch without `iSEE`.")
  } else { 
    message("Need to settle on a canonical 'first grouping' and autoNMF")
    browser()
    iSEE(sce,
         initial=list(UMAP=new("ReducedDimensionPlot",
                               ColorByColumnData = "ExperimentalGroup",
                               ColorBy = "Column data", 
                               Type = "UMAP"),
                      Feature=new("RowDataTable"),
                      Abundance=new("FeatureAssayPlot", 
                                    Assay = "MzR", 
                                    XAxis = "Column data", 
                                    XAxisColumnData = "ExperimentalGroup", 
                                    YAxisFeatureSource = "Feature", 
                                    YAxisFeatureDynamicSource = TRUE, 
                                    ColorByColumnData = "ExperimentalGroup", 
                                    ColorBy = "Column data") 
                      )
          )
  }

}
