#' add elution method name to rownames() to disambiguate results
#' 
#' This function will first search in metadata(sce)[[name]] for a value. 
#' If none is found, the function will search for unique(colData(sce)[[name]]).
#' If neither works, the function will return an error. 
#'
#' @param sce     a SingleCellExperiment with information about elution method.
#' @param name    name of column or item indicating elution method ("method").
#' @param sep     separator between method and peak name (":").
#' 
#' @return        a SingleCellExperiment with renamed rownames (or an error). 
#'
#' @import        SingleCellExperiment
#' 
#' @export
addMethodToPeaks <- function(sce, name="method", sep=":") {

  if (name %in% names(metadata)) {
    m <- metadata(sce)[[name]]
    if (is.null(m)) {
      stop("Error: metadata(sce)$", name, " is NULL.")
    }
  } else if (name %in% names(colData(sce))) {
    m <- unique(colData(sce)[[name]])
    if (length(m) > 1) {
      stop("Error: more than one value for sce$", name, ".")
    }
  } else {
    stop("Error: no colData or metadata item named ", name, " was found.")
  }

  if (!all(grepl(m, rownames(sce)))) {
    rownames(sce) <- paste(m, rownames(sce), sep=sep)
  }
  return(sce)

}
