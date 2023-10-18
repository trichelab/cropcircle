#' parse MNP metabio standardized column names
#' 
#' e.g. 'Exp_D9-8W-LN_2_C18neg-MC00138-S028'
#'
#' @param cn            column names
#' @param sep1          separator 1 ("_")
#' @param sep2          separator 2 ("-")
#' @param full          attempt to extract full sample names? (TRUE)
#' @param expCon        value of $condition for non-control specimens ("Exp")
#' 
#' @return              a data.frame (coercible to DataFrame for SCE objects)
#' 
#' @details 
#' In the provided 'Exp_D9-8W-LN_2_C18neg-MC00138-S028', splitting on '_':
#'   token 1 is condition (Exp, TC, or CAL)
#'   token 2 is experimental factors, further split by '-':
#'     factor1:'D9'
#'     factor2:'8W'
#'     factor3:'LN'
#'   token 3 is replicate number (2)
#'   token 4 is identifiers, further split by '-':
#'     method:'C18neg'
#'     experiment:'MC00138'
#'     sample:'S028'
#' 
#' names(colData) for a default SCE (with fullSample=FALSE) would then be 
#' 
#' c("condition", "method", "experiment", "sample",
#'   "factor1", "factor2", "factor3", "replicate")
#'
#' If full == TRUE, then a further check is performed to extract complete
#' sample names for each control (i.e., $condition != expCon) sample. 
#'
#' c("condition", "method", "experiment", "sample", 
#'   "factor1", "factor2", "factor3", "replicate", 
#'   "originalName", "fullSample")
#'
#' This is the default behavior, but it may not always be desirable.
#' 
#' @examples  
#' cn <- c("Exp_D9-8W-LN_1_C18neg-MC00138-S027", 
#'         "Exp_D9-8W-LN_3_C18neg-MC00298-S029",
#'         "TC_Blank-conditioning_1_C30pos-MC00138-B3",
#'         "TC_QC_1_PHneg-MC00138-B2",
#'         "Exp_GFP-16w-NC_4_T3pos-MC00138-S053",
#'         "TC_QC_3_T3pos-MC00298-B1")
#'         
#' parseColnames(cn, fullSample=FALSE)
#' parseColnames(cn, fullSample=TRUE) # default
#'
#' @export
#'
parseColnames <- function(cn, sep1="_", sep2="-", full=TRUE, expCon="Exp"){

  if (is(cn, "data.frame")) cn <- names(cn)
  if (is(cn, "matrix")) cn <- colnames(cn)
  tokens <- do.call(rbind, strsplit(cn, sep1, fixed=TRUE))
  colnames(tokens) <- c("condition", "factors", "replicate", "ids")
  factors <- do.call(rbind, strsplit(tokens[, "factors"], sep2, fixed=TRUE))
  colnames(factors) <- paste0("factor", seq_len(ncol(factors)))
  ids <- do.call(rbind, strsplit(tokens[, "ids"], sep2, fixed=TRUE))
  colnames(ids) <- c("method", "experiment", "sample")
  cdat <- cbind(condition=tokens[, "condition"],
                as.data.frame(ids), 
                as.data.frame(factors), 
                replicate=tokens[,"replicate"])
  rownames(cdat) <- cn

  # rename better
  if (fullSample) { 
    
    cdat$originalName <- rownames(cdat)
    stopifnot("experiment" %in% colnames(cdat))
    cdat$experiment <- as.character(cdat$experiment)
    stopifnot("condition" %in% colnames(cdat))
    cdat$fullSample <- 
      with(cdat, 
           ifelse(condition == expCon, 
                  paste0(experiment,":",sample),
                  paste0(experiment,":",sample,"_",factor1,"_",replicate)))
 
  }

  # either way
  return(cdat)

}
