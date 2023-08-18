#' parse MNP metabio standardized column names
#' 
#' e.g. 'Exp_D9-8W-LN_2_C18neg-MC00138-S028'
#' splitting on '_':
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
#' names(colData) for a default SCE would then be 
#' 
#' c("experiment","method","sample","replicate","factor1","factor2","factor3")
#' 
#' @param cn    column names
#' @param sep1  separator 1 ("_")
#' @param sep2  separator 2 ("-")
#' 
#' @return      a data.frame (coercible to DataFrame for SCE objects)
#' 
#' @examples  
#' cn <- c("Exp_D9-8W-LN_1_C18neg-MC00138-S027", 
#'         "Exp_D9-8W-LN_2_C30neg-MC00138-S028")
#'         "Exp_D9-8W-LN_3_C18pos-MC00138-S029")
#'         "Exp_D9-8W-LN_4_308pos-MC00138-S030")
#' parseColnames(cn) 
#'
#' @export
parseColnames <- function(cn, sep1="_", sep2="-") {

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
  return(cdat)

}
