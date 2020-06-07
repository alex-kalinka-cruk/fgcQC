#' read_analysis_config_json
#'
#' Extracts sample and comparison data from an analysis config JSON file to be used as input for the AZ-CRUK CRISPR analysis pipeline.
#'
#' @param file A path to a valid analysis config JSON file.
#'
#' @return A list containing all the sections in the JSON file.
#' @author Alex T. Kalinka, \email{alex.kalinka@@cancer.org.uk}
#' @importFrom jsonlite fromJSON
#' @export
read_analysis_config_json <- function(file){
  expected_sections <- c("samples","comparisons","meta","general","qc")

  tryCatch(jlist <- jsonlite::fromJSON(file),
           error = function(e) stop("read_analysis_config_json: unable to process file",file,":",e))

  miss <- setdiff(expected_sections, names(jlist))
  if(length(miss) > 0)
    stop(paste("read_analysis_config_json: analysis config is missing expected sections:",paste(miss,collapse=",")))

  # Do we have a 'date_transduced' column in the qc section?
  if(!"date_transduced" %in% colnames(jlist$qc))
    stop(paste("read_analysis_config_json: unable to find a 'date_transduced' column in the qc section"))

  return(jlist)
}
