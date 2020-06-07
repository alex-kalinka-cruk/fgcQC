#' mock_config_qc_section
#'
#' Mock an analysis config write to a temporary JSON file, and return the temporary file path.
#'
#' @param config A list of config elements retrieved by reading in the analysis config JSON.
#' @param sample_ids A character vector naming the sample IDs (indexes) to be mocked.
#' @param sample_names A character vector of sample names extracted from the analysis config comparisons.
#'
#' @return A file path to a temporary mock analysis config JSON file.
#' @importFrom dplyr slice n mutate
#' @importFrom magrittr %<>%
#' @importFrom jsonlite write_json
#' @export
mock_config_qc_section <- function(config, sample_ids, sample_names){
  tryCatch({
    num_samps <- length(sample_ids)
    if(!"qc" %in% names(config)){
      qc <- fgcQC::config_qc_metrics %>%
        dplyr::slice(rep(1:dplyr::n(), each = num_samps)) %>%
        dplyr::mutate(name = sample_names, indexes = sample_ids)
      config$qc <- qc
    }else if(!"date_transduced" %in% colnames(config$qc)){
      config$qc %<>%
        dplyr::mutate(date_transduced = "NA")
    }
    # Are we missing a 'label' column in the config 'samples'?
    if(!"label" %in% colnames(config$samples)){
      config$samples %<>%
        dplyr::mutate(label = paste("S",1:dplyr::n(),sep=""))
    }
    config_file <- tempfile("mock_analysis_config.")
    jsonlite::write_json(config, config_file)
  },
  error = function(e) stop(paste("mock_config_qc_section: unable to mock qc section in analysis config:",e))
  )
  return(config_file)
}
