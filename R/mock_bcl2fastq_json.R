#' mock_bcl2fastq_json
#'
#' Mock bcl2fastq data for a specified number of samples with given sample IDs (to match an analysis config), write to a temporary file, and return the temp file path.
#'
#' @param sample_ids A character vector naming the sample IDs to be mocked.
#'
#' @return A file path to a temporary mock 'bcl2fastq' JSON file.
#' @importFrom dplyr slice n mutate
#' @importFrom magrittr %<>%
#' @importFrom jsonlite write_json
#' @export
mock_bcl2fastq_json <- function(sample_ids){
  tryCatch({
    num_samps <- length(sample_ids)
    mock_bcl <- fgcQC::mock_bcl2fastq_template
    mock_bcl$ConversionResults$DemuxResults[[1]] %<>%
      # Expand number of samples of template bcl2fastq data.
      dplyr::slice(rep(1:dplyr::n(), each = num_samps)) %>%
      # Replace sample IDs.
      dplyr::mutate(SampleId = sample_ids)
    mock_bcl_file <- tempfile("mock_bcl2fastq.")
    jsonlite::write_json(mock_bcl, mock_bcl_file)
  },
  error = function(e) stop(paste("mock_bcl2fastq_json: unable to mock bcl2fastq data:",e))
  )
  return(mock_bcl_file)
}
