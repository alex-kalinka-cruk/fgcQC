#' mock_missing_FGC_data
#'
#' Function to write out mock data that is missing from the inputs to `QC_fgc_crispr_data`. The mock data is written to temporary files that are deleted when the `R` session ends.
#'
#' @param analysis_config A path to a valid analysis config JSON file to be used as input for the AZ-CRUK CRISPR analysis pipeline.
#' @param combined_counts A path to a valid combined counts csv file (produced by the AZ-CRUK CRISPR analysis pipeline).
#' @param bcl2fastq A character string giving one or more paths to valid `bcl2fastq2` summary output JSON files (paths separated by commas). `NULL` if missing.
#' @param library A valid path to a library tsv file in which the first column gives the sgRNA sequence and the second column gives the sgRNA ID (produced by the AZ-CRUK CRISPR reference data generation pipeline). `NULL` if missing.
#'
#' @return A named list of temporary file names pointing to mocked data.
#' @importFrom jsonlite fromJSON
#' @export
mock_missing_FGC_data <- function(analysis_config, combined_counts, bcl2fastq, library){
  ret <- list()
  if(is.null(analysis_config))
    stop("mock_missing_FGC_data: 'analysis_config' must not be NULL")
  if(is.null(combined_counts))
    stop("mock_missing_FGC_data: 'combined_counts' must not be NULL")

  ## Samples from analysis config.
  tryCatch(config <- jsonlite::fromJSON(analysis_config),
                     error = function(e) stop("mock_missing_FGC_data: unable to read analysis config JSON",analysis_config,":",e))
  samples <- config$samples$name
  indexes <- config$samples$indexes
  if(all(indexes=="") || all(is.na(indexes)) || all(is.null(indexes)) || all(is.nan(indexes))){
    # Require dummy sample indexes to be populated in bcl2fastq file.
    indexes <- paste("i",1:length(indexes),sep="")
  }

  # Do we need to mock a 'qc' section in the config?
  if(!"qc" %in% names(config)){
    analysis_config <- fgcQC::mock_config_qc_section(config, indexes, samples)
  }

  if(is.null(bcl2fastq)){
    bcl2fastq <- fgcQC::mock_bcl2fastq_json(indexes, samples)
  }
  if(is.null(library)){
    library <- fgcQC::mock_library_cleanr_tsv(combined_counts)
  }
  ret <- list(analysis_config = analysis_config, bcl2fastq = bcl2fastq, library = library)
  return(ret)
}
