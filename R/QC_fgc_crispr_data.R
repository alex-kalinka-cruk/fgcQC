# Helper functions.
.file_exists <- function(file){
  if(!file.exists(file))
    stop(paste("unable to find",file))
}

.check_qc_inputs <- function(analysis_config, combined_counts, bagel_ctrl_plasmid, bcl2fastq, library, comparison, output){
  .file_exists(analysis_config)
  .file_exists(combined_counts)
  .file_exists(bagel_ctrl_plasmid)
  .file_exists(library)
  if(!is.character(bcl2fastq))
    stop("'bcl2fastq' must be a character string pointing to one or more bcl2fastq2 output JSON files separated by commas")
  b2f <- strsplit(bcl2fastq,",")
  lapply(b2f, .file_exists)
  if(!is.character(comparison))
    stop("'comparison' must be a character string naming a comparison in the analysis config JSON file")
  if(!is.character(output) & !is.null(output))
    stop("'output' must be a character string naming an output file for QC results or NULL")
}

#' QC_fgc_crispr_data
#'
#' Produce Quality Control (QC) metrics from FGC CRISPR screen data inputs: an analysis config JSON file, a combined counts csv file, a Bagel results tsv file (Control vs Plasmid), a gRNA library tsv file ('cleanr.tsv'), and one or more bcl2fastq2 output JSON files. QC output will be generated for a single named comparison.
#'
#' @param analysis_config A path to a valid analysis config JSON file to be used as input for the AZ-CRUK CRISPR analysis pipeline.
#' @param combined_counts A path to a valid combined counts csv file (produced by the AZ-CRUK CRISPR analysis pipeline).
#' @param bagel_ctrl_plasmid A path to a valid Bagel tsv results file for Control vs Plasmid.
#' @param bcl2fastq A character string giving one or more paths to valid `bcl2fastq2` summary output JSON files (paths separated by commas).
#' @param library A valid path to a library tsv file in which the first column gives the sgRNA sequence and the second column gives the sgRNA ID (produced by the AZ-CRUK CRISPR reference data generation pipeline).
#' @param comparison A character string naming a single comparison to extract QC data for (should correspond to the comparison name used in the analysis config JSON file).
#' @param output A character string giving an output file name for the csv results. If `NULL`, do not write out any results.
#'
#' @return A data frame containing QC metrics as columns and samples as rows; this data will also be written out to the `output` file, if not `NULL`.
#' @importFrom dplyr mutate
#' @export
QC_fgc_crispr_data <- function(analysis_config, combined_counts, bagel_ctrl_plasmid, bcl2fastq, library, comparison, output){
  .check_qc_inputs(analysis_config, combined_counts, bagel_ctrl_plasmid, bcl2fastq, library, comparison, output)
  # Analysis config.
  json_list <- fgcQC::read_analysis_config_json(analysis_config)
  json_comp <- fgcQC::extract_analysis_comparisons(json_list)
  if(!comparison %in% json_comp$comparisons)
    stop(paste("unable to find",comparison,"in",analysis_config))
  return(json_list)
}
