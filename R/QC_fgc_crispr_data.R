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


.read_bcl2fastq <- function(b2f){
  b2f <- strsplit(bcl2fastq,",")
  ret <- NULL
  for(i in 1:length(b2f)){
    tryCatch(
      ret <- rbind(ret, fgcQC::extract_b2f_json(b2f[i])),
      error = function(e) stop(paste("unable to build bcl2fastq2 data frame:",e))
    )
  }
  return(ret)
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
#' @param comparison_name A character string naming a single comparison to extract QC data for (should correspond to the comparison name used in the analysis config JSON file).
#' @param output A character string giving an output file name for the csv results. If `NULL`, do not write out any results.
#'
#' @return A data frame containing QC metrics as columns and samples as rows; this data will also be written to the `output` file, if not `NULL`.
#' @importFrom dplyr mutate select filter
#' @export
QC_fgc_crispr_data <- function(analysis_config, combined_counts, bagel_ctrl_plasmid, bcl2fastq, library, comparison_name, output){
  # File path checks.
  .check_qc_inputs(analysis_config, combined_counts, bagel_ctrl_plasmid, bcl2fastq, library, comparison_name, output)

  # Analysis config.
  json_list <- fgcQC::read_analysis_config_json(analysis_config)
  comparisons <- fgcQC::extract_analysis_comparisons(json_list)
  if(!comparison_name %in% comparisons$name)
    stop(paste("unable to find comparison",comparison_name,"in",analysis_config))
  comparisons %<>%
    dplyr::filter(comparison == comparison_name)

  # CI sequencing bcl2fastq2 output.
  b2f <- .read_bcl2fastq(bcl2fastq) %>%
    dplyr::mutate(screen_type = comparisons$type[1],
                  screen_goal = comparisons$goal[1],
                  SampleName = json_list$samples$name[match(SampleId, json_list$samples$indexes)],
                  SampleLabel = json_list$samples$label[match(SampleId, json_list$samples$indexes)],
                  SampleClass = comparisons$class[match(SampleName, comparisons$sample)]) %>%
    # Keep only samples in the given comparison.
    dplyr::filter(SampleName %in% comparisons$sample) %>%
    dplyr::select(Flowcell,RunId,RunDate,read_lengths,LaneNumber,TotalClustersRaw,TotalClustersPF,ReadsPF_percent,
                  Non_Demultiplexed_Reads_percent,Q30_bases_samples_percent,
                  SampleId,IndexSequence,SampleName,SampleLabel,SampleClass,NumberReads,Q30_bases_percent,
                  Average_base_quality,Index_OneBaseMismatch_percent,Trimmed_bases_percent,Sample_Representation)


}
