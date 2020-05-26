# Helper functions.
.file_exists <- function(file){
  if(!file.exists(file))
    stop(paste("unable to find",file))
}


.check_qc_inputs <- function(analysis_config, combined_counts, bagel_ctrl_plasmid, bcl2fastq, library, comparison, output, norm_method){
  .file_exists(analysis_config)
  .file_exists(combined_counts)
  .file_exists(bagel_ctrl_plasmid)
  .file_exists(library)
  if(!is.character(bcl2fastq))
    stop("'bcl2fastq' must be a character string pointing to one or more bcl2fastq2 output JSON files separated by commas")
  b2f <- strsplit(bcl2fastq,",")
  lapply(b2f, .file_exists)
  if(!is.character(comparison))
    stop("'comparison_name' must be a character string naming a comparison in the analysis config JSON file")
  if(!is.character(output) & !is.null(output))
    stop("'output' must be a character string naming an output file for QC results or NULL")
  if(!norm_method %in% c("median_ratio","relative"))
    stop("'norm_method' should be one of 'median_ratio' or 'relative'")
}


.read_bcl2fastq <- function(bcl2fastq){
  b2f <- unlist(strsplit(bcl2fastq,","))
  ret <- NULL
  for(i in 1:length(b2f)){
    tryCatch(
      ret <- rbind(ret, fgcQC::extract_b2f_json(b2f[i])),
      error = function(e) stop(paste("unable to read bcl2fastq2 data:",e))
    )
  }
  return(ret)
}


.check_comparison <- function(comparisons, comparison_name, analysis_config){
  if(!comparison_name %in% comparisons$comparison)
    stop(paste("unable to find comparison",comparison_name,"in",analysis_config))
  if(!"plasmid" %in% comparisons$class)
    stop(paste("unable to find any 'plasmid' samples in",comparison_name))
  if(!"control" %in% comparisons$class)
    stop(paste("unable to find any 'control' samples in",comparison_name))
  if(comparisons$goal[1] != "lethality"){
    if(!"treatment" %in% comparisons$class)
      stop(paste("unable to find any 'treatment' samples in",comparison_name))
  }
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
#' @param norm_method A character string naming a normalization method for the count data. Can be `median_ratio` or `relative`. Defaults to `median_ratio`.
#'
#' @return A data frame containing QC metrics as columns and samples as rows; this data will also be written to the `output` file, if not `NULL`.
#' @importFrom dplyr mutate select filter right_join
#' @export
QC_fgc_crispr_data <- function(analysis_config, combined_counts, bagel_ctrl_plasmid, bcl2fastq, library,
                               comparison_name, output, norm_method = "median_ratio"){
  ## File path checks.
  .check_qc_inputs(analysis_config, combined_counts, bagel_ctrl_plasmid, bcl2fastq, library, comparison_name, output, norm_method)

  ## Analysis config.
  json_list <- fgcQC::read_analysis_config_json(analysis_config)
  comparisons <- fgcQC::extract_analysis_comparisons(json_list)

  comparisons %<>%
    dplyr::filter(comparison == comparison_name)
  .check_comparison(comparisons, comparison_name, analysis_config)

  # Prep qc data from analysis config.
  qc_config <- json_list$qc %>%
    dplyr::select(-indexes, -label) %>%
    dplyr::filter(name %in% comparisons$sample)

  ## CI sequencing bcl2fastq2 output.
  tryCatch({
    b2f <- .read_bcl2fastq(bcl2fastq) %>%
      dplyr::mutate(screen_type = comparisons$type[1],
                    screen_goal = comparisons$goal[1],
                    SampleName = json_list$samples$name[match(SampleId, json_list$samples$indexes)],
                    SampleLabel = json_list$samples$label[match(SampleId, json_list$samples$indexes)]) %>%
      # Keep only samples in the given comparison.
      dplyr::filter(SampleName %in% comparisons$sample) %>%
      dplyr::right_join(qc_config, by = c("SampleName" = "name")) %>%
      dplyr::mutate(SampleClass = comparisons$class[match(SampleName, comparisons$sample)],
                    SampleLabel = ifelse(SampleClass == "plasmid",SampleName,SampleLabel)) %>%
      dplyr::select(slx_id,Flowcell,RunId,RunDate,read_lengths,LaneNumber,TotalClustersRaw,TotalClustersPF,ReadsPF_percent,
                    Non_Demultiplexed_Reads_percent,Q30_bases_samples_percent,
                    virus_batch,plasmid_batch,cas_activity,minimum_split_cell_number,
                    SampleId,IndexSequence,SampleName,SampleLabel,SampleClass,NumberReads,Q30_bases_percent,
                    Average_base_quality,Index_OneBaseMismatch_percent,Trimmed_bases_percent,Sample_Representation,
                    cell_population_doublings,index_plate,index_plate_well_id,library_dna_yield_ng.ul)
  },
  error = function(e) stop(paste("unable to build bcl2fastq data frame:",e))
  )

  return(b2f)

  # Besides the plasmid (usually won't be sequenced), are any samples missing?
  miss_samps <- setdiff(comparisons$sample[comparisons$class!="plasmid"], b2f$SampleName)
  if(length(miss_samps) > 0){
    stop(paste("could not find the following comparison samples in the sequencing output (bcl2fastq):",miss_samps))
  }

  # Output determined by screen type.
  if(comparisons$type[1] == "n"){
    ## Read and normalize counts.
    tryCatch(counts <- read.delim(combined_counts, sep="\t", header=T, stringsAsFactors = F),
             error = function(e) stop(paste("unable to read combined counts file",combined_counts,":",e)))
    if(norm_method == "median_ratio"){
      counts_norm <- fgcQC::normalize_library_depth_median_ratio(counts)
    }else{
      counts_norm <- fgcQC::normalize_library_depth_relative(counts, 2e7)
    }

    ## Log fold change at gRNA and gene levels.
    control_samp <- comparisons$sample[comparisons$class == "control"]
    plasmid_samp <- comparisons$sample[comparisons$class == "plasmid"]
    # 1. Control vs Plasmid.
    lfc.ctrl_pl <- fgcQC::calc_log2_fold_change_gRNAs(counts, ref = plasmid_samp, comp = control_samp)
    lfc.ctrl_pl.genes <- fgcQC::calc_log2_fold_change_genes(lfc.ctrl_pl)
  }else{

  }

}
