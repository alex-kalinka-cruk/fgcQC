#' merge_bcl2fastq_config_qc
#'
#' Merge bcl2fastq sequencing metrics with the QC section in the analysis config.
#'
#' @param b2f A data frame of `bcl2fastq` sequencing metrics.
#' @param qc_config A data frame of QC metrics taken from the 'qc' section of the analysis config JSON.
#' @param comparisons A data frame of comparisons taken from the 'comparisons' section of the analysis config JSON.
#' @return A data frame.
#' @importFrom dplyr right_join mutate select everything
#' @export
merge_bcl2fastq_config_qc <- function(b2f, qc_config, comparisons){
  tryCatch({
    qc_metrics <- b2f %>%
      # Fold in QC data from analysis config.
      dplyr::right_join(qc_config, by = c("SampleName" = "name")) %>%
      dplyr::mutate(SampleClass = comparisons$class[match(SampleName, comparisons$sample)],
                    SampleLabel = ifelse(SampleClass == "plasmid",SampleName,SampleLabel)) %>%
      dplyr::select(slx_id,screen_type,screen_goal,
                    Flowcell:Q30_bases_samples_percent,
                    virus_batch,plasmid_batch,cas_activity,minimum_split_cell_number,
                    SampleId,SampleName,SampleLabel,SampleClass,
                    dplyr::everything())
    return(qc_metrics)
  },
  error = function(e) stop(paste("merge_bcl2fastq_config_qc: unable to merge bcl2fastq data with config qc data:",e))
  )
}
