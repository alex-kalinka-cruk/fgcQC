#' summarise_samples_bcl2fastq
#'
#' Summarises sequencing metrics for samples sequenced across multiple lanes and/or flowcells. Minimum values are taken for numeric metrics and count metrics are summed.
#'
#' @param data A data frame of sequencing metrics derived ultimately from `fgcQC::extract_b2f_json`.
#'
#' @return A summarised data frame with one row per sample.
#' @importFrom dplyr group_by ungroup summarise select
#' @importFrom magrittr %<>%
#' @export
summarise_samples_bcl2fastq <- function(data){
  tryCatch({
    data %<>%
      dplyr::group_by(SampleId) %>%
      dplyr::summarise(Flowcell = paste(Flowcell,collapse=";"), RunId = paste(RunId,collapse=";"),
                       RunDate = paste(RunDate,collapse=";"), read_lengths = read_lengths[1],
                       LaneNumber = paste(LaneNumber,collapse=";"),
                       ReadsPF_percent = min(ReadsPF_percent, na.rm=T),
                       Non_Demultiplexed_Reads_percent = min(Non_Demultiplexed_Reads_percent, na.rm=T),
                       Q30_bases_samples_percent = min(Q30_bases_samples_percent, na.rm=T),
                       SampleName = SampleName[1], IndexSequence = IndexSequence[1],
                       NumberReads = sum(NumberReads, na.rm=T),
                       Index_OneBaseMismatch_percent = min(Index_OneBaseMismatch_percent, na.rm=T),
                       Q30_bases_percent = min(Q30_bases_percent, na.rm=T),
                       Average_base_quality = min(Average_base_quality, na.rm=T),
                       Sample_Representation = min(Sample_Representation, na.rm=T)) %>%
      dplyr::ungroup() %>%
      dplyr::select(Flowcell:Q30_bases_samples_percent, SampleId, SampleName:Sample_Representation)
  },
  error = function(e) stop(paste("unable to summarise bcl2fastq2 data:",e))
  )
  return(data)
}
