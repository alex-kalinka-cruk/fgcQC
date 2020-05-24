.sanity_check_b2f_out <- function(data){

}


#' extract_b2f_json
#'
#' Extracts Sequencing QC metrics from a `bcl2fastq2` 'Stats.json' file.
#'
#' @param path A valid path to a 'Stats.json' file produced by `bcl2fastq2`.
#'
#' @return A list containing the following elements:
#' `flowcell` - a data frame containing flowcell-level data.
#' `run_info` - a data frame containing read lengths.
#' `summary` - a data frame containing run-level QC metrics.
#' `samples` - a data frame containing sample-level QC metrics.
#'
#' @importFrom jsonlite fromJSON
#' @importFrom dplyr as_tibble select mutate
#' @importFrom tidyr unnest
#' @importFrom magrittr %<>%
#' @export
#' @examples \dontrun{jj <- extract_b2f_json('Stats.json')}
extract_b2f_json <- function(path){
  ret <- list()
  tryCatch({
    all <- jsonlite::fromJSON(path) %>%
      dplyr::as_tibble()
    # Sanity check input.
    # Flowcell info.
    flowcell <- all[,1:3]
    # RunInfo (read config).
    run_info <- as.data.frame(all$ReadInfosForLanes$ReadInfos) %>%
      dplyr::mutate(RunId = flowcell$RunId)
    # Flowcell summary stats.
    summ <- all$ConversionResults %>%
      dplyr::select(-DemuxResults, -Undetermined) %>%
      dplyr::mutate(RunId = flowcell$RunId)
    # Sample summary stats.
    samps <- all$ConversionResults$DemuxResults[[1]] %>%
      dplyr::mutate(RunId = flowcell$RunId) %>%
      tidyr::unnest(c(IndexMetrics, ReadMetrics),
                    names_repair="universal") %>%
      dplyr::mutate(MismatchCounts = MismatchCounts[,1],
                    Index_OneBaseMismatch_percent = 100*(NumberReads-MismatchCounts)/NumberReads,
                    Q30_bases_percent = 100*YieldQ30/Yield...8,
                    Average_base_quality = QualityScoreSum/Yield...8,
                    Trimmed_bases_percent = 100*TrimmedBases/Yield...8,
                    Sample_Representation = 100*NumberReads/summ$TotalClustersPF)
    not_demux_count <- all$ConversionResults$Undetermined$NumberReads
    # Flowcell-level summary.
    summ %<>%
      dplyr::mutate(ReadsPF_percent = 100*TotalClustersPF/TotalClustersRaw,
                    Non_Demultiplexed_Reads_percent = 100*not_demux_count/summ$TotalClustersPF,
                    Q30_bases_samples_percent = 100*sum(samps$YieldQ30)/sum(samps$Yield...8))
    ret$flowcell <- flowcell
    ret$run_info <- run_info
    ret$summary <- summ
    ret$samples <- samps
  },
  error = function(e) stop(paste("unable to extract bcl2fastq2 'Stats.json' output:",e))
  )
  return(ret)
}
