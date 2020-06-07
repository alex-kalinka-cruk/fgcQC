# Helper function.
.sanity_check_b2f_out <- function(b2f_list){
  exp_names <- c("Flowcell","RunNumber","RunId","ReadInfosForLanes","ConversionResults","UnknownBarcodes")
  if(!setequal(names(b2f_list), exp_names))
    stop(paste(".sanity_check_b2f_out: bcl2fastq2 expected sections:",exp_names,"found instead:",names(b2f_list)))
}


#' extract_b2f_json
#'
#' Extracts Sequencing QC metrics from a `bcl2fastq2` 'Stats.json' file.
#'
#' @param path A valid path to a 'Stats.json' file produced by `bcl2fastq2`.
#'
#' @return A data frame containing flowcell and sample-level sequencing metrics.
#' @author Alex T. Kalinka, \email{alex.kalinka@@cancer.org.uk}
#'
#' @importFrom jsonlite fromJSON
#' @importFrom dplyr as_tibble select mutate inner_join right_join rename
#' @importFrom tidyr unnest
#' @importFrom magrittr %<>%
#' @export
extract_b2f_json <- function(path){
  tryCatch({
    all <- jsonlite::fromJSON(path)
    # Sanity check input.
    .sanity_check_b2f_out(all)
    all %<>%
      dplyr::as_tibble()
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
                    Sample_Representation = 100*NumberReads/summ$TotalClustersPF,
                    SampleId = gsub("_","-",SampleId))
    not_demux_count <- all$ConversionResults$Undetermined$NumberReads
    # Flowcell-level summary.
    summ %<>%
      dplyr::mutate(ReadsPF_percent = 100*TotalClustersPF/TotalClustersRaw,
                    Non_Demultiplexed_Reads_percent = 100*not_demux_count/summ$TotalClustersPF,
                    Q30_bases_samples_percent = 100*sum(samps$YieldQ30)/sum(samps$Yield...8))

    flowcell %<>%
      dplyr::mutate(RunDate = gsub("^(..)(..)(..)$","20\\1-\\2-\\3",unlist(strsplit(RunId,"_"))[1]),
                    read_lengths = paste(sort(run_info$NumCycles,decreasing = T),collapse="-")) %>%
      dplyr::inner_join(summ, by = "RunId") %>%
      dplyr::right_join(samps, by = "RunId") %>%
      dplyr::rename(Yield_bases_total = Yield, Yield_bases_sample = Yield...8)
  },
  error = function(e) stop(paste("unable to extract bcl2fastq2 'Stats.json' output:",e))
  )
  return(flowcell)
}
