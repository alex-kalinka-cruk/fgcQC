#' assemble_counts_QC
#'
#' Assembles counts-based QC metrics.
#'
#' @param qc_metrics A data frame containing only sequencing-based QC metrics (rows as samples, QC metrics as columns).
#' @param counts A data frame of counts for each sample in the study (samples as columns, gRNAs as rows).
#' @param counts_norm A data frame of normalized counts for each sample in the study (samples as columns, gRNAs as rows).
#' @param library A data frame containing the library file in which the first column gives the sgRNA sequence and the second column gives the sgRNA ID.
#' @param plasmid_samp A character vector naming one or more plasmid samples.
#' @param control_samp A character vector naming one or more control samples.
#' @param treat_samp A character vector naming one or more treatment samples. If `NULL` then there are no treatment samples.
#'
#' @return A modified `qc_metrics` data frame containing all counts-based QC metrics (both experiment and sample-level).
#' @author Alex T. Kalinka, \email{alex.kalinka@@cancer.org.uk}
#' @importFrom dplyr left_join
#' @importFrom tibble add_column
#' @importFrom magrittr %<>%
#' @export
assemble_counts_QC <- function(qc_metrics, counts, counts_norm, library,
                               plasmid_samp, control_samp, treat_samp){
  tryCatch({
    qc_metrics %<>%
      ### QC for count data.
      ## Zero and low count plasmid gRNAs.
      tibble::add_column(fgcQC::calc_low_zero_count_plasmid_gRNAs(counts, plasmid_samp),
                         .before = "SampleId") %>%
      ## Percent reads matching gRNAs.
      fgcQC::calc_percent_matching_gRNAs(counts) %>%
      ## Gini coefficients.
      dplyr::left_join(fgcQC::calc_gini_coefficient_counts(counts), by = "SampleName") %>%
      ## distance correlation between gRNA counts and GC content.
      dplyr::left_join(fgcQC::calc_dcorr_GC_content_counts(counts_norm, library), by = "SampleName") %>%
      ## inefficient gRNA (GCC and TT) count ratios.
      dplyr::left_join(fgcQC::calc_inefficient_gRNA_count_ratios(counts_norm, library), by = "SampleName") %>%
      ## GICC.
      tibble::add_column(fgcQC::calc_GICC_gene_sets(counts_norm, fgcQC::crispr_gene_sets$essential,
                                                    plasmid_samp, control_samp, treat_samp),
                         .before = "SampleId") %>%
      ## dispersion shrinkage estimator for Treat-Ctrl.
      tibble::add_column(fgcQC::calc_dispersion_adj_gRNA(counts, control_samp, treat_samp),
                         .before = "SampleId") %>%
      ## Mahalanobis distance to plasmid.
      tibble::add_column(fgcQC::calc_mahalanobis_dist_plasmid(counts_norm, plasmid_samp,
                                                              control_samp, treat_samp),
                         .before = "SampleId")
  },
  error = function(e) stop(paste("assemble_counts_QC: unable to build counts-based QC metrics:",e))
  )
  return(qc_metrics)
}
