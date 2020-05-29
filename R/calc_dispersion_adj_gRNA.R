#' calc_dispersion_adj_gRNA
#'
#' Calculates adjusted dispersion estimates at the gRNA level using a shrinkage estimator and returns the mean across all gRNAs for the Control vs Treatment comparison.
#'
#' @param counts A data frame of normalized counts for each sample in the study (samples as columns, gRNAs as rows).
#' @param control_sample_id A character vector naming one or more control samples.
#' @param treat_sample_id A character vector naming one or more treatment samples. If `NULL` then there are no treatment samples.
#'
#' @importFrom DSS newSeqCountSet estNormFactors estDispersion
#' @references Wu, H. et al. 2013. A new shrinkage estimator for dispersion improves differential expression detection in RNA-seq data. Biostatistics 14: 232-243.
#' @export
calc_dispersion_adj_gRNA <- function(counts, control_sample_id, treat_sample_id){
  tryCatch({
    if(!is.null(treat_sample_id)){
      counts <- as.matrix(counts[,c(control_sample_id, treat_sample_id)])
      design <- rep(0,length(c(control_sample_id, treat_sample_id)))
      design[colnames(counts) %in% treat_sample_id] <- 1
      colnames(counts) <- NULL
      seq_data <- DSS::newSeqCountSet(counts, design)
      seq_data <- DSS::estNormFactors(seq_data)
      seq_data <- DSS::estDispersion(seq_data)
      return(data.frame(dispersion_adj_gRNA.treat_ctrl = mean(seq_data@dispersion, na.rm=T)))
    }else{
      return(data.frame(dispersion_adj_gRNA.treat_ctrl = NA))
    }
  },
  error = function(e) stop(paste("calc_dispersion_adj_gRNA: unable to calculate dispersion estimates for treatment-control:",e))
  )
}
