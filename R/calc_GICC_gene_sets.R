# Helper function.
.run_gicc <- function(data, group1, group2){
  tryCatch({
    data <- data[,c(group1, group2)]
    samps <- colnames(data)
    samps[samps %in% group1] <- "g1"
    samps[samps %in% group2] <- "g2"
    td <- data.frame(id = samps, t(data), stringsAsFactors = F)
    gicc <- fgcQC::GICC(td)
  },
  error = function(e) stop(paste("GICC: unable to calculate GICC for:",group1,"and",group2,":",e))
  )
  return(gicc$GICC$GICC.1)
}


#' calc_GICC_gene_sets
#'
#' Calculates GICC for normalized counts given a list of gene sets.
#'
#' @param counts A data frame of normalized counts for each sample in the study (samples as columns, gRNAs as rows).
#' @param gene_sets A named list of gene sets.
#' @param plasmid_sample_id A character vector naming one or more plasmid samples.
#' @param control_sample_id A character vector naming one or more control samples.
#' @param treat_sample_id A character vector naming one or more treatment samples. If `NULL` then there are no treatment samples. Defaults to `NULL`.
#'
#' @return A data frame with columns named `GICC.<comp>.<gene_set>` where '<comp' is one of 'ctrl_plasmid', 'treat_plasmid', or 'treat_ctrl', and '<gene_set>' is the name of each gene set provided in the `gene_sets` list.
#' @importFrom dplyr filter do select
#' @importFrom tibble tibble
#' @importFrom rlang :=
#' @export
calc_GICC_gene_sets <- function(counts, gene_sets, plasmid_sample_id, control_sample_id, treat_sample_id = NULL){
  ret <- data.frame(dummy = NA)
  for(i in 1:length(gene_sets)){
    tryCatch({
      # Control vs Plasmid.
      gcol <- sym(paste("GICC","ctrl_plasmid",names(gene_sets)[i],sep="."))
      gicc_cvp <- counts[,c("gene", plasmid_sample_id, control_sample_id)] %>%
        dplyr::filter(gene %in% gene_sets[[i]]) %>%
        dplyr::do(tibble::tibble(!!gcol := .run_gicc(., plasmid_sample_id, control_sample_id)))
      ret <- cbind(ret, gicc_cvp)
      gcol_tp <- sym(paste("GICC","treat_plasmid",names(gene_sets)[i],sep="."))
      gcol_tc <- sym(paste("GICC","treat_ctrl",names(gene_sets)[i],sep="."))
      if(!is.null(treat_sample_id)){
        # Treatment vs Plasmid.
        gicc_tvp <- counts[,c("gene", plasmid_sample_id, treat_sample_id)] %>%
          dplyr::filter(gene %in% gene_sets[[i]]) %>%
          dplyr::do(tibble::tibble(!!gcol_tp := .run_gicc(., plasmid_sample_id, treat_sample_id)))
        # Treatment vs Control.
        gicc_tvc <- counts[,c("gene", treat_sample_id, control_sample_id)] %>%
          dplyr::filter(gene %in% gene_sets[[i]]) %>%
          dplyr::do(tibble::tibble(!!gcol_tc := .run_gicc(., treat_sample_id, control_sample_id)))
        ret <- cbind(ret, gicc_tvp, gicc_tvc)
      }else{
        gicc_tvp <- tibble::tibble(!!gcol_tp := NA)
        gicc_tvc <- tibble::tibble(!!gcol_tc := NA)
        ret <- cbind(ret, gicc_tvp, gicc_tvc)
      }
    },
    error = function(e) stop(paste("calc_GICC_gene_sets: unable to process gene sets:",e))
    )
  }
  return(ret %>% dplyr::select(-dummy))
}
