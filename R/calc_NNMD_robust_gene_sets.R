#' calc_NNMD_robust_gene_sets
#'
#' Calculates a Robust equivalent of the Null-Normalized Mean Difference (NNMD) for sets of genes (based on median difference divided by median absolute deviation of reference set).
#'
#' @param lfc_ctrl_pl A data frame of gRNA log fold change data for Control vs Plasmid.
#' @param lfc_treat_pl A data frame of gRNA log fold change data for Treatment vs Plasmid. If `NULL`, then ignored.
#' @param gene_sets A named list of gene sets.
#'
#' @return A data frame with columns named `NNMD_robust.<comp>.<gene_set>` where '<comp' is one of 'ctrl_plasmid' or 'treat_plasmid', and '<gene_set>' is the name of each gene set provided in the `gene_sets` list.
#' @author Alex T. Kalinka, \email{alex.kalinka@@cancer.org.uk}
#' @importFrom dplyr filter do select sym
#' @importFrom tibble tibble
#' @importFrom rlang :=
#' @export
calc_NNMD_robust_gene_sets <- function(lfc_ctrl_pl, lfc_treat_pl, gene_sets){
  ret <- data.frame(dummy = NA)
  noness <- gene_sets$hart_nonessential
  gene_sets$hart_nonessential <- NULL
  for(i in 1:length(gene_sets)){
    tryCatch({
      # Control vs Plasmid.
      ncol <- dplyr::sym(paste("NNMD_robust","ctrl_plasmid",names(gene_sets)[i],sep="."))
      nnmd_cvp <- lfc_ctrl_pl %>%
        dplyr::do(tibble::tibble(!!ncol := fgcQC::NNMD_robust(., "gene", "log2FC", gene_sets[[i]], noness)))
      ret <- cbind(ret, nnmd_cvp)
      ncol_tp <- dplyr::sym(paste("NNMD_robust","treat_plasmid",names(gene_sets)[i],sep="."))
      if(!is.null(lfc_treat_pl)){
        # Treatment vs Plasmid.
        nnmd_tvp <- lfc_treat_pl %>%
          dplyr::do(tibble::tibble(!!ncol := fgcQC::NNMD_robust(., "gene", "log2FC", gene_sets[[i]], noness)))
      }else{
        nnmd_tvp <- tibble::tibble(!!ncol_tp := NA)
        ret <- cbind(ret, nnmd_tvp)
      }
    },
    error = function(e) stop(paste("calc_NNMD_robust_gene_sets: unable to process gene sets:",e))
    )
  }
  return(ret %>% dplyr::select(-dummy))
}
