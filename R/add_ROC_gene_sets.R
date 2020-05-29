#' add_bagel_ROC_gene_sets
#'
#' Adds True Positive and False Positive columns to Bagel output for a list of gene sets and calculates AUROC for each gene set.
#'
#' @param bagel_ctrl_pl A data frame containing Bagel Bayes Factors results for Control vs Plasmid.
#' @param bagel_treat_pl A data frame containing Bagel Bayes Factors results for Treatment vs Plasmid. Ignored if `NULL`.
#' @param gene_sets A named list of gene sets.
#'
#' @return A list with the following elements:
#' `bagel_ctrl_pl` - The original data frame with columns `True_Positive_Rate.<gene_set>` and `False_Positive_Rate.<gene_set>` where '<gene_set>' is the name of a gene set.
#' `bagel_treat_pl` - The original data frame with columns `True_Positive_Rate.<gene_set>` and `False_Positive_Rate.<gene_set>` where '<gene_set>' is the name of a gene set. `NULL` if not given.
#' `AUROC` - A data frame containing columns `AUROC.<comp>.<gene_set>` where '<comp>' is either 'ctrl_plasmid' or 'treat_plasmid' and <gene_set>' is the name of a gene set.
#' @importFrom dplyr mutate sym rename select
#' @importFrom magrittr %<>%
#' @importFrom tibble tibble
#' @importFrom rlang :=
#' @export
add_bagel_ROC_gene_sets <- function(bagel_ctrl_pl, bagel_treat_pl, gene_sets){
  ret_auc <- data.frame(dummy = NA)
  noness <- gene_sets$hart_nonessential
  gene_sets$hart_nonessential <- NULL
  for(i in 1:length(gene_sets)){
    tryCatch({
      ## Control vs Plasmid.
      bagel_ctrl_pl %<>%
        dplyr::mutate(TP = GENE %in% gene_sets[[i]]) %>%
        fgcQC::add_ROC("BF") %>%
        dplyr::rename(!!paste("True_Positive_Rate",names(gene_sets)[i],sep=".") := True_Positive_Rate,
                      !!paste("False_Positive_Rate",names(gene_sets)[i],sep=".") := False_Positive_Rate)
      auroc <- bagel_ctrl_pl %>%
        fgcQC::calc_AUC("BF") %>%
        dplyr::rename(!!paste("AUROC.ctrl_plasmid",names(gene_sets)[i],sep=".") := AUC)
      ret_auc <- cbind(ret_auc, auroc)
      if(!is.null(bagel_treat_pl)){
        # Treatment vs Plasmid.
        bagel_treat_pl %<>%
          dplyr::mutate(TP = GENE %in% gene_sets[[i]]) %>%
          fgcQC::add_ROC("BF") %>%
          dplyr::rename(!!paste("True_Positive_Rate",names(gene_sets)[i],sep=".") := True_Positive_Rate,
                        !!paste("False_Positive_Rate",names(gene_sets)[i],sep=".") := False_Positive_Rate)
        auroc <- bagel_treat_pl %>%
          fgcQC::calc_AUC("BF") %>%
          dplyr::rename(!!paste("AUROC.treat_plasmid",names(gene_sets)[i],sep=".") := AUC)
        ret_auc <- cbind(ret_auc, auroc)
      }else{
        auroc <- tibble::tibble(!!paste("AUROC.treat_plasmid",names(gene_sets)[i],sep=".") := NA)
        ret_auc <- cbind(ret_auc, auroc)
      }
    },
    error = function(e) stop(paste("add_bagel_ROC_gene_sets: unable to process gene sets:",e))
    )
  }
  ret_auc %<>%
    dplyr::select(-dummy)
  ret <- list()
  ret$bagel_ctrl_pl <- bagel_ctrl_pl %>% dplyr::select(-TP)
  ret$bagel_treat_pl <- bagel_treat_pl %>% dplyr::select(-TP)
  ret$AUROC <- ret_auc
  return(ret)
}
