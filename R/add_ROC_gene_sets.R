#' add_bagel_ROC_gene_sets
#'
#' Adds True Positive and False Positive columns to Bagel output for a list of gene sets and calculates AUROC for each gene set: all relative to the Hart non-essential genes.
#'
#' @param bagel_ctrl_pl A data frame containing Bagel Bayes Factors results for Control vs Plasmid.
#' @param bagel_treat_pl A data frame containing Bagel Bayes Factors results for Treatment vs Plasmid. Ignored if `NULL`.
#' @param gene_sets A named list of gene sets.
#'
#' @return A list with the following elements:
#' `bagel_ctrl_plasmid` - The original data frame with columns `gene_set`, `True_Positive_Rate.<gene_set>`, `False_Positive_Rate.<gene_set>` where '<gene_set>' is the name of a gene set.
#' `bagel_treat_plasmid` - The original data frame with columns `gene_set`, `True_Positive_Rate.<gene_set>` and `False_Positive_Rate.<gene_set>` where '<gene_set>' is the name of a gene set. `NULL` if not given.
#' `AUROC` - A data frame containing columns `AUROC.<comp>.<gene_set>` where '<comp>' is either 'ctrl_plasmid' or 'treat_plasmid' and <gene_set>' is the name of a gene set.
#' @author Alex T. Kalinka, \email{alex.kalinka@@cancer.org.uk}
#' @importFrom dplyr mutate sym rename select filter
#' @importFrom magrittr %<>%
#' @importFrom tibble tibble
#' @importFrom rlang :=
#' @export
add_bagel_ROC_gene_sets <- function(bagel_ctrl_pl, bagel_treat_pl, gene_sets){
  ret_cp <- ret_tp <- NULL
  ret_auc <- data.frame(dummy = NA)
  noness <- gene_sets$hart_nonessential
  gene_sets$hart_nonessential <- NULL
  for(i in 1:length(gene_sets)){
    tryCatch({
      ## Control vs Plasmid.
      td <- bagel_ctrl_pl %>%
        dplyr::filter(GENE %in% c(gene_sets[[i]], noness)) %>%
        dplyr::mutate(TP = GENE %in% gene_sets[[i]],
                      gene_set = names(gene_sets)[i]) %>%
        fgcQC::add_ROC("BF")
      ret_cp <- rbind(ret_cp, td %>% dplyr::select(-TP))
      auroc <- td %>%
        fgcQC::calc_AUC("BF") %>%
        dplyr::rename(!!paste("AUROC.ctrl_plasmid",names(gene_sets)[i],sep=".") := AUC)
      ret_auc <- cbind(ret_auc, auroc)
      if(!is.null(bagel_treat_pl)){
        # Treatment vs Plasmid.
        td <- bagel_treat_pl %>%
          dplyr::filter(GENE %in% c(gene_sets[[i]], noness)) %>%
          dplyr::mutate(TP = GENE %in% gene_sets[[i]],
                        gene_set = names(gene_sets)[i]) %>%
          fgcQC::add_ROC("BF")
        ret_tp <- rbind(ret_tp, td %>% dplyr::select(-TP))
        auroc <- td %>%
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
  ret$bagel_ctrl_plasmid <- ret_cp
  ret$bagel_treat_plasmid <- ret_tp
  ret$AUROC <- ret_auc
  return(ret)
}
