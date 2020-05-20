.get_auc <- function(scores, TP){
  # Taking all pairs of TP and FP, what is probability TP have a higher score than FP?
  score_fp <- scores[!TP]
  score_tp <- scores[TP]
  o <- outer(score_tp, score_fp, "-")
  return(mean((o>0) + 0.5*(o==0)))
}

#' calc_AUC
#'
#' Calculates AUROC in a data frame containing a 'TP' column.
#'
#' @param data A data frame containing a 'TP' column.
#' @param score_col A character string naming the column to be used as a score column.
#' @param group A character string naming a column by which to group 'data'. If `NULL` no grouping, defaults to `NULL`.
#'
#' @return A data frame.
#' @importFrom dplyr mutate summarise
#' @importFrom magrittr %<>%
#' @export
calc_AUC <- function(data, score_col, group = NULL){
  if(!"TP" %in% colnames(data))
    stop("expecting to find a 'TP' column in 'data'")
  if(!score_col %in% colnames(data))
    stop(paste("expecting to find a",score_col,"column in 'data'"))
  if(!is.null(group)){
    if(!group %in% colnames(data))
      stop(paste("expecting to find a",group,"column"))
    grp <- sym(group)
  }
  sc <- sym(score_col)
  if(!is.null(group)){
    data %<>%
      group_by(!!grp)
  }
  data %<>%
      summarise(AUC = .get_auc(!!sc, TP))
  return(data)
}
