#' add_ROC
#'
#' Adds TPR and FPR columns to a data frame containing TP and a score metric (e.g. Bayes factors).
#'
#' @param data A data frame containing a TP and score column.
#' @param score_col A character string naming the column to be used as a score column.
#'
#' @return A data frame containing the following columns: `True_Positive_Rate` and `False_Positive_Rate`.
#' @importFrom dplyr arrange desc mutate filter select sym
#' @importFrom magrittr %<>%
#' @export
add_ROC <- function(data, score_col){
  if(!"TP" %in% colnames(data))
    stop("expecting to find a 'TP' column in 'data'")
  if(!score_col %in% colnames(data))
    stop(paste("expecting to find a",score_col,"column in 'data'"))

  tryCatch({
    sc <- dplyr::sym(score_col)
    data %<>%
      dplyr::arrange(dplyr::desc(!!sc)) %>%
      dplyr::mutate(True_Positive_Rate = cumsum(TP)/sum(TP),
                    False_Positive_Rate = cumsum(!TP)/sum(!TP))
  },
  error = function(e) stop(paste("add_ROC: unable to calculate ROC:",e))
  )
  return(data)
}
