#' calc_gini_coefficient_counts
#'
#' Calculates the Gini coefficient of dispersion for sample gRNA count distributions.
#'
#' @param counts A data frame of counts for each sample in the study (samples as columns, gRNAs as rows).
#'
#' @return A data frame containing a `SampleName` column, and a column named `gini_coefficient_counts`.
#' @importFrom dplyr mutate
#' @importFrom ineq ineq
#' @export
calc_gini_coefficient_counts <- function(counts){
  tryCatch({
    gc <- apply(counts[,3:ncol(counts)], 2, ineq, type = "Gini")
    ret <- data.frame(SampleName = names(gc),
                      gini_coefficient_counts = gc)
    rownames(ret) <- NULL
  },
  error = function(e) stop(paste("unable to calculate Gini coefficients:",e))
  )
  return(ret)
}
