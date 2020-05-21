.geom_mean <- function(...){
  x <- c(unlist(list(...)))
  return(exp(sum(log(x[x > 0]), na.rm=T) / length(x)))
}


#' normalize_library_depth_median_ratio
#'
#' Normalizes library depth from a counts file using the median-ratio method [1].
#'
#' @param counts A data frame of counts for each sample in the study (samples as columns, gRNAs as rows).
#'
#' @return A data frame in which the sample columns are library-depth normalized.
#' @importFrom dplyr mutate summarise select
#' @importFrom purrr pmap_dbl
#' @export
normalize_library_depth_median_ratio <- function(counts){
  tryCatch({
    size_factors <- counts %>%
      dplyr::mutate(geom_mean = purrr::pmap_dbl(
        .l = dplyr::select(., -sgRNA, -gene),
        .f = .geom_mean
      )) %>%
      dplyr::summarise()
  },
  error = function(e) stop(paste("unable to calculate median ratio normalization:",e))
  )
  return(size_factors)
}
