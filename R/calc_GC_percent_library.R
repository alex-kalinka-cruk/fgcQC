#' calc_GC_percent_library
#'
#' Calculates the GC percent of sgRNAs in a given library.
#'
#' @param library A data frame containing the library file in which the first column gives the sgRNA sequence and the second column gives the sgRNA ID.
#'
#' @return A data frame with a column named `GC_percent`
#' @importFrom dplyr mutate rowwise ungroup filter
#' @importFrom magrittr %<>%
calc_GC_percent_library <- function(library){
  tryCatch({
    library %<>%
      dplyr::rowwise() %>%
      dplyr::mutate(GC_percent = 100*(sum(unlist(strsplit(V1,""))=="G") + sum(unlist(strsplit(V1,""))=="C"))/guide_len) %>%
      dplyr::ungroup()
    },
    error = function(e) stop(paste("unable to calculate GC percent for library:",e))
  )
  return(library)
}
