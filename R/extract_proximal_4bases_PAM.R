#' extract_proximal_4bases_PAM
#'
#' Extracts the proximal 4 bases upstream of the PAM in a sgRNA library data frame into "GCC", "TT", or "Other" (Graf et al. 2019).
#'
#' @param library A data frame containing the library file in which the first column gives the sgRNA sequence and the second column gives the sgRNA ID.
#'
#' @return A data frame with a column named `proximal_4bases_pam`.
#' @author Alex T. Kalinka, \email{alex.kalinka@@cancer.org.uk}
#' @importFrom dplyr mutate rowwise ungroup case_when
#' @importFrom magrittr %<>%
#' @references Graf, R. et al. (2019) sgRNA sequence motifs blocking efficient CRISPR/Cas9-mediated gene editing. Cell Rep 26: 1098-1103.
#' @export
extract_proximal_4bases_PAM <- function(library){
  tryCatch({
    guide_len <- nchar(library$V1[1])
    library %<>%
      dplyr::rowwise() %>%
      dplyr::mutate(pam_upstream4 = substr(V1, guide_len-3, guide_len),
      proximal_4bases_pam = dplyr::case_when(grepl("TT",pam_upstream4) ~ "TT",
                                      grepl("GCC",pam_upstream4) ~ "GCC",
                                      TRUE ~ "Other")) %>%
      dplyr::ungroup()
  },
  error = function(e) stop(paste("extract_proximal_4bases_PAM: unable to extract the proximal 4 bases of the PAM:",e))
  )
  return(library)
}
