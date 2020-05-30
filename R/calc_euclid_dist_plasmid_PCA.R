# Helper function.
.euclid_dist_2D <- function(q1, q2, p1, p2){
  return(sqrt((q1-p1)^2+(q2-p2)^2))
}


#' calc_euclid_dist_plasmid_PCA
#'
#' Calculates a 2D-Euclidean distance from a set of given samples to the plasmid sample for the first 2 PCs of a PCA derived from the normalized counts. The minimum distance is returned for each comparison.
#'
#' @param counts A data frame of normalized counts for each sample in the study (samples as columns, gRNAs as rows).
#' @param plasmid_sample_id A character vector naming the plasmid sample(s).
#' @param control_sample_id A character vector naming the control sample(s).
#' @param treat_sample_id A character vector naming the treatment sample(s). Ignored if `NULL`.
#'
#' @return A data frame with the following column(s): `euclid_dist_plasmid_PCA.<comp>` where '<comp>' is one of 'control' or 'treat'.
#' @importFrom dplyr mutate filter summarise
#' @importFrom tibble rownames_to_column
#' @export
calc_euclid_dist_plasmid_PCA <- function(counts, plasmid_sample_id, control_sample_id, treat_sample_id){
  tryCatch({
    # Control samples.
    pca.ctrl <- as.data.frame(scale(prcomp(as.matrix(t(scale(counts[,c(plasmid_sample_id, control_sample_id)],scale=F))))$x[,1:2])) %>%
      tibble::rownames_to_column("sample_id") %>%
      dplyr::mutate(sample_type = ifelse(sample_id %in% plasmid_sample_id,"plasmid","control"))
    plasmid_median <- pca.ctrl %>%
      dplyr::filter(sample_type == "plasmid") %>%
      dplyr::summarise(PC1 = median(PC1, na.rm=T), PC2 = median(PC2, na.rm=T))
    ret <- pca.ctrl %>%
      dplyr::filter(sample_type == "control") %>%
      dplyr::mutate(euclid_dist_2d = .euclid_dist_2D(plasmid_median$PC1, plasmid_median$PC2,
                                                     PC1, PC2)) %>%
      dplyr::summarise(euclid_dist_plasmid_PCA.control = min(euclid_dist_2d))

    if(!is.null(treat_sample_id)){
      # Treatment samples.
      pca.treat <- as.data.frame(scale(prcomp(as.matrix(t(scale(counts[,c(plasmid_sample_id, treat_sample_id)],scale=F))))$x[,1:2])) %>%
        tibble::rownames_to_column("sample_id") %>%
        dplyr::mutate(sample_type = ifelse(sample_id %in% plasmid_sample_id,"plasmid","treat"))
      plasmid_median <- pca.treat %>%
        dplyr::filter(sample_type == "plasmid") %>%
        dplyr::summarise(PC1 = median(PC1, na.rm=T), PC2 = median(PC2, na.rm=T))
      ret <- cbind(ret, pca.treat %>%
        dplyr::filter(sample_type == "treat") %>%
        dplyr::mutate(euclid_dist_2d = .euclid_dist_2D(plasmid_median$PC1, plasmid_median$PC2,
                                                       PC1, PC2)) %>%
        dplyr::summarise(euclid_dist_plasmid_PCA.treat = min(euclid_dist_2d))
      )
    }
  return(ret)
  },
  error = function(e) stop(paste("calc_euclid_dist_plasmid_PCA: unable to calculate 2D Euclidean distance to plasmid:",e))
  )
}
