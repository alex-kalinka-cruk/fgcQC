#' calc_mahalanobis_dist_plasmid
#'
#' Calculates the Mahalanobis distance ratio for control vs plasmid and treatment vs plasmid.
#'
#' @param counts A data frame of normalized counts for each sample in the study (samples as columns, gRNAs as rows).
#' @param plasmid_sample_id A character vector naming the plasmid sample(s).
#' @param control_sample_id A character vector naming the control sample(s).
#' @param treat_sample_id A character vector naming the treatment sample(s). Ignored if `NULL`.
#'
#' @return A data frame with the following column(s): `mahalanobis_dist_ratio.<comp>` where '<comp>' is one of 'ctrl_plasmid' or 'treat_plasmid'.
#' @author Alex T. Kalinka, \email{alex.kalinka@@cancer.org.uk}
#' @importFrom dplyr filter
#' @importFrom tibble rownames_to_column
#' @export
calc_mahalanobis_dist_plasmid <- function(counts, plasmid_sample_id, control_sample_id, treat_sample_id){
  tryCatch({
    counts.ctrl <- as.matrix(t(counts[,c(plasmid_sample_id, control_sample_id)]))
    counts.ctrl <- counts.ctrl[,which(apply(counts.ctrl, 2, var) != 0)]
    mdist <- as.data.frame(as.matrix(dist(scale(prcomp(counts.ctrl)$x))))
    # Control samples.
    dist.ctrl <- mdist %>%
      tibble::rownames_to_column("sample_id") %>%
      dplyr::filter(sample_id %in% control_sample_id)
    plasm_min <- min(c(unlist(dist.ctrl[,plasmid_sample_id])), na.rm = T)
    ctrl_max <- max(c(unlist(dist.ctrl[,control_sample_id])), na.rm = T)
    ret <- data.frame(mahalanobis_dist_ratio.ctrl_plasmid = ctrl_max/plasm_min)

    if(!is.null(treat_sample_id)){
      counts.treat <- as.matrix(t(counts[,c(plasmid_sample_id, treat_sample_id)]))
      counts.treat <- counts.treat[,which(apply(counts.treat, 2, var) != 0)]
      mdist <- as.data.frame(as.matrix(dist(scale(prcomp(counts.treat)$x))))
      # Treatment samples.
      dist.treat <- mdist %>%
        tibble::rownames_to_column("sample_id") %>%
        dplyr::filter(sample_id %in% treat_sample_id)
      plasm_min <- min(c(unlist(dist.treat[,plasmid_sample_id])), na.rm = T)
      treat_max <- max(c(unlist(dist.treat[,treat_sample_id])), na.rm = T)
      ret <- cbind(ret, data.frame(mahalanobis_dist_ratio.treat_plasmid = treat_max/plasm_min))
    }
  return(ret)
  },
  error = function(e) stop(paste("calc_mahalanobis_dist_plasmid: unable to calculate Mahalanobis distance ratio to plasmid:",e))
  )
}
