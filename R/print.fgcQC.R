#' print.fgcQC
#'
#' S3 method for generic function `print` for printing summaries of objects of class `fgcQC` to the console.
#'
#' @param x An object of class `fgcQC`.
#' @param ... Other named arguments to be passed to `print`.
#'
#' @return Prints a summary of the object to the console.
#' @export
print.fgcQC <- function(x, ...){
  # S3 method for generic function "print".
  # x is an "fgcQC" object.
  if("plasmid" %in% x$qc_metrics$SampleClass)
    x$qc_metrics <- x$qc_metrics[!x$qc_metrics$SampleClass == "plasmid",]
  qc <- x$qc_metrics
  out <- paste("####  fgcQC summary  ####",
            "\nSLX_ID: ", qc$slx_id[1],
            "\nScreen Type: ", qc$screen_type[1],
            "\nScreen Goal: ", qc$screen_goal[1],
            "\nNumber of Flowcells: ", length(unique(c(unlist(strsplit(qc$Flowcell,";"))))))
  if(qc$screen != "a"){
    out <- paste(out,
                 "\nAUROC for pan-cancer Sanger (Control vs Plasmid): ", round(qc$AUROC.ctrl_plasmid.pan_cancer_Sanger[1],3),
                 "\nAUROC for Hart essentials (Control vs Plasmid): ", round(qc$AUROC.ctrl_plasmid.hart_essential[1],3),
                 "\nAUROC for Moderately negative (Control vs Plasmid): ", round(qc$AUROC.ctrl_plasmid.moderately_negative[1],3),
                 "\nAUROC for Weakly negative (Control vs Plasmid): ", round(qc$AUROC.ctrl_plasmid.weakly_negative[1],3),
                 sep="")
  }
  cat(out)
}
