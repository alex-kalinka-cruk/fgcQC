#' assemble_logfc_QC
#'
#' Assembles logFC-based QC metrics.
#'
#' @param qc_metrics A data frame containing only sequencing-based QC metrics (rows as samples, QC metrics as columns).
#' @param screen_goal A character string naming the screen goal.
#' @param lfc.ctrl_pl A data frame of gRNA logFC data for Control vs Plasmid (gRNAs as rows).
#' @param lfc.treat_pl A data frame of gRNA logFC data for Treatment vs Plasmid (gRNAs as rows). Ignored, if `NULL`.
#' @param lfc.ctrl_pl.repl A data frame of gRNA logFC data for replicate-based Control vs Plasmid (gRNAs as rows).
#' @param lfc.treat_pl.repl A data frame of gRNA logFC data for replicate-based Treatment vs Plasmid (gRNAs as rows).
#' @param library A data frame containing the library file in which the first column gives the sgRNA sequence and the second column gives the sgRNA ID.
#'
#' @return A modified `qc_metrics` data frame containing all logFC-based QC metrics (both experiment and sample-level).
#' @author Alex T. Kalinka, \email{alex.kalinka@@cancer.org.uk}
#' @importFrom magrittr %<>%
#' @importFrom tibble add_column
#' @export
assemble_logfc_QC <- function(qc_metrics, screen_goal, lfc.ctrl_pl, lfc.treat_pl, lfc.ctrl_pl.repl, lfc.treat_pl.repl, library){
  tryCatch({
    qc_metrics %<>%
      ## NNMD.
      tibble::add_column(fgcQC::calc_NNMD_gene_sets(lfc.ctrl_pl, lfc.treat_pl, fgcQC::crispr_gene_sets$essential),
                         .before = "SampleId") %>%
      ## NNMD_robust.
      tibble::add_column(fgcQC::calc_NNMD_robust_gene_sets(lfc.ctrl_pl, lfc.treat_pl, fgcQC::crispr_gene_sets$essential),
                         .before = "SampleId") %>%
      ## distance correlation between gRNA logFC and GC content.
      tibble::add_column(fgcQC::calc_dcorr_GC_content_logfc(lfc.ctrl_pl, library, "ctrl_plasmid"),
                         .before = "SampleId") %>%
      ## inefficient gRNA (GCC and TT) logFC ratios.
      tibble::add_column(fgcQC::calc_inefficient_gRNA_logfc(lfc.ctrl_pl, library, "ctrl_plasmid"),
                         .before = "SampleId") %>%
      ## replicate logFC correlations.
      tibble::add_column(fgcQC::calc_replicate_logfc_corr(lfc.ctrl_pl.repl, "ctrl_plasmid"),
                         .before = "SampleId") %>%
      tibble::add_column(fgcQC::calc_replicate_logfc_corr(lfc.treat_pl.repl, "treat_plasmid"),
                         .before = "SampleId")

    if(screen_goal != "lethality"){
      qc_metrics %<>%
        ## distance correlation between gRNA logFC and GC content.
        tibble::add_column(fgcQC::calc_dcorr_GC_content_logfc(lfc.treat_pl, library, "treat_plasmid"),
                           .before = "SampleId") %>%
        ## inefficient gRNA (GCC and TT) logFC ratios.
        tibble::add_column(fgcQC::calc_inefficient_gRNA_logfc(lfc.treat_pl, library, "treat_plasmid"),
                           .before = "SampleId")
    }else{
      qc_metrics %<>%
        ## distance correlation between gRNA logFC and GC content.
        tibble::add_column(data.frame(distcorr_GC_content_logfc.treat_plasmid = NA),
                           .before = "SampleId") %>%
        ## inefficient gRNA (GCC and TT) logFC ratios.
        tibble::add_column(data.frame(log2FC_GCC_diff.treat_plasmid = NA, log2FC_TT_diff.treat_plasmid = NA),
                           .before = "SampleId")
    }

    # Add NNMD for Control gRNAs (non-targeting guides), if they exist.
    if(any(grepl("^Control",library$sgRNA))){
      ctrl_genes <- unique(library$gene[grepl("^Control",library$gene)])
      qc_metrics %<>%
        ## NNMD.
        tibble::add_column(fgcQC::calc_NNMD_gene_sets(lfc.ctrl_pl, lfc.treat_pl, list(control_guides = ctrl_genes)),
                           .before = "SampleId") %>%
        ## NNMD_robust.
        tibble::add_column(fgcQC::calc_NNMD_robust_gene_sets(lfc.ctrl_pl, lfc.treat_pl, list(control_guides = ctrl_genes)),
                           .before = "SampleId")
    }else{
      qc_metrics %<>%
        tibble::add_column(data.frame(NNMD.ctrl_plasmid.control_guides = NA, NNMD.treat_plasmid.control_guides = NA,
                                      NNMD_robust.ctrl_plasmid.control_guides = NA, NNMD_robust.treat_plasmid.control_guides = NA),
                           .before = "SampleId")
    }
  },
  error = function(e) stop(paste("assemble_logfc_QC: unable to build logFC QC metrics:",e))
  )
  return(qc_metrics)
}
