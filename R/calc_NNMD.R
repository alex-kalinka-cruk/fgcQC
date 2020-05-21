#' calc_NNMD
#'
#' Calculates Null-Normalized Mean Difference (NNMD; also known as Glass's delta) for the log fold change distributions of essential and non-essential genes - defined as the mean difference divided by the standard deviation of non-essential genes (the reference set).
#'
#' @param data A data fram containing log fold change data (normalized in terms of library depth and any potential copy-number effects) for a set of genes.
#' @param gene_col A character string indicating which column contains gene names.
#' @param lfc_col A character string indicating which column contains $log_2$ fold change data for each gene.
#' @param essential_genes A vector of essential gene names.
#' @param nonessential_genes A vector of non-essential gene names.
#'
#' @return A single numeric value.
#' @importFrom dplyr filter select mutate summarise
#' @importFrom magrittr %<>%
#' @export
calc_NNMD <- function(data, gene_col, lfc_col, essential_genes, nonessential_genes){
  if(!gene_col %in% colnames(data))
    stop(paste("expecting to find a",gene_col,"column in 'data'"))
  if(!lfc_col %in% colnames(data))
    stop(paste("expecting to find a",lfc_col,"column in 'data'"))
  gc <- sym(gene_col)
  lfc <- sym(lfc_col)
  data %<>%
    dplyr::filter(!!gc %in% c(essential_genes,nonessential_genes) & ! abs(!!lfc) == Inf) %>%
    dplyr::mutate(TP = !!gc %in% essential_genes, Lfc = !!lfc) %>%
    dplyr::summarise(NNMD.lfc = (mean(Lfc[TP], na.rm=T) - mean(Lfc[!TP], na.rm=T))/sd(Lfc[!TP], na.rm=T))
  return(data$NNMD.lfc)
}
