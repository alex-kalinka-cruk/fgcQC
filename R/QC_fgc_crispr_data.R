# Helper functions.
.print_progress <- function(msg){
  msg <- paste(date(), "*** fgcQC:",msg,"\n")
  cat(msg)
}


.file_exists <- function(file){
  if(!file.exists(file))
    stop(paste(".check_qc_inputs: unable to find",file))
}


.check_qc_inputs <- function(analysis_config, combined_counts, bagel_ctrl_plasmid, bcl2fastq, library, comparison,
                             output, output_R_object, norm_method, bagel_treat_plasmid){
  .file_exists(analysis_config)
  .file_exists(combined_counts)
  if(is.null(library))
    stop(".check_qc_inputs: a path to a valid cleanr library file is required")
  .file_exists(library)
  if(is.null(bagel_ctrl_plasmid)){
    warning(".check_qc_inputs: no 'bagel_ctrl_plasmid' data available")
  }else{
    .file_exists(bagel_ctrl_plasmid)
  }

  if(!is.null(bagel_treat_plasmid))
    .file_exists(bagel_treat_plasmid)
  if(!is.character(bcl2fastq) || is.null(bcl2fastq))
    stop(".check_qc_inputs: 'bcl2fastq' must be a comma-separated character string pointing to one or more bcl2fastq2 output JSON files")
  b2f <- unlist(strsplit(bcl2fastq,","))
  sapply(b2f, .file_exists)
  if(!is.character(comparison))
    stop(".check_qc_inputs: 'comparison_name' must be a character string naming a comparison in the analysis config JSON file")
  if(!is.character(output) & !is.null(output))
    stop(".check_qc_inputs: 'output' must be a character string naming an output file for QC results or NULL")
  if(!is.character(output_R_object) & !is.null(output_R_object))
    stop(".check_qc_inputs: 'output_R_object' must be a character string naming an output file for the QC R object or NULL")
  if(!norm_method %in% c("median_ratio","relative"))
    stop(".check_qc_inputs: 'norm_method' should be one of 'median_ratio' or 'relative'")
}


.read_bcl2fastq <- function(bcl2fastq){
  b2f <- unlist(strsplit(bcl2fastq,","))
  ret <- NULL
  for(i in 1:length(b2f)){
    tryCatch(
      ret <- rbind(ret, fgcQC::extract_b2f_json(b2f[i])),
      error = function(e) stop(paste(".read_bcl2fastq: problem reading bcl2fastq2 data:",e))
    )
  }
  return(ret)
}


.check_comparison <- function(comparisons, comparison_name, analysis_config){
  if(!comparison_name %in% comparisons$comparison)
    stop(paste(".check_comparison: unable to find comparison",comparison_name,"in",analysis_config))
  if(!"plasmid" %in% comparisons$class)
    stop(paste(".check_comparison: unable to find any 'plasmid' samples in",comparison_name))
  if(!"control" %in% comparisons$class)
    stop(paste(".check_comparison: unable to find any 'control' samples in",comparison_name))
  if(comparisons$goal[1] != "lethality"){
    if(!"treatment" %in% comparisons$class)
      stop(paste(".check_comparison: unable to find any 'treatment' samples in",comparison_name))
  }
}


.prep_bcl2fastq <- function(bcl2fastq, comparisons, json_list){
  tryCatch({
    b2f <- .read_bcl2fastq(bcl2fastq) %>%
      fgcQC::summarise_samples_bcl2fastq() %>%
      dplyr::mutate(screen_type = comparisons$type[1],
                    screen_goal = comparisons$goal[1],
                    SampleName = json_list$samples$name[match(SampleId, json_list$samples$indexes)],
                    SampleLabel = json_list$samples$label[match(SampleId, json_list$samples$indexes)]) %>%
      # Keep only samples in the given comparison.
      dplyr::filter(SampleName %in% comparisons$sample)
    return(b2f)
  },
  error = function(e) stop(paste(".prep_bcl2fastq: unable to prep bcl2fastq data:",e))
  )
}


.make_dummy_bagel <- function(gene_set, type){
  tryCatch({
    gene_set$hart_nonessential <- NULL
    ret <- data.frame(dummy = NA)
    for(genes in names(gene_set)){
      ret %<>%
        dplyr::mutate(!!paste(type,"ctrl_plasmid",genes,sep=".") := NA,
                      !!paste(type,"treat_plasmid",genes,sep=".") := NA)
      if(type == "AUPrRc"){
        ret %<>%
          dplyr::mutate(!!paste("Sensitivity_FDR_10pct.ctrl_plasmid",genes,sep=".") := NA,
                        !!paste("Sensitivity_FDR_5pct.ctrl_plasmid",genes,sep=".") := NA,
                        !!paste("Sensitivity_FDR_10pct.treat_plasmid",genes,sep=".") := NA,
                        !!paste("Sensitivity_FDR_5pct.treat_plasmid",genes,sep=".") := NA)
      }
    }
    ret %<>% dplyr::select(-dummy)
    if(type == "AUROC"){
      return(list(AUROC = ret))
    }else{
      return(list(AUPrRc = ret))
    }
  },
  error = function(e) stop(paste(".make_dummy_bagel: unable to build dummy Bagel QC output:",e))
  )
}


.add_masks_mocked_columns <- function(data, b2f, lib){
  tryCatch({
    cols2mask <- dplyr::case_when((is.null(b2f) && is.null(lib)) ~ list(mask=c(fgcQC::mask_mocked_columns$bcl2fastq,
                                                                     fgcQC::mask_mocked_columns$library)),
                                  (is.null(b2f) && !is.null(lib)) ~ list(mask=fgcQC::mask_mocked_columns$bcl2fastq),
                                  (!is.null(b2f) && is.null(lib)) ~ list(mask=fgcQC::mask_mocked_columns$library),
                                  (!is.null(b2f) && !is.null(lib)) ~ list(mask=NULL))
    if(is.null(unlist(cols2mask)))
      return(data)
    data[,unlist(cols2mask)] <- NA
    return(data)
  },
  error = function(e) stop(paste(".add_masks_mocked_columns: unable to mask mocked columns:",e))
  )
}


#' QC_fgc_crispr_data
#'
#' Produce Quality Control (QC) metrics from the following FGC CRISPR screen data inputs: an analysis config JSON file, a combined counts csv file, a Bagel results tsv file (Control vs Plasmid), a gRNA library tsv file ('cleanr.tsv'), and one or more bcl2fastq2 output JSON files ('Stats.json'). QC output will be returned for a single named comparison.
#'
#' @param analysis_config A path to a valid analysis config JSON file to be used as input for the AZ-CRUK CRISPR analysis pipeline.
#' @param combined_counts A path to a valid combined counts csv file (produced by the AZ-CRUK CRISPR analysis pipeline).
#' @param bagel_ctrl_plasmid A path to a valid Bagel tsv results file for Control vs Plasmid. If `NULL`, no such file is available.
#' @param bagel_treat_plasmid A path to a valid Bagel tsv results file for Treatment vs Plasmid. If `NULL`, no such file is available.
#' @param bcl2fastq A character string giving one or more paths to valid `bcl2fastq2` summary output JSON files (paths separated by commas). May be `NULL` if `mock_missing_data` argument is `TRUE`.
#' @param library A valid path to a library tsv file in which the first column gives the sgRNA sequence and the second column gives the sgRNA ID (produced by the AZ-CRUK CRISPR reference data generation pipeline). May be `NULL` if `mock_missing_data` argument is `TRUE`.
#' @param comparison_name A character string naming a single comparison to extract QC data for (should correspond to the comparison name used in the analysis config JSON file).
#' @param output A character string giving an output file name for the csv results. If `NULL`, do not write out any results.
#' @param output_R_object A character string giving an output file name for the returned R object (useful for trouble-shooting). If `NULL`, do not save the R object.
#' @param norm_method A character string naming a normalization method for the count data. Can be `median_ratio` or `relative`. Defaults to `median_ratio`.
#' @param mock_missing_data A logical indicating whether any missing inputs should be mocked or not. Defaults to `FALSE`.
#'
#' @return A list containing the following elements:
#' * `qc_metrics` - A data frame containing QC metrics as columns and samples as rows; this data will also be written to the `output` file, if not `NULL`.
#' * `comparisons` - A data frame of samples belonging to the focal comparison.
#' * `seq_metrics` - A data frame of CI sequencing metrics at both flowcell and sample levels.
#' * `log2FC` - A list containing normalized counts and logFC data frames at both the gRNA and gene level.
#' * `bagel_ROC` - A list containing Bagel Bayes Factor data with `True_Positive_Rate` and `False_Positive_Rate` columns for specific `gene_sets`.
#' * `bagel_PrRc` - A list containing Precision-Recall data for different sample comparisons.
#' @details There are 7 main steps in the function:
#' 1. Check inputs (mock missing data, if needed).
#' 2. Read data.
#' 3. QC for sequencing metrics (merge with qc data in analysis config).
#' 4. Normalize counts and calculate logFC data.
#' 5. QC for counts and logFC data.
#' 6. QC for Bagel Bayes Factors.
#' 7. Wrap up and write out results (mask any mocked columns).
#' @md
#' @author Alex T. Kalinka, \email{alex.kalinka@@cancer.org.uk}
#' @importFrom dplyr mutate select filter left_join inner_join case_when
#' @importFrom tibble add_column
#' @importFrom magrittr %<>%
#' @export
QC_fgc_crispr_data <- function(analysis_config, combined_counts, bagel_ctrl_plasmid, bagel_treat_plasmid,
                               bcl2fastq, library, comparison_name, output, output_R_object,
                               norm_method = "median_ratio", mock_missing_data = FALSE){
  ### 1. Prep.
  .print_progress("Checking data inputs")
  ## Do we need to mock missing inputs?
  if(mock_missing_data){
    mask_bcl2fastq <- bcl2fastq
    mask_library <- library
    mock_file_paths <- fgcQC::mock_missing_FGC_data(analysis_config, combined_counts, bcl2fastq, library)
    analysis_config <- mock_file_paths$analysis_config
    bcl2fastq <- mock_file_paths$bcl2fastq
    library <- mock_file_paths$library
  }

  ## File path checks.
  .check_qc_inputs(analysis_config, combined_counts, bagel_ctrl_plasmid, bcl2fastq, library, comparison_name,
                   output, output_R_object, norm_method, bagel_treat_plasmid)

  ## Analysis config.
  json_list <- fgcQC::read_analysis_config_json(analysis_config)
  comparisons <- fgcQC::extract_analysis_comparisons(json_list)

  comparisons %<>%
    # Limit to samples in given comparison.
    dplyr::filter(comparison == comparison_name)
  .check_comparison(comparisons, comparison_name, analysis_config)

  # Prep qc data from analysis config.
  qc_config <- json_list$qc %>%
    dplyr::select(-indexes, -label) %>%
    dplyr::filter(name %in% comparisons$sample)
  ##-----------------------------------


  ### 2. Read data.
  .print_progress("Reading data")
  # counts.
  tryCatch(counts <- read.delim(combined_counts, sep="\t", header=T, stringsAsFactors = F),
           error = function(e) stop(paste("QC_fgc_crispr_data: unable to read combined counts file",combined_counts,":",e)))
  # library file.
  tryCatch(library <- read.delim(library, skip = 1, header=F, stringsAsFactors = F),
           error = function(e) stop(paste("QC_fgc_crispr_data: unable to read library file",library,":",e)))
  # Bagel results.
  if(!is.null(bagel_ctrl_plasmid)){
    tryCatch(bagel_ctrl_plasmid <- read.delim(bagel_ctrl_plasmid, sep="\t", header=T, stringsAsFactors = F),
             error = function(e) stop(paste("QC_fgc_crispr_data: unable to read bagel ctrl-vs-plasmid file",bagel_ctrl_plasmid,":",e)))
  }else{
    bagel_ctrl_plasmid <- NULL
  }
  if(!is.null(bagel_treat_plasmid)){
    tryCatch(bagel_treat_plasmid <- read.delim(bagel_treat_plasmid, sep="\t", header=T, stringsAsFactors = F),
             error = function(e) stop(paste("QC_fgc_crispr_data: unable to read bagel treatment-vs-plasmid file",
                                            bagel_treat_plasmid,":",e)))
  }else{
    bagel_treat_plasmid <- NULL
  }
  ##-----------------------------------


  ### 3. QC for sequencing metrics (bcl2fastq2 output).
  .print_progress("QC for sequencing metrics")
  b2f <- .prep_bcl2fastq(bcl2fastq, comparisons, json_list)

  if(any(!b2f$SampleName %in% qc_config$name))
    stop(paste("QC_fgc_crispr_data: analysis config QC data missing for:",b2f$SampleName[!b2f$SampleName %in% qc_config$name]))

  ## Merge with QC section in analysis config.
  qc_metrics <- fgcQC::merge_bcl2fastq_config_qc(b2f, qc_config, comparisons)

  # Are any samples missing?
  miss_samps <- setdiff(comparisons$sample, qc_metrics$SampleName)
  if(length(miss_samps) > 0)
    stop(paste("QC_fgc_crispr_data: could not find the following comparison samples in the sequencing output (bcl2fastq):",miss_samps))
  ##-----------------------------------


  ### 4. Normalize counts and calculate logFC data for QC metric calculations.
  .print_progress("Normalizing and calculating logFC data")
  # Output determined by screen type.
  if(comparisons$type[1] != "a"){
    ## Normalize counts.
    if(norm_method == "median_ratio"){
      counts_norm <- fgcQC::normalize_library_depth_median_ratio(counts)
    }else{
      counts_norm <- fgcQC::normalize_library_depth_relative(counts, 2e7)
    }

    ## Log2 fold change at gRNA and gene levels.
    control_samp <- comparisons$sample[comparisons$class == "control"]
    plasmid_samp <- comparisons$sample[comparisons$class == "plasmid"]
    # 1. Control vs Plasmid.
    lfc.ctrl_pl <- fgcQC::calc_log2_fold_change_gRNAs(counts_norm, ref = plasmid_samp, comp = control_samp)
    lfc.ctrl_pl.genes <- fgcQC::calc_log2_fold_change_genes(lfc.ctrl_pl)

    # Add 'treatment' samples if screen goal is not 'lethality'.
    if(comparisons$goal[1] != "lethality"){
      treat_samp <- comparisons$sample[comparisons$class == "treatment"]
      # 2. Treatment vs Plasmid.
      lfc.treat_pl <- fgcQC::calc_log2_fold_change_gRNAs(counts_norm, ref = plasmid_samp, comp = treat_samp)
      lfc.treat_pl.genes <- fgcQC::calc_log2_fold_change_genes(lfc.treat_pl)
    }else{
      treat_samp <- NULL
      lfc.treat_pl <- NULL
    }

    # logfc for replicate logfc correlations.
    if(length(control_samp) > 1){
      lfc.ctrl_pl.repl <- fgcQC::calc_log2_fold_change_gRNAs(counts_norm, ref = plasmid_samp, comp = control_samp[1]) %>%
        dplyr::inner_join(fgcQC::calc_log2_fold_change_gRNAs(counts_norm, ref = plasmid_samp, comp = control_samp[2]),
                          by = "sgRNA", suffix = c(".repl_1",".repl_2"))
    }else{
      lfc.ctrl_pl.repl <- NULL
    }

    if(!is.null(treat_samp)){
      if(length(treat_samp) > 1){
        lfc.treat_pl.repl <- fgcQC::calc_log2_fold_change_gRNAs(counts_norm, ref = plasmid_samp, comp = treat_samp[1]) %>%
          dplyr::inner_join(fgcQC::calc_log2_fold_change_gRNAs(counts_norm, ref = plasmid_samp, comp = treat_samp[2]),
                            by = "sgRNA", suffix = c(".repl_1",".repl_2"))
      }else{
        lfc.treat_pl.repl <- NULL
      }
    }else{
      lfc.treat_pl.repl <- NULL
    }
    ##-----------------------------------


    ### 5. QC metric calculations for counts and logFC data.
    .print_progress("QC for counts and logFC data")
    qc_metrics %<>%
      ### QC for count data.
      fgcQC::assemble_counts_QC(counts, counts_norm, library, plasmid_samp, control_samp, treat_samp) %>%
      ### QC for logFC data.
      fgcQC::assemble_logfc_QC(comparisons$goal[1], lfc.ctrl_pl, lfc.treat_pl, lfc.ctrl_pl.repl, lfc.treat_pl.repl, library)
    ##-----------------------------------


    ### 6. QC for Bagel binary classification data.
    .print_progress("QC for Bagel Bayes Factors")
    if(!is.null(bagel_ctrl_plasmid)){
      bagel_roc <- fgcQC::add_bagel_ROC_gene_sets(bagel_ctrl_plasmid, bagel_treat_plasmid,
                                                fgcQC::crispr_gene_sets$essential)
      bagel_prroc <- fgcQC::add_bagel_PRROC_gene_sets(bagel_ctrl_plasmid, bagel_treat_plasmid,
                                                fgcQC::crispr_gene_sets$essential)
    }else{
      bagel_roc <- .make_dummy_bagel(fgcQC::crispr_gene_sets$essential, "AUROC")
      bagel_prroc <- .make_dummy_bagel(fgcQC::crispr_gene_sets$essential, "AUPrRc")
    }

    qc_metrics %<>%
      ## AUROC.
      tibble::add_column(bagel_roc$AUROC, .before = "SampleId") %>%
      ## AUPrRc.
      tibble::add_column(bagel_prroc$AUPrRc, .before = "SampleId")
    ##-----------------------------------
  }


  ### 7. Wrap up, write out results, and prep return object.
  .print_progress("Wrapping up and saving QC results")
  # Mask mocked columns, if needed.
  if(mock_missing_data){
    qc_metrics <- .add_masks_mocked_columns(qc_metrics, mask_bcl2fastq, mask_library)
  }

  # Write out QC metrics to a csv file.
  if(!is.null(output)){
    tryCatch(
      write.csv(qc_metrics, file = output, row.names=F, quote=F),
      error = function(e) stop(paste("QC_fgc_crispr_data: unable to write QC metric results to:",output,":",e))
    )
  }

  # Build return object.
  ret <- list()
  ret$qc_metrics <- qc_metrics
  ret$comparisons <- comparisons
  ret$seq_metrics <- b2f
  if(comparisons$type[1] != "a"){
    ret$log2FC <- list()
    ret$log2FC$control_vs_plasmid.gRNA <- lfc.ctrl_pl
    ret$log2FC$control_vs_plasmid.gene <- lfc.ctrl_pl.genes
    ret$bagel_ROC <- list()
    ret$bagel_ROC$bagel_ctrl_plasmid <- bagel_roc$bagel_ctrl_plasmid
    ret$bagel_PrRc <- list()
    ret$bagel_PrRc$bagel_ctrl_plasmid <- bagel_prroc$bagel_ctrl_plasmid
    if(comparisons$goal[1] != "lethality"){
      ret$log2FC$treatment_vs_plasmid.gRNA <- lfc.treat_pl
      ret$log2FC$treatment_vs_plasmid.gene <- lfc.treat_pl.genes
      ret$bagel_ROC$bagel_treat_plasmid <- bagel_roc$bagel_treat_plasmid
      ret$bagel_PrRc$bagel_treat_plasmid <- bagel_prroc$bagel_treat_plasmid
    }
  }else{
    ret$log2FC <- NA
    ret$bagel_ROC <- NA
    ret$bagel_PrRc <- NA
  }
  class(ret) <- "fgcQC"

  # Save R object.
  if(!is.null(output_R_object)){
    if(!grepl("\\.rds$",output_R_object)){
      output_R_object <- paste(output_R_object,".rds",sep="")
    }
    tryCatch(saveRDS(ret, file = output_R_object, compress="xz", version=2),
             error = function(e) stop(paste("QC_fgc_crispr_data: unable to save R object")))
  }
  ##-----------------------------------

  return(ret)
}
