context("test-logfc.R")
library(dplyr)
set.seed(4)


## Mock data.
# Mock counts data frame.
counts_mock <- data.frame(sgRNA = letters[1:10],
                          gene = c(rep("A",5),rep("B",5)),
                          plasmid = c(85,22,0,154,42,87,32,0,29,30),
                          ctrl_1 = c(122,54,89,22,45,99,5,15,0,2),
                          ctrl_2 = c(151,61,15,11,51,8,2,12,1,10),
                          treat_1 = c(80,52,512,10,12,18,55,88,1,0),
                          treat_2 = c(0,51,129,12,88,13,12,8,9,2))
lfc_mock <- counts_mock %>%
  fgcQC::normalize_library_depth_median_ratio() %>%
  fgcQC::calc_log2_fold_change_gRNAs(ref = "plasmid", comp = c("ctrl_1","ctrl_2"))

lfc_mock.repl <- counts_mock %>%
  fgcQC::calc_log2_fold_change_gRNAs(ref = "plasmid", comp = "ctrl_1") %>%
  dplyr::inner_join(fgcQC::calc_log2_fold_change_gRNAs(counts_mock, ref = "plasmid", comp = "ctrl_2"),
                    by = "sgRNA", suffix = c(".repl_1",".repl_2"))

# Mock library file.
make_guides <- function(guide_len, num_guides){
  bases <- c("A","T","C","G")
  guides <- rep(NA,num_guides)
  guides <- sapply(guides, function(x) paste(sample(bases, guide_len, replace=T),collapse=""))
  return(guides)
}
library_mock <- data.frame(V1 = make_guides(20,10),
                           V2 = letters[1:10], stringsAsFactors = F)

## Count Normalization and logFC calculation tests.
lfc.fgc_0008 <- counts.fgc_0008 %>%
  fgcQC::normalize_library_depth_median_ratio() %>%
  fgcQC::calc_log2_fold_change_gRNAs(ref = "Yusa_v3_KO_Plasmid_FGC_batch_1",
                                     comp = c("FGC0008_01_01_04","FGC0008_01_01_05")) %>%
  dplyr::inner_join(cleanr.lfc.fgc_0008, by = c("sgRNA" = "SEQID"))


## Tests.
# Median-ratio normalization and log2FC calculation.
test_that("correlation of logFC with CRISPRcleanR 'avgFC' output is > 0.99",{
  expect_true(cor(lfc.fgc_0008$log2FC, lfc.fgc_0008$avgFC) > 0.99)
})

# NNMD,
test_that("NNMD calculation correct",{
  nnmd <- fgcQC::NNMD(lfc_mock, "gene", "log2FC", "A", "B")
  expect_equal(nnmd, 1.166517, tolerance = 1e-5)
})

# NNMD_robust,
test_that("NNMD calculation correct",{
  nnmd <- fgcQC::NNMD_robust(lfc_mock, "gene", "log2FC", "A", "B")
  expect_equal(nnmd, 1.812433, tolerance = 1e-5)
})

## Distance correlation GC content vs lfc.
test_that("distance correlation for GC content vs logFC is working correctly",{
  dc_lfc <- fgcQC::calc_dcorr_GC_content_logfc(lfc_mock, library_mock, "ctrl_plasmid")$distcorr_GC_content_logfc.ctrl_plasmid
  range_dcorr <- sapply(dc_lfc, dplyr::between, -1, 1)
  expect_true(all(range_dcorr))
})

## Inefficient gRNA (GCC and TT) count ratios.
test_that("inefficient gRNA logFC diffs are correct",{
  library_mock_ineff_guides <- library_mock
  library_mock_ineff_guides$V1[1] <- paste(paste(rep("A",17),collapse=""),"GCC",sep="")
  library_mock_ineff_guides$V1[2:3] <- paste(paste(rep("A",18),collapse=""),"TT",sep="")
  icr <- fgcQC::calc_inefficient_gRNA_logfc(lfc_mock, library_mock_ineff_guides,
                                            "ctrl_plasmid")
  expect_equal(icr$log2FC_GCC_diff.ctrl_plasmid, 1.057514, tolerance = 1e-5)
  expect_equal(icr$log2FC_TT_diff.ctrl_plasmid, 1.974297, tolerance = 1e-5)
})

## Correlation of replicate logFC values.
test_that("replicate logFC correlations are correct",{
  corr.lfc <- fgcQC::calc_replicate_logfc_corr(lfc_mock.repl, "ctrl_plasmid")
  expect_equal(corr.lfc$repl_log2FC_pearson_corr.ctrl_plasmid,
               0.8380559, tolerance = 1e-5)
})

