context("test-bagel.R")
library(dplyr)


## Mock data with known AUROC of 0.825.
# Taken from: https://blog.revolutionanalytics.com/2016/11/calculating-auc.html
category <- c(1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0)
scores <- rev(seq_along(category))
scores[9:10] <- mean(scores[9:10])
data.auc <- data.frame(score=scores, TP=category) %>%
  dplyr::mutate(TP = ifelse(TP==1,T,F))


## AUROC.
test_that("AUROC gives correct result",{
  auc <- fgcQC::calc_AUC(data.auc, "score")
  expect_equal(auc$AUC, 0.825)
})

## Bagel FGC_0008 AUC.
test_that("AUROC for fgc_0008 is correct",{
  auc <- fgcQC::add_bagel_ROC_gene_sets(bagel.fgc_0008, NULL,
                                        fgcQC::crispr_gene_sets$essential)
  # Independently calculated AUC of 0.89 in original 'qc_essentiality.html' output.
  expect_equal(auc$AUROC$AUROC.ctrl_plasmid.hart_essential, 0.89, tolerance = 1e-2)
})

## Bagel FGC_0008 AUPrRc.
test_that("AUPrRc for fgc_0008 is working",{
  auprrc <- fgcQC::add_bagel_PRROC_gene_sets(bagel.fgc_0008, NULL,
                                             fgcQC::crispr_gene_sets$essential)
  expect_true(auprrc$AUPrRc$AUPrRc.ctrl_plasmid.pan_cancer_Sanger > 0.99)
  expect_true(auprrc$AUPrRc$AUPrRc.ctrl_plasmid.weakly_negative < 0.55)
  expect_true(auprrc$AUPrRc$Sensitivity_FDR_10pct.ctrl_plasmid.hart_essential > 0.75)
  expect_true(auprrc$AUPrRc$Sensitivity_FDR_5pct.ctrl_plasmid.hart_essential < 0.75)
})

