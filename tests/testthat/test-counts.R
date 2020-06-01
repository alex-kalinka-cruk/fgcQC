context("test-counts.R")
library(dplyr)

## Tests related to counts QC functions.

## Mock data.
# Mock counts data frame.
counts_mock <- data.frame(sgRNA = letters[1:10],
                          gene = c(rep("A",5),rep("B",5)),
                          plasmid = c(85,22,0,154,42,87,32,0,29,30),
                          ctrl_1 = c(122,54,89,22,45,99,5,15,0,2),
                          ctrl_2 = c(151,61,15,11,51,8,2,12,1,10),
                          treat_1 = c(80,52,512,10,12,18,55,88,1,0),
                          treat_2 = c(0,51,129,12,88,13,12,8,9,2))
# > colSums(counts_mock[,3:ncol(counts_mock)])
#plasmid  ctrl_1  ctrl_2 treat_1 treat_2
#481     453     322     828     324

# Mock seq metrics (bcl2fastq) data frame.
b2f_mock <- data.frame(SampleName = c("plasmid","ctrl_1","ctrl_2","treat_1","treat_2"),
                       NumberReads = c(598,879,384,981,654))
#> 100*colSums(counts_mock[,3:ncol(counts_mock)])/b2f_mock$N
#plasmid   ctrl_1   ctrl_2  treat_1  treat_2
#80.43478 51.53584 83.85417 84.40367 49.54128

# Mock library file.
make_guides <- function(guide_len, num_guides){
  bases <- c("A","T","C","G")
  guides <- rep(NA,num_guides)
  guides <- sapply(guides, function(x) paste(sample(bases, guide_len, replace=T),collapse=""))
  return(guides)
}
library_mock <- data.frame(sgRNA = make_guides(20,10),
                           gene = c(rep("A",5),rep("B",5)))


### Tests
## Zero and low count guides in plasmid.
test_that("zero and low count plasmid guide percentages are correct",{
  qc_zeros <- fgcQC::calc_low_zero_count_plasmid_gRNAs(counts_mock, "plasmid")
  # expect zero percentage = 20%, low count percentage = 40%.
  expect_equal(qc_zeros$percent_zero_count_plasmid_gRNAs, 20)
  expect_equal(qc_zeros$percent_low_count_plasmid_gRNAs, 40)
})

## Percent matching gRNAs.
test_that("percent matching gRNAs are correct",{
  perc_match <- fgcQC::calc_percent_matching_gRNAs(b2f_mock, counts_mock)
  pm <- perc_match$percent_reads_matching_gRNAs
  names(pm) <- NULL
  expect_equal(pm, c(80.43478,51.53584,83.85417,84.40367,49.54128), tolerance = 1e5)
})

## Gini coefficients.
test_that("Gini coefficients produced and correct",{
  gc <- fgcQC::calc_gini_coefficient_counts(counts_mock)$gini_coefficient_counts
  range_corr <- dplyr::between(gc,0,1)
  expect_true(all(range_corr))
})

## Distance correlation GC content vs counts.
test_that("distance correlation for GC content vs counts is working correctly",{
  dc_counts <- fgcQC::calc_dcorr_GC_content_counts(counts_mock)
})

