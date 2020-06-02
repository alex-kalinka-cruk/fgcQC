Summary
=======

`fgcQC` is an `R` package that provides implementations of CRISPR QC
metric calcluations pertaining to sequencing, laboratory, and analysis
metrics - both standard and novel metrics are included. It is designed
to ingest data produced by the [AZ-CRUK CRISPR analysis
pipeline](https://bitbucket.astrazeneca.com/projects/DA/repos/az-cruk-crispr-pipeline/browse).

Dependencies
============

-   `R` v3.6.1 (not tested on any other `R` versions).
-   CRAN `R` packages: `dplyr`, `tibble`, `rlang`, `purrr`, `jsonlite`,
    `ineq`, `energy`, `PRROC`
-   Bioconductor `R` packages: `DSS`

Design
======

-   The main entry point is `QC_fgc_crispr_data()` and this function
    takes an analysis config JSON file, count data, bagel output,
    sequencing metrics (bcl2fastq), and a `cleanr.tsv` gRNA library
    file.
-   QC data is produced for a **single comparison** in the analysis
    config.
-   QC metrics are saved to a `csv` (*non-confidential*) file and an `R`
    object containing the intermediate data is saved to an `rds` file
    (*confidential* - to be used for troubleshooting or more in-depth
    analyses and plotting).
-   Each row of the output `csv` file corresponds to a single sample.
-   All the QC calculation functions are exported and available to be
    used stand-alone within the `R` package.
-   Lethality screens are handled by setting treatment file paths to
    `NULL` (see ‘Usage’ section below).

Assumptions
===========

-   A `qc` section exists in the analysis config and contains a
    `date_transduced` variable for each sample.
-   Sample Labels from the analysis config do not contain any
    **confidential** information.
-   Index barcodes are unique within a single SLX\_ID even if sequencing
    has occurred across more than one flowcell.
-   Processing individual analysis comparisons one at a time is
    feasible.
-   Non-targeting guides are named ‘Control’ in the gene column of the
    library file.
-   CRISPR-i and CRISPR-a screens will only return sequencing QC
    metrics.

Usage
=====

In the below command, multiple `bcl2fastq` `Stats.json` files can be
provided in a single comma-separated string. For ‘lethality’ screens,
set the argument `bagel_treat_plasmid` to `NULL`.

    qc <- QC_fgc_crispr_data(analysis_config = "path/to/analysis_config.json",
                             combined_counts = "path/to/combined_counts.csv",
                             bagel_ctrl_plasmid = "path/to/bagel_Control_vs_Plasmid.bf",
                             bagel_treat_plasmid = "path/to/bagel_Treatment_vs_Plasmid.bf",
                             bcl2fastq = "path/to/first/Stats.json,path/to/second/Stats.json",
                             library = "path/to/gRNA-library-file/cleanr.tsv",
                             comparison_name = "Single_Comparison_name",
                             output = "qc-out.csv",
                             output_R_object = "qc-out")

A summary of the return object can be printed to the console:

    qc
    ####  fgcQC summary  ####
    SLX_ID: SLX-19037
    Screen Type: n
    Screen Goal: lethality
    Number of Flowcells: 1
    AUROC for pan-cancer Sanger (Control vs Plasmid): 0.996
    AUROC for Hart essentials (Control vs Plasmid): 0.887
    AUROC for Moderately negative (Control vs Plasmid): 0.932
    AUROC for Weakly negative (Control vs Plasmid): 0.634

Unit Tests
==========

All unit tests can be run from the root of the `R` package directory:

    devtools::test("path/to/fgcQC")

QC metrics
==========

QC metrics to the left of the `SampleId` column are experiment-wide,
comparison-wide, or calculated for pairs of samples (e.g. control vs
plasmid). QC metrics to the right of the `SampleId` column are
sample-specific metrics. In the below tables, `<comp>` is one of
`ctrl_plasmid`, `treat_plasmid`, or `treat_ctrl`, and `<gene_set>` is
the name of a gene set, e.g. `hart_essential`.

Standard
--------

<table>
<colgroup>
<col style="width: 47%" />
<col style="width: 52%" />
</colgroup>
<thead>
<tr class="header">
<th>QC metric</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><code>AUROC.&lt;comp&gt;.&lt;gene_set&gt;</code></td>
<td>Area Under the Receiver Operating Curve for Bagel Bayes Factors.</td>
</tr>
<tr class="even">
<td><code>AUPrRc.&lt;comp&gt;.&lt;gene_set&gt;</code></td>
<td>Area Under the Precision-Recall curve for Bagel Bayes Factors.</td>
</tr>
<tr class="odd">
<td><code>Sensitivity_FDR_10pct.&lt;comp&gt;.&lt;gene_set&gt;</code></td>
<td>Sensitivity at 10% FDR for Bagel Bayes Factors.</td>
</tr>
<tr class="even">
<td><code>NNMD.&lt;comp&gt;.&lt;gene_set&gt;</code></td>
<td>log2 fold change Null-Normalized Mean Difference relative to plasmid.</td>
</tr>
<tr class="odd">
<td><code>NNMD_robust.&lt;comp&gt;.&lt;gene_set&gt;</code></td>
<td>log2 fold change Null-Normalized Mean Difference relative to plasmid defined using the median and median absolute deviation.</td>
</tr>
<tr class="even">
<td><code>repl_log2FC_pearson_corr.&lt;comp&gt;</code></td>
<td>Pearson correlation of log-fold change values for replicate pairs belonging to the same sample type. Only calculated for the first two replicates when there are &gt;2 total replicates.</td>
</tr>
<tr class="odd">
<td><code>cas_activity</code></td>
<td>Cas-9 activity.</td>
</tr>
<tr class="even">
<td><code>virus_batch</code></td>
<td>Batch ID of virus.</td>
</tr>
<tr class="odd">
<td><code>plasmid_batch</code></td>
<td>Batch ID of plasmid.</td>
</tr>
<tr class="even">
<td><code>minimum_split_cell_number</code></td>
<td>Cell density when cell populations are split.</td>
</tr>
<tr class="odd">
<td><code>cell_population_doublings</code></td>
<td>Doubling rate of cell lines.</td>
</tr>
<tr class="even">
<td><code>library_dna_yield_ng.ul</code></td>
<td>Amplified library DNA yield.</td>
</tr>
<tr class="odd">
<td><code>index_plate</code></td>
<td>The index plate ID for prepping the library.</td>
</tr>
<tr class="even">
<td><code>index_plate_well_id</code></td>
<td>The well ID for the index plate.</td>
</tr>
<tr class="odd">
<td><code>percent_low_count_plasmid_gRNAs</code></td>
<td>Percentage of plasmid gRNAs with less than 30 reads.</td>
</tr>
<tr class="even">
<td><code>percent_zero_plasmid_gRNAs</code></td>
<td>Percentage of plasmid gRNAs with zero-count reads.</td>
</tr>
</tbody>
</table>

Experimental
------------

<table>
<colgroup>
<col style="width: 47%" />
<col style="width: 52%" />
</colgroup>
<thead>
<tr class="header">
<th>QC metric</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><code>GICC.&lt;comp&gt;.&lt;gene_set&gt;</code></td>
<td>Generalized Intraclass Correlation Coefficient of normalized counts - multivariate measure of reproducibility (fits a linear model and returns the proportion of the total variance attributable to between-sample variance). Values range between 0 and 1 - higher values indicate more reproducible data.</td>
</tr>
<tr class="even">
<td><code>mahalanobis_dist_ratio.&lt;comp&gt;</code></td>
<td>Mahalanobis distance ratio of control or treatment samples relative to the plasmid for normalized counts. Lower values indicate the control or treatment samples cluster away from the plasmid. Values <span class="math inline">≥</span> 1 indicate control or treatment samples are close to the plasmid and should be investigated.</td>
</tr>
<tr class="odd">
<td><code>dispersion_adj_gRNA.treat_ctrl</code></td>
<td>An empirical-Bayesian shrinkage estimator for the dispersion of gRNA counts - the mean across gRNAs is returned. Only calculated for Treatment-vs-Control.</td>
</tr>
<tr class="even">
<td><code>distcorr_GC_content_counts</code></td>
<td>A ‘distance correlation’ (allows for non-linear, non-monotonic relationships) for gRNA GC content and normalized counts</td>
</tr>
<tr class="odd">
<td><code>distcorr_GC_content_logfc.&lt;comp&gt;</code></td>
<td>A ‘distance correlation’ (allows for non-linear, non-monotonic relationships) for gRNA GC content and log fold changes</td>
</tr>
<tr class="even">
<td><code>gini_coefficient_counts</code></td>
<td>The Gini coefficient for count data. Values range between 0 and 1 with high values indicating most reads belong to a small number of gRNAs.</td>
</tr>
<tr class="odd">
<td><code>norm_counts_GCC_ratio</code></td>
<td>gRNAs with GCC in the 4 bases upstream of the PAM have been found to be inefficient. This metric captures the ratio of median normalized counts for gRNAs with GCC versus the median of all gRNAs.</td>
</tr>
<tr class="even">
<td><code>norm_counts_TT_ratio</code></td>
<td>gRNAs with TT in the 4 bases upstream of the PAM have been found to be inefficient. This metric captures the ratio of median normalized counts for gRNAs with TT versus the median of all gRNAs.</td>
</tr>
<tr class="odd">
<td><code>log2FC_GCC_diff.&lt;comp&gt;</code></td>
<td>gRNAs with GCC in the 4 bases upstream of the PAM have been found to be inefficient. This metric captures the logFC difference for gRNAs with GCC versus all gRNAs.</td>
</tr>
<tr class="even">
<td><code>log2FC_TT_diff.&lt;comp&gt;</code></td>
<td>gRNAs with TT in the 4 bases upstream of the PAM have been found to be inefficient. This metric captures the logFC difference for gRNAs with TT versus all gRNAs.</td>
</tr>
<tr class="odd">
<td><code>NNMD.&lt;comp&gt;.control_guides</code></td>
<td>NNMD for non-targeting control gRNAs. Expected to be &gt; 0 since they do not incur the cost of dsDNA cuts.</td>
</tr>
<tr class="even">
<td><code>NNMD_robust.&lt;comp&gt;.control_guides</code></td>
<td>NNMD_robust for non-targeting control gRNAs. Expected to be &gt; 0 since they do not incur the cost of dsDNA cuts.</td>
</tr>
</tbody>
</table>

Sequencing
----------

When more than one flowcell has been used, the worst example of each
flowcell-level QC metric is reported (min or max depending on the
metric).

<table>
<colgroup>
<col style="width: 47%" />
<col style="width: 52%" />
</colgroup>
<thead>
<tr class="header">
<th>QC metric</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><code>ReadsPF_percent</code></td>
<td>The percentage of reads Passing Filter (PF) for the experiment as a whole. Should be high.</td>
</tr>
<tr class="even">
<td><code>Non_Demultiplexed_Reads_percent</code></td>
<td>The percentage of non-demultiplexed reads for the experiment as a whole.</td>
</tr>
<tr class="odd">
<td><code>Sample_Representation</code></td>
<td>The percentage of reads assigned to each sample.</td>
</tr>
<tr class="even">
<td><code>percent_reads_matching_gRNAs</code></td>
<td>The percentage of reads that match a gRNA sequence in each sample.</td>
</tr>
<tr class="odd">
<td><code>Q30_bases_percent</code></td>
<td>The percentage of reads with quality scores in excess of 30 for each sample. Should be high.</td>
</tr>
<tr class="even">
<td><code>Average_base_quality</code></td>
<td>Average base quality score for each sample.</td>
</tr>
<tr class="odd">
<td><code>Index_OneBaseMismatch_percent</code></td>
<td>The percentage of 1-base index mismatches per sample. Should be low.</td>
</tr>
</tbody>
</table>

Gene sets
=========

-   `pan_cancer_Sanger` - the core-fitness genes defined by Sanger.
-   `hart_essential` - the Hart essentials defined in 2017.
-   `moderately_negative` - 245 genes defined as having logFC between
    -1.2 and -0.5 for ≥ 90% of the DepMap cell lines.
-   `weakly_negative` - 340 genes defined as having logFC between -0.8
    and -0.2 for ≥ 90% of the DepMap cell lines.
-   `hart_nonessential` - the Hart non-essential set.

Bugs, Issues, or Requests
=========================

Please contact [Alex Kalinka](mailto:alex.kalinka@cancer.org.uk)
