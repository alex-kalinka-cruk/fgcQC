% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarise_samples_bcl2fastq.R
\name{summarise_samples_bcl2fastq}
\alias{summarise_samples_bcl2fastq}
\title{summarise_samples_bcl2fastq}
\usage{
summarise_samples_bcl2fastq(data)
}
\arguments{
\item{data}{A data frame of sequencing metrics derived ultimately from `fgcQC::extract_b2f_json`.}
}
\value{
A summarised data frame with one row per sample.
}
\description{
Summarises sequencing metrics for samples sequenced across multiple lanes and/or flowcells. Minimum values are taken for numeric metrics and count metrics are summed.
}
\author{
Alex T. Kalinka, \email{alex.kalinka@cancer.org.uk}
}
