% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_dcorr_GC_content_counts.R
\name{calc_dcorr_GC_content_counts}
\alias{calc_dcorr_GC_content_counts}
\title{calc_dcorr_GC_content_counts}
\usage{
calc_dcorr_GC_content_counts(counts, library)
}
\arguments{
\item{counts}{A data frame of normalized counts for each sample in the study (samples as columns, gRNAs as rows).}

\item{library}{A data frame containing the library file in which the first column gives the sgRNA sequence and the second column gives the sgRNA ID.}
}
\value{
A data frame with a `SampleName` column and a `distcorr_GC_content_counts` column.
}
\description{
Calculates the distance correlation between normalized counts and the GC content of sgRNAs in a given library.
}
\references{
Szekely, G.J., Rizzo, M.L., and Bakirov, N.K. (2007), Measuring and Testing Dependence by Correlation of Distances, Annals of Statistics, Vol. 35 No. 6, pp. 2769-2794.
}
\author{
Alex T. Kalinka, \email{alex.kalinka@cancer.org.uk}
}
