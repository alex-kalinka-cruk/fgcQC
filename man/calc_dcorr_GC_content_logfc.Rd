% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_dcorr_GC_content_logfc.R
\name{calc_dcorr_GC_content_logfc}
\alias{calc_dcorr_GC_content_logfc}
\title{calc_dcorr_GC_content_logfc}
\usage{
calc_dcorr_GC_content_logfc(logfc, library, col_suffix)
}
\arguments{
\item{logfc}{A data frame of log2 fold change data for a sample contrast (e.g. Control vs Plasmid; gRNAs as rows) as generated by `fgcQC::calc_log2_fold_change_gRNAs`.}

\item{library}{A data frame containing the library file in which the first column gives the sgRNA sequence and the second column gives the sgRNA ID.}

\item{col_suffix}{A character string providing a name to go at the end of the output column name (e.g. `ctrl_plasmid`).}
}
\value{
A data frame with a `SampleName` column and a `distcorr_GC_content_logfc` column.
}
\description{
Calculates the distance correlation between the log2 fold change and the GC content of sgRNAs in a given library.
}
\references{
Szekely, G.J., Rizzo, M.L., and Bakirov, N.K. (2007), Measuring and Testing Dependence by Correlation of Distances, Annals of Statistics, Vol. 35 No. 6, pp. 2769-2794.
}
\author{
Alex T. Kalinka, \email{alex.kalinka@cancer.org.uk}
}
