% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mock_config_qc_section.R
\name{mock_config_qc_section}
\alias{mock_config_qc_section}
\title{mock_config_qc_section}
\usage{
mock_config_qc_section(config, sample_ids, sample_names)
}
\arguments{
\item{config}{A list of config elements retrieved by reading in the analysis config JSON.}

\item{sample_ids}{A character vector naming the sample IDs (indexes) to be mocked.}

\item{sample_names}{A character vector of sample names extracted from the analysis config comparisons.}
}
\value{
A file path to a temporary mock analysis config JSON file.
}
\description{
Mock an analysis config write to a temporary JSON file, and return the temporary file path.
}
