% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_analysis_config_json.R
\name{read_analysis_config_json}
\alias{read_analysis_config_json}
\title{read_analysis_config_json}
\usage{
read_analysis_config_json(file)
}
\arguments{
\item{file}{A path to a valid analysis config JSON file.}
}
\value{
A list containing all the sections in the JSON file.
}
\description{
Extracts sample and comparison data from an analysis config JSON file to be used as input for the AZ-CRUK CRISPR analysis pipeline.
}
\author{
Alex T. Kalinka, \email{alex.kalinka@cancer.org.uk}
}
