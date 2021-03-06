% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/specific_ceRNAInteraction.R
\name{get_specific_ceRNAInteractions}
\alias{get_specific_ceRNAInteractions}
\title{Find All Possible ceRNA Interactions Between Identifiers}
\usage{
get_specific_ceRNAInteractions(
  disease_name = NULL,
  ensg_number = NULL,
  gene_symbol = NULL,
  pValue = 0.05,
  pValueDirection = "<",
  limit = 100,
  offset = NULL
)
}
\arguments{
\item{disease_name}{Name of the specific cancer type/dataset. If default is set, all available datasets with corresponding informations are shown.
Fuzzy search available.}

\item{ensg_number}{A vector of ensg number(s). If ensg number is set, gene symbol and gene type must be NULL. One of the three identifiers must be provided.}

\item{gene_symbol}{A vector of gene symbol(s). If gene symbol is set, ensg number and gene type must be NULL. One of the three identifiers must be provided.}

\item{pValue}{Threshold of the FDR adjusted p-value. Default is 0.05.}

\item{pValueDirection}{Direction of the FDR adjusted p-value threshold (<, >). Must be set if pValue is set. Possible values are: "<", ">".}

\item{limit}{Number of results that should be shown. Default value is 100 and can be up to 1000.
For more results please use batches, the provided offset parameter or download the whole dataset.}

\item{offset}{Starting point from where results should be shown.}
}
\value{
A data_frame containing all interactions between genes of interest.
}
\description{
Get all interactions between the given identifiers (ensg_number or gene_symbol).
}
\examples{
get_specific_ceRNAInteractions(disease_name = "pancancer",
             gene_symbol = c("PTENP1","VCAN","FN1"))
\dontrun{
# Do not use both identifiers at the same time
get_specific_ceRNAInteractions(disease_name = "pancancer",
                               ensg_number = c("ENSG00000115414","ENSG00000038427"),
                               gene_symbol = c("VCAN","FN1"))
}
}
