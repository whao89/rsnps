% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NCBI_snp_query2.R
\name{NCBI_snp_query2}
\alias{NCBI_snp_query2}
\title{Query NCBI's dbSNP for information on a set of SNPs}
\usage{
NCBI_snp_query2(SNPs)
}
\arguments{
\item{SNPs}{A vector of SNPs (rs numbers).}
}
\description{
Query NCBI's dbSNP for information on a set of SNPs
}
\examples{
\dontrun{
SNPs <- c("rs332", "rs420358", "rs1837253", "rs1209415715", "rs111068718")
NCBI_snp_query2(SNPs)
# NCBI_snp_query2("123456") ## invalid: must prefix with 'rs'
NCBI_snp_query2("rs420358")
NCBI_snp_query2("rs332") # warning, merged into new one
NCBI_snp_query2("rs121909001") 
NCBI_snp_query2("rs1837253")
NCBI_snp_query2("rs1209415715") # no data available
NCBI_snp_query2("rs111068718") # chromosomal information may be unmapped
}
}

