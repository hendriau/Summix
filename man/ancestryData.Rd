% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ancestryData-data.R
\docType{data}
\name{ancestryData}
\alias{ancestryData}
\title{Example allele frequency data}
\format{
Chromosome, SNP, base pair, reference and alternate alleles, 
 reference allele frequencies, observed allele frequencies
}
\usage{
data(ancestryData)
}
\description{
reference data is 1000 Genomes and NAM. 
1000 Genomes data was downloaded from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/ on May 31, 2018
The IAM Affymetrix 6.0 data were downloaded from ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130711_native_american_admix_train (data accessed October 2018) and had been previously harmonized with the 1000 Genomes data.
Observed data is from gnomAD. gnomAD v2 data was downloaded from https://gnomad.broadinstitute.org/downloads on Oct. 11, 2018
}
\examples{
data("ancestryData")
}
\keyword{datasets}
