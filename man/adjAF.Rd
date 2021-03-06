% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adjAF.R
\name{adjAF}
\alias{adjAF}
\title{Ancestry adjusted allele frequencies}
\usage{
adjAF(data, reference, observed, pi.target, pi.observed)
}
\arguments{
\item{data}{dataframe of unadjusted allele frequency for observed group, K-1 reference ancestry allele frequencies for N SNPs}

\item{reference}{character vector of the column names for K-1 reference ancestry groups. The name of the last reference ancestry group is not included as that group is not used to estimate the adjusted allele frequencies.}

\item{observed}{character value for the column name of observed data group}

\item{pi.target}{numeric vector of the mixture proportions for K reference ancestry groups in the target sample or subject. The order must match the order of the reference columns with the last entry matching the missing reference group.}

\item{pi.observed}{numeric vector of the mixture proportions for K reference ancestry groups for the observed group. The order must match the order of the reference columns with the last entry matching the missing reference group.}
}
\value{
pi: table of input reference ancestry groups, pi.observed, and pi.target

observed.data: name of the data column for the observed group from which adjusted ancestry allele frequency is estimated

Nsnps: number of SNPs for which adjusted AF is estimated

adjusted.AF: data frame of original data with an appended column of adjusted allele frequencies
}
\description{
Adjusts allele frequencies for heterogeneous populations in genetic data given proportion of reference ancestry groups
}
\examples{
data(ancestryData)
tmp.aa<-adjAF(data   = ancestryData,
    reference   = c("ref_AF_eur_1000G"),
    observed    = "gnomad_AF_afr",
    pi.target   = c(0, 1),
    pi.observed = c(.15, .85))
tmp.aa$adjusted.AF[1:5,]
}
\seealso{
\code{\link{summix}} for estimating the proportion of reference ancestry groups and \url{https://github.com/hendriau/Summix} for further documentation
}
\author{
Gregory Matesi, \email{gregory.matesi@ucdenver.edu}

Audrey Hendricks, \email{audrey.hendricks@ucdenver.edu}
}
\keyword{ancestry}
\keyword{genetics}
\keyword{mixture}
