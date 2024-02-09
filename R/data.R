#' @title ancestryData
#' @description Sample dataset containing reference and observed allele frequencies to be used for examples within the Summix package.
#' @format A data frame with 1000 rows (representing individual SNPs) and 10 columns:
#' \describe{
#'   \item{POS}{Position of SNP on given chromosome.}
#'   \item{REF}{Reference allele}
#'   \item{ALT}{Alternate allele}
#'   \item{CHROM}{Chromosome}
#'   \item{reference_AF_afr}{Allele frequency column of the African reference ancestry.}
#'   \item{reference_AF_eas}{Allele frequency column of the East Asian reference ancestry.}
#'   \item{reference_AF_eur}{Allele frequency column of the European reference ancestry.}
#'   \item{reference_AF_iam}{Allele frequency column of the Indigenous American reference ancestry.}
#'   \item{reference_AF_sas}{Allele frequency column of the South Asian reference ancestry.}
#'   \item{gnomad_AF_afr}{Allele frequency column of the observed gnomAD v3.1.2 African/African American population.}
#' }
#' @source <https://gnomad.broadinstitute.org/downloads#v3>
"ancestryData"
