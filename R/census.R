#' Census sample data
#'
#' Sample from a synthetic population census dataset
#'
#' @format A data frame with 3000 rows and 8 variables:
#' \describe{
#'   \item{CP_CD}{Hierarchical variable 1 (3 levels: 21, 22, 23)}
#'   \item{CDW_CD}{Hierarchical variable 2 (≈80 levels)}
#'   \item{ZONE_CD}{Hierarchical variable 3 (≈100 levels)}
#'   \item{RPRSNTV_SEXDSTN}{Sex (2 levels)}
#'   \item{AGE_FACTOR}{ Age (6 levels)}
#'   \item{BR_ACT}{ Act (3 levels)}
#'   \item{BR_JOJIK}{ Jojik (2 levels)}
#'   \item{ORG_FORM_CD}{ Organization (5 levels)}
#' }
#' @details 'CP_CD','CDW_CD','ZONE_CD' are hierachical key variables and
#' 'RPRSNTV_SEXDSTN','AGE_FACTOR','BR_ACT','BR_JOJIK','ORG_FORM_CD' are key variables.
#' @source Internal sample derived from the synthetic population census dataset.
#' @keywords datasets
#' @usage data(census)
"census"
