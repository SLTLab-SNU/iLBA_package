#' Small Cell Adjustment (SCA)
#'
#' Applies SCA algorithm to cell frequencies below a threshold B.
#' This function is intended for internal use and is not meant to be called directly by users.
#'
#' @param freq Numeric vector of cell frequencies.
#' @param B Threshold for masking (default = 3).
#' @return Numeric vector of masked cell frequencies.
#' @keywords internal
#' @references Hundepool, A., Domingo-Ferrer, J., Franconi, L., Giessing, S., Nordholt, E. S., Spicer, K., & de Wolf, P.-P. (2012). *Statistical Disclosure Control*. Wiley.
SCA <- function(freq, B = 3) {
  if (!is.vector(freq)) freq <- as.vector(freq)
  if (!sum(freq < B)) return(freq)
  n <- length(freq)
  changer <- ifelse(runif(n) < freq / B, B, 0)
  return(ifelse(freq < B, changer, freq))
}
