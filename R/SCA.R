#' Small Cell Adjustment (SCA)
#'
#' Applies a random adjustment to frequencies below a threshold B.
#'
#' @param freq Numeric vector of cell frequencies.
#' @param B Threshold for adjustment (default = 3).
#' @return A numeric vector of adjusted frequencies.
#' @export
SCA <- function(freq, B = 3) {
  if (!is.vector(freq)) freq <- as.vector(freq)
  if (!sum(freq < B)) return(freq)
  n <- length(freq)
  changer <- ifelse(runif(n) < freq / B, B, 0)
  return(ifelse(freq < B, changer, freq))
}
