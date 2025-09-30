#' Small Cell Adjustment (SCA)
#'
#' Applies the SCA algorithm to cell frequencies below a threshold.
#' Internal; assumes inputs are already validated
#' (freq: integer, nonnegative, NA-free; mask_thr: integer scalar >= 3).
#'
#' @param freq Integer vector of cell frequencies.
#' @param mask_thr Integer scalar: masking threshold.
#' @return Integer vector of masked cell frequencies.
#' @keywords internal
#' @details
#' For cells with `freq < mask_thr`: returns mask_thr with probability freq / mask_thr, else 0.
#' For cells with `freq >= mask_thr`: returns the original freq.
#' Results are randomized; use `set.seed()` upstream for reproducibility.
#' The parameter `mask_thr` corresponds to **B** in the iLBA literature.
#' @references Hundepool et al. (2012) *Statistical Disclosure Control*. Wiley.
apply_SCA <- function(freq, mask_thr) {
  freq <- as.integer(freq)

  # Bernoulli draw for small cells
  u <- stats::runif(length(freq))
  changer <- ifelse(u < (freq / mask_thr), mask_thr, 0L)

  # Apply SCA rule
  as.integer(ifelse(freq < mask_thr, changer, freq))
}
