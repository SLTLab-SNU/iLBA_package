#' Information Loss Bounded Aggregation (iLBA)
#'
#' Applies the iLBA algorithm when aggregating cell frequencies.
#' Internal; assumes inputs are already validated
#' (iLBA_inputs: length-3 integer vector, mask_thr: integer scalar >= 3).
#'
#' @param iLBA_inputs Integer vector of length 3: c(n_small, n_SCA_thr, sum_small)
#'   - n_small: number of original cells <= mask_thr
#'   - n_SCA_thr: number of cells equal to mask_thr after SCA
#'   - sum_small: true aggregated frequency from full frequency table
#' @param mask_thr Integer scalar: masking threshold.
#' @return Integer vector of length 3:
#' \describe{
#'   \item{sum_masked}{Masked frequency after applying iLBA.}
#'   \item{type1}{Indicator (0/1): whether a disclosure risk control at the lower bound was applied.}
#'   \item{type2}{Indicator (0/1): whether a disclosure risk control at the upper bound was applied.}
#' }
#' @keywords internal
#' @details
#' For `n_small > 1`, the base masked sum is the center of an interval `[quotient * mask_thr + 1, (quotient + 1) * mask_thr]`,
#' then adjusted by `+/-mask_thr` if disclosure risk occurs (type1/type2 flags).
#' For `n_small <= 1`, the masked sum is `n_SCA_thr * mask_thr`.
#' `type1` indicates a lower disclosure risk control; `type2` indicates an upper disclosure risk control.
#' `mask_thr` corresponds to **B** in the iLBA literature.
#' @references Park, M.J., Kim, H.J., & Kwon, S. (2024). *Journal of the Korean Statistical Society*. Springer Nature.
apply_iLBA <- function(iLBA_inputs, mask_thr) {
  # Divide iLBA_inputs
  iLBA_inputs <- as.integer(iLBA_inputs)
  n_small <- iLBA_inputs[1L]
  n_SCA_thr <- iLBA_inputs[2L]
  sum_small <- iLBA_inputs[3L]

  # Initialize masked sum as n_SCA_thr * mask_thr
  sum_masked <- n_SCA_thr * mask_thr
  type1 <- 0L
  type2 <- 0L

  # Apply the iLBA algorithm
  if (n_small > 1L) {
    # Define quotient = floor((sum_small - 1) / mask_thr)
    quotient <- (sum_small - 1L) %/% mask_thr

    # Compute default masked sum for n_small > 1
    sum_masked <- quotient * mask_thr + (mask_thr %/% 2L) + 1L

    # Compute lower and upper bounds inferred by SCA
    lbd_SCA <- n_SCA_thr
    ubd_SCA <- n_SCA_thr + n_small * (mask_thr - 1L)

    # Compute lower and upper bounds inferred by the system
    lbd_sys <- quotient * mask_thr + 1L
    ubd_sys <- (quotient + 1L) * mask_thr

    # Compute masked sum for the cases where disclosure risk control is needed
    if (quotient < 0L) {
      sum_masked <- 0L
    } else if (lbd_sys < lbd_SCA) {
      sum_masked <- sum_masked + mask_thr
      type1 <- 1L
    } else if (ubd_sys > ubd_SCA) {
      sum_masked <- sum_masked - mask_thr
      type2 <- 1L
    }

    # Set masked sum as mask_thr instead of [mask_thr/2] + 1 for quotient = 0
    if (sum_masked > 0L && sum_masked < mask_thr) sum_masked <- mask_thr
  }

  as.integer(c(sum_masked, type1, type2))
}
