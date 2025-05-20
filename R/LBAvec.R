#' LBA Vector Masking Calculation
#'
#' Computes masked cell values based on SCA logic and system bounds.
#'
#' @param x Numeric vector of length 3: c(K, k, f_i)
#'   - K: number of original cells ≤ B
#'   - k: number of cells equal to B after SCA
#'   - f_i: true frequency from upper table
#' @param B Threshold parameter (default = 3).
#' @return Numeric vector c(Masked, type1, type2).
#' @export
LBAvec <- function(x, B = 3) {
  nLessEqOrigin <- x[1] # K
  nEqSCA <- x[2]        # k
  SumSCA <- x[3]        # f_i
  Masked <- nEqSCA * B
  type1 <- type2 <- 0
  if (nLessEqOrigin > 1) {
    a <- (SumSCA - 1) %/% B
    Masked <- a * B + trunc(B / 2) + 1
    lbdSCA <- nEqSCA
    ubdSCA <- nEqSCA + nLessEqOrigin * (B - 1)
    lbdSystem <- a * B + 1
    ubdSystem <- (a + 1) * B
    if (a < 0) {
      Masked <- 0
    } else if (lbdSystem < lbdSCA) {
      Masked <- Masked + B
      type1 <- 1
    } else if (ubdSystem > ubdSCA) {
      Masked <- Masked - B
      type2 <- 1
    }
    if (Masked > 0 && Masked < B) Masked <- B
  }
  return(c(Masked, type1, type2))
}
