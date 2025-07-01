#' Return masked frequency for a specific cell
#'
#' Computes the masked count via iLBA algorithm for a specified cell
#' defined by values of hierarchical key variables and key variables.
#'
#' @param hkey.level Integer indicating the number of hierarchical key levels to use (1 indicates top level).
#' @param key Character vector of key variable names used for the aggregated cell frequency.
#' @param hkey.value A vector of values for the hierarchical key variables (length must equal `hkey.level`).
#' @param key.value A vector of values for the key variables (in the same order as `key`).
#' @param input.path String path to load the RDS object produced by `savefulltb()` (default = "fulltable.rds").
#'
#' @importFrom magrittr %>%
#' @export
returnaggfreq <- function(hkey.level, key, hkey.value, key.value, input.path = "fulltable.rds") {

  # Check if input file exists
  if (!file.exists(input.path)) {
    stop("The specified input file does not exist: ", shQuote(input.path))
  }

  # Load saved full table
  fulltable <- readRDS(input.path)
  fulltb <- fulltable$fulltb
  B <- fulltable$B
  hkey_full <- fulltable$hkey
  hkey_values_full <- fulltable$hkey_values
  key_full <- fulltable$key
  key_values_full <- fulltable$key_values

  # Validate hkey.level
  if (hkey.level < 1 || hkey.level > length(hkey_full)) {
    stop(sprintf("Invalid hkey.level: %d (must be between 1 and %d)", hkey.level, length(hkey_full)))
  }

  # Validate hkey.value length
  if (length(hkey.value) != hkey.level) {
    stop(sprintf("`hkey.value` must be a vector of length %d, one value per hierarchical level", hkey.level))
  }

  # Define target hierarchical keys
  target <- hkey_full[seq_len(hkey.level)]

  # Validate each hierarchical key value
  for (i in seq_along(target)) {
    key_name <- target[i]
    if (!(hkey.value[i] %in% hkey_values_full[[key_name]])) {
      stop(sprintf("Invalid hkey value: '%s' not found in variable '%s'", hkey.value[i], key_name))
    }
  }

  # Validate key variables
  invalid_key <- setdiff(key, key_full)
  if (length(invalid_key) > 0) {
    stop("The following key variables do not exist in the full table: ", paste(invalid_key, collapse = ", "))
  }
  if (length(key.value) != length(key)) {
    stop(sprintf("`key.value` must be length %d (same as `key`)", length(key)))
  }
  for (i in seq_along(key)) {
    if (!(key.value[i] %in% key_values_full[[key[i]]])) {
      stop(sprintf("Invalid key value: '%s' not found in variable '%s'", key.value[i], key[i]))
    }
  }

  # Subset matching cell
  condition <- rep(TRUE, nrow(fulltb))
  # apply hierarchical key filters
  for (i in seq_along(target)) {
    condition <- condition & (fulltb[[ target[i] ]] == hkey.value[i])
  }
  # apply other key filters
  for (i in seq_along(key)) {
    condition <- condition & (fulltb[[ key[i] ]] == key.value[i])
  }
  subtb <- fulltb[condition, ]

  # Return 0 if no match
  if (nrow(subtb) == 0) {
    message("No matching cell found. Check if hkey.value and key.value are valid.")
    return(invisible(0))
  }

  # Compute masked value components
  SumGrOrigin   <- sum(subtb$N[subtb$N > B])
  nLessEqOrigin <- sum(subtb$N <= B)
  nEqSCA        <- sum(subtb$N_SCA == B)
  SumOrigin     <- sum(subtb$N[subtb$N <= B])

  # Run iLBA core
  result_vec <- iLBA(c(nLessEqOrigin, nEqSCA, SumOrigin), B = B)
  names(result_vec) <- c("Masked", "Type1", "Type2")
  N.Masked <- SumGrOrigin + result_vec["Masked"]

  # Construct result table
  out <- data.frame(
    matrix(c(hkey.value, key.value, N.Masked), nrow = 1, byrow = TRUE),
    stringsAsFactors = FALSE
  )
  colnames(out) <- c(target, key, "N.Masked")

  # Output
  cat("N.Masked for specified cell\n")
  cat("Hierarchical keys: ", paste(sprintf("%s = %s", target, hkey.value), collapse = ", "), "\n")
  cat("Key vars: ", paste(sprintf("%s = %s", key, key.value), collapse = ", "), "\n")
  cat("N.Masked: ", N.Masked, "\n\n")
  print(out, row.names = FALSE)

  return(invisible(N.Masked))
}
