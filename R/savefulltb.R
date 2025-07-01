#' Save full frequency table with true and masked cell frequencies
#'
#' Computes a full frequency table from the microdata, applies SCA masking,
#' and saves the results as an RDS file.
#'
#' @param data Data.frame or data.table containing the raw microdata.
#' @param hkey Character vector of hierarchical key variable names.
#' @param key Character vector of key variable names (default = all column names of data except hkey).
#' @param B Threshold for masking (default = 3).
#' @param rank Numeric vector specifying the rank of hierarchical key variables (default = c(1,...,length(hkey)).
#' @param key.threshold Maximum number of unique values allowed for key variables (default = 100). Variables with more unique values will be removed.
#' @param output.path String path to save the resulting RDS object (default = "fulltable.rds").
#'
#' This function saves the following metadata:
#' \describe{
#'   \item{fulltb}{A `data.table` with raw counts (`N`) and masked counts (`N_SCA`).}
#'   \item{B}{The masking threshold used for SCA.}
#'   \item{hkey}{Ordered hierarchical key variable names.}
#'   \item{hkey_values}{Named list of sorted unique values for each hierarchical key.}
#'   \item{key}{Final set of key variables after threshold filtering.}
#'   \item{key_values}{Named list of sorted unique values for each key.}
#' }
#' @importFrom magrittr %>%
#' @export
#'
savefulltb <- function(data, hkey, key = NULL, B = 3, rank = NULL, key.threshold = 100, output.path = "fulltable.rds") {

  dt <- data.table::as.data.table(data)

  # Remove rows with missing values
  na_rows <- dt[, sum(!stats::complete.cases(.SD))]
  if (na_rows > 0) {
    dt <- dt[stats::complete.cases(dt)]
    cat(paste(na_rows, "rows with missing values have been removed.\n"))
  }

  # Check that all hkey variables exist in data
  invalid_hkeys <- setdiff(hkey, names(dt))
  if (length(invalid_hkeys) > 0) {
    stop("The following hkey variables do not exist in data: ", paste(invalid_hkeys, collapse = ", "))
  }

  # If key is not provided, use all columns except hkey
  if (is.null(key)) {
    key <- setdiff(names(dt), hkey)
  } else {
    # Check that all key variables exist in data
    invalid_keys <- setdiff(key, names(dt))
    if (length(invalid_keys) > 0) {
      stop("The following key variables do not exist in data: ", paste(invalid_keys, collapse = ", "))
    }
  }

  # Sort hkey variables by rank if provided
  if (!is.null(rank) && length(rank) != length(hkey)) {
    stop("Length of 'rank' must be equal to length of 'hkey'")
  }
  hkey_ordered <- if (!is.null(rank)) hkey[order(rank)] else hkey

  # Remove key variables with too many unique values
  unique_counts <- dt[, sapply(.SD, data.table::uniqueN), .SDcols = key]
  to_remove <- names(unique_counts[unique_counts > key.threshold])
  if (length(to_remove) > 0) {
    dt[, (to_remove) := NULL]
    key <- setdiff(key, to_remove)
    cat("Removed key variables with >", key.threshold, "unique values: ",
        paste(to_remove, collapse = ", "), "\n")
  }

  # Create full frequency table with masked counts
  fulltb <- dt[, .(N = .N), by = c(hkey_ordered, key)][, N_SCA := SCA(N, B)]

  # Extract sorted unique values for hkey and key variables
  hkey_values <- lapply(dt[, ..hkey_ordered], function(x) sort(unique(x)))
  names(hkey_values) <- hkey_ordered
  key_values <- lapply(dt[, ..key], function(x) sort(unique(x)))
  names(key_values) <- key

  # Assemble result list
  result <- list(
    fulltb = fulltb,
    B = B,
    hkey = hkey_ordered,
    hkey_values = hkey_values,
    key = key,
    key_values = key_values
  )

  # Save result to RDS
  saveRDS(result, file = output.path)
  cat("Full table saved to", shQuote(output.path), "\n")
  cat("Hierarchical key variables:\n")
  for (i in seq_along(hkey_ordered)) {
    cat(i, ":", hkey_ordered[i], "\n")
  }
  cat("Key variables:", paste(key, collapse = ', '), "\n")
  cat("B:", B, "\n")
}
