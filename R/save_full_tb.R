#' Save finest level frequency table with true and masked counts
#'
#' Generates a finest level frequency table from a microdata set, applies
#' SCA masking to the cell counts, and saves an RDS file containing the
#' masked finest level frequency table, the masking threshold, and the
#' variables used and attribute information.
#'
#' @param data A `data.frame` or `data.table` containing a raw microdata set.
#' @param hkey Vector of hierarchical key variable names used to define the finest level table.
#' @param key Vector of key variable names used to define the finest level table. If `NULL` (default),
#'   all columns in `data` except those in `hkey` are used.
#' @param mask_thr Integer. Masking threshold for SCA/iLBA (default `5`).
#'   This corresponds to \eqn{K} in the iLBA literature.
#' @param hkey_rank Optional numeric vector specifying the order of hierarchical key variables in `hkey`.
#'   Its length must match `hkey`. If `hkey_rank` is provided, `hkey` is reordered accordingly.
#'   For example, if `hkey` is specified as `c('province', 'town', 'county')`,
#'   the adequate input for 'hkey_rank' is `c(1,3,2)`.
#' @param key_thr Integer. Maximum allowed number of unique values for each
#'   `key` variable (default `100`). Variables exceeding this limit are removed.
#' @param output_path String specifying the path where the resulting RDS object
#'   is saved (default `"full_tb.rds"`).
#'
#' @return
#' Invisibly returns a list that is also saved to `output_path`. The list contains:
#' \describe{
#'   \item{full_tb}{A `data.table` containing original counts (`N`)
#'                 and SCA-masked counts (`N_masked`).}
#'   \item{mask_thr}{The masking threshold used for SCA.}
#'   \item{hkey}{Ordered hierarchical key variable names.}
#'   \item{hkey_values}{A named list of sorted unique values for each `hkey`.}
#'   \item{key}{Key variable names.}
#'   \item{key_values}{A named list of sorted unique values for each `key`.}
#' }
#'
#' @details
#' Rows with missing values are removed only if they occur in the selected
#' `hkey` or `key` columns. The finest level frequency table is computed over
#' `c(hkey, key)`, with counts stored in `N`. SCA masking is then applied
#' using `apply_SCA(N, mask_thr)` to produce `N_masked`.
#'
#' @examples
#' save_full_tb(
#' data = census,
#' hkey = c("LA1","LA2","LA3","OA"),
#' key = c("gender", "age", "edu", "mar", "htype"),
#' mask_thr = 5,
#' output_path = "full_tb.rds"
#' )
#'
#'
#'
#' @import data.table
#' @export
save_full_tb <- function(
    data,
    hkey,
    key = NULL,
    mask_thr = 5L,
    hkey_rank = NULL,
    key_thr = 100L,
    output_path = "full_tb.rds"
) {
  # Validate inputs
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame or data.table.", call. = FALSE)
  }
  if (!is.character(hkey) || length(hkey) < 1L) {
    stop("`hkey` must be a non-empty character vector.", call. = FALSE)
  }
  mask_thr <- as.integer(mask_thr)
  if (length(mask_thr) != 1L || is.na(mask_thr) || mask_thr < 3L) {
    stop("`mask_thr` must be a single integer >= 3.", call. = FALSE)
  }
  key_thr <- as.integer(key_thr)
  if (length(key_thr) != 1L || is.na(key_thr) || key_thr < 1L) {
    stop("`key_thr` must be a single integer >= 1.", call. = FALSE)
  }
  if (!is.character(output_path) || length(output_path) != 1L) {
    stop("`output_path` must be a single string path.", call. = FALSE)
  }

  data <- data.table::as.data.table(data)

  # Check existence of input hkey in data
  invalid_hkeys <- setdiff(hkey, names(data))
  if (length(invalid_hkeys) > 0L) {
    stop("The following `hkey` variables do not exist in `data`: ",
         paste(invalid_hkeys, collapse = ", "), call. = FALSE)
  }

  # Decide key variables
  if (is.null(key)) {
    key <- setdiff(names(data), hkey)
  } else {
    if (!is.character(key)) stop("`key` must be a character vector.", call. = FALSE)
    invalid_keys <- setdiff(key, names(data))
    if (length(invalid_keys) > 0L) {
      stop("The following `key` variables do not exist in `data`: ",
           paste(invalid_keys, collapse = ", "), call. = FALSE)
    }
  }

  # Handle hkey rank
  if (!is.null(hkey_rank)) {
    if (length(hkey_rank) != length(hkey)) {
      stop("Length of `hkey_rank` must equal length of `hkey`.", call. = FALSE)
    }
    hkey_ordered <- hkey[order(hkey_rank)]
  } else {
    hkey_ordered <- hkey
  }

  # Drop over-threshold keys
  if (length(key) > 0L) {
    unique_counts <- stats::setNames(sapply(key, function(v) data.table::uniqueN(data[[v]])), key)
    to_remove <- names(unique_counts)[unique_counts > key_thr]
    if (length(to_remove) > 0L) {
      key <- setdiff(key, to_remove)
      message("Removed key variables with > ", key_thr,
              " unique values: ", paste(to_remove, collapse = ", "))
    }
  }

  # Remove rows with NA in the remaining selected columns (hkey + key)
  cols_used <- c(hkey_ordered, key)
  rows_wo_na <- stats::complete.cases(data[, cols_used, with = FALSE])
  na_count <- sum(!rows_wo_na)
  if (na_count > 0L) {
    data <- data[rows_wo_na]
    message(na_count, " rows with missing values in selected columns have been removed.")
  }

  # Build full frequency table and apply SCA
  full_tb <- data[, .(N = .N), by = cols_used][, N_masked := apply_SCA(N, mask_thr)]

  # Save hkey_values and key_values
  hkey_values <- lapply(data[, ..hkey_ordered], function(x) sort(unique(x)))
  names(hkey_values) <- hkey_ordered

  key_values <- if (length(key) > 0L) {
    out <- lapply(data[, ..key], function(x) sort(unique(x)))
    names(out) <- key
    out
  } else {
    stats::setNames(vector("list", 0L), character())
  }

  # Assemble as a list and save in rds format
  result <- list(
    full_tb = full_tb,
    mask_thr = mask_thr,
    hkey = hkey_ordered,
    hkey_values = hkey_values,
    key = key,
    key_values = key_values
  )

  saveRDS(result, file = output_path)

  cat("Hierarchical key variables: ",
          paste0(seq_along(hkey_ordered), ":", hkey_ordered, collapse = " | "), "\n")
  cat("Key variables: ", if (length(key)) paste(key, collapse = ", ") else "(none)", "\n")
  cat("Masking threshold: ", mask_thr, "\n")
  message("Full table saved to ", shQuote(output_path))

  invisible(result)
}
