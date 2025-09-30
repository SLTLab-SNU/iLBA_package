#' Return masked frequency for a specific cell
#'
#' Computes the masked count via the iLBA algorithm for a specified cell
#' defined by values of hierarchical key variables and key variables.
#'
#' @param hkey_level Integer indicating the number of hierarchical key levels to use (1 = top level).
#' @param key Character vector of key variable names used for the aggregated cell frequency.
#' @param hkey_value Vector of values for the hierarchical key variables (length must equal `hkey_level`).
#' @param key_value Vector of values for the key variables (same order as `key`).
#' @param input_path String path to the RDS object produced by `save_full_tb()` (default `"full_tb.rds"`).
#' @return Integer scalar: masked count (`N_masked`) for the specified cell.
#' @importFrom magrittr %>%
#' @export
get_agg_freq <- function(
    hkey_level,
    key,
    hkey_value,
    key_value,
    input_path = "full_tb.rds"
) {
  # Load saved full table
  if (!file.exists(input_path)) {
    stop("Input file does not exist: ", shQuote(input_path), call. = FALSE)
  }
  full_tb_list <- readRDS(input_path)
  full_tb           <- full_tb_list$full_tb
  mask_thr          <- full_tb_list$mask_thr
  hkey_full         <- full_tb_list$hkey
  key_full          <- full_tb_list$key
  hkey_values_full  <- full_tb_list$hkey_values
  key_values_full   <- full_tb_list$key_values

  # Validate hkey_level
  if (length(hkey_level) != 1L || is.na(hkey_level) || hkey_level < 1L || hkey_level > length(hkey_full)) {
    stop("`hkey_level` must be a single integer in [1, ", length(hkey_full), "].", call. = FALSE)
  }

  # Validate hkey_value
  if (length(hkey_value) != hkey_level) {
    stop("`hkey_value` must be length ", hkey_level, " (one per hierarchical level).", call. = FALSE)
  }
  target_hkey <- hkey_full[seq_len(hkey_level)]
  for (i in seq_along(target_hkey)) {
    nm <- target_hkey[i]
    if (!(hkey_value[i] %in% hkey_values_full[[nm]])) {
      stop("Invalid `hkey_value`: ", shQuote(hkey_value[i]), " not found in variable ", shQuote(nm), ".", call. = FALSE)
    }
  }

  # Validate key and key_value
  invalid_key <- setdiff(key, key_full)
  if (length(invalid_key) > 0L) {
    stop("The following `key` variables do not exist in the full table: ",
         paste(invalid_key, collapse = ", "), call. = FALSE)
  }
  if (length(key_value) != length(key)) {
    stop("`key_value` must be length ", length(key), " (same as `key`).", call. = FALSE)
  }
  for (i in seq_along(key)) {
    nm <- key[i]
    if (!(key_value[i] %in% key_values_full[[nm]])) {
      stop("Invalid `key_value`: ", shQuote(key_value[i]), " not found in variable ", shQuote(nm), ".", call. = FALSE)
    }
  }

  # Subset the matching rows in the full table
  cond <- rep(TRUE, nrow(full_tb))
  for (i in seq_along(target_hkey)) {
    cond <- cond & (full_tb[[ target_hkey[i] ]] == hkey_value[i])
  }
  for (i in seq_along(key)) {
    cond <- cond & (full_tb[[ key[i] ]] == key_value[i])
  }
  subtb <- full_tb[cond, ]

  # If no match, return 0 with a message
  if (nrow(subtb) == 0L) {
    message("No matching cell found. Check `hkey_value` and `key_value`.")
    return(0L)
  }

  # Summaries for iLBA (NA-free by design from save_full_tb())
  sum_big   <- sum(subtb$N[subtb$N >  mask_thr])
  n_small   <- sum(subtb$N <= mask_thr)
  n_SCA_thr <- sum(subtb$N_masked == mask_thr)
  sum_small <- sum(subtb$N[subtb$N <= mask_thr])

  # Apply iLBA on the small-cell part
  iLBA_inputs <- c(n_small, n_SCA_thr, sum_small)
  res <- apply_iLBA(as.integer(iLBA_inputs), mask_thr = as.integer(mask_thr))
  names(res) <- c("sum_masked", "type1", "type2")

  # Final masked count for the cell
  N_masked <- as.integer(sum_big + res["sum_masked"])
  N_masked
}
