#' Save masked aggregated frequency table and information loss distribution
#'
#' Applies the iLBA algorithm to a previously saved finest level frequency table to generate
#' a masked aggregated frequency table and the corresponding information loss distribution.
#' Both results are saved as CSV files.
#'
#' @param hkey_level Integer indicating the hierarchical level to aggregate at
#'   (e.g., `1` corresponds to the coarsest level).
#' @param key Vector of key variable names used to define the aggregated frequency table.
#' @param input_path String path to the RDS produced by `save_full_tb()`
#'   (default `"full_tb.rds"`).
#' @param output_tb_path String path to save the aggregated masked table
#'   (CSV; default `"agg_tb.csv"`).
#' @param output_iL_path String path to save the information loss distribution
#'   (CSV; default `"info_loss.csv"`).
#'
#' @return
#' (Invisibly) returns a list with two elements:
#' \describe{
#'   \item{agg_tb}{Aggregated masked table as a `data.table`.}
#'   \item{info_loss}{Distribution of information loss as a `data.frame`.}
#' }
#' Also saves both objects to CSV files.
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
#' save_agg_tb(
#' hkey_level = 3,
#' key = c("gender","age","htype"),
#' input_path = "full_tb.rds",
#' output_tb_path = "agg_tb.csv",
#' output_iL_path = "info_loss.csv"
#' )
#'
#'
#'
#' @importFrom magrittr %>%
#' @export
save_agg_tb <- function(
    hkey_level,
    key,
    input_path = "full_tb.rds",
    output_tb_path = "agg_tb.csv",
    output_iL_path = "info_loss.csv"
) {
  # Load full table object
  if (!file.exists(input_path)) {
    stop("Input file does not exist: ", shQuote(input_path), call. = FALSE)
  }
  full_tb_list <- readRDS(input_path)
  full_tb <- full_tb_list$full_tb
  mask_thr <- full_tb_list$mask_thr
  hkey_full <- full_tb_list$hkey
  key_full <- full_tb_list$key

  # Validate inputs
  if (length(hkey_level) != 1L || is.na(hkey_level) || hkey_level < 1L) {
    stop("`hkey_level` must be a single integer >= 1.", call. = FALSE)
  }
  if (hkey_level > length(hkey_full)) {
    stop("`hkey_level` (", hkey_level, ") exceeds number of hierarchical keys (",
         length(hkey_full), ").", call. = FALSE)
  }
  invalid_key <- setdiff(key, key_full)
  if (length(invalid_key) > 0L) {
    stop("The following `key` variables do not exist in the full table: ",
         paste(invalid_key, collapse = ", "), call. = FALSE)
  }

  # Compute inputs needed for applying iLBA
  hkey <- hkey_full[seq_len(hkey_level)]
  summary_data <- full_tb %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(hkey, key)))) %>%
    dplyr::summarise(
      sum_big = sum(N[N > mask_thr]),
      n_small = sum(N <= mask_thr),
      n_SCA_thr = sum(N_masked == mask_thr),
      sum_small = sum(N[N <= mask_thr]),
      .groups = "drop"
    )

  # Apply iLBA per group
  iLBA_inputs <- as.data.frame(summary_data[, c("n_small", "n_SCA_thr", "sum_small")])
  iLBA_result <- t(apply(
    iLBA_inputs, 1,
    function(v) apply_iLBA(as.integer(v), mask_thr = as.integer(mask_thr))
  ))
  colnames(iLBA_result) <- c("sum_masked", "type1", "type2")

  agg_tb <- data.frame(summary_data, iLBA_result, check.names = FALSE)
  agg_tb$N_masked <- agg_tb$sum_big + agg_tb$sum_masked
  loss <- agg_tb$N_masked - agg_tb$sum_big - agg_tb$sum_small

  # Summarize information loss distribution
  info_loss <- data.frame(Loss = loss) %>%
    dplyr::group_by(Loss) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(perc = round(n / sum(n) * 100, 2))
  info_loss <- rbind(info_loss, data.frame(Loss = "Total", n = sum(info_loss$n), perc = 100))

  # Retain columns and sort
  agg_tb <- agg_tb[, c(hkey, key, "N_masked", "type1", "type2"), drop = FALSE]
  agg_tb <- data.table::as.data.table(agg_tb)
  data.table::setorderv(agg_tb, c(hkey, key))

  # Preview
  cat("Header of aggregated masked table via iLBA\n")
  print(utils::head(agg_tb), row.names = FALSE)
  cat("\nDistribution of Information Loss\n")
  print(as.data.frame(info_loss), row.names = FALSE)

  # Save results in CSV format
  utils::write.csv(agg_tb, file = output_tb_path, row.names = FALSE)
  message("\nAggregated table saved to ", shQuote(output_tb_path))
  utils::write.csv(info_loss, file = output_iL_path, row.names = FALSE)
  message("Information loss saved to ", shQuote(output_iL_path))

  # Return both results invisibly
  invisible(list(
    agg_tb = agg_tb,
    info_loss = info_loss
  ))
}
