#' Save aggregated masked table and distribution of information loss
#'
#' Aggregates and masks counts from a previously saved full frequency table via iLBA algorithm,
#' computes the distribution of information loss, and saves the results as a CSV file.
#'
#' @param hkey.level Integer indicating the level of hierarchical key at which the aggregation is performed (1 indicates the top level).
#' @param key Character vector of key variable names used in the aggregated table.
#' @param input.path String path to load the RDS object produced by `savefulltb()` (default = "fulltable.rds").
#' @param output.table.path String path to save the aggregated masked table in CSV format (default = "aggtable.csv").
#' @param output.infoloss.path String path to save the information loss distribution in CSV format (default = "infoloss.csv").
#'
#' @return
#' Saves two CSV files —
#' (1) the aggregated masked table and (2) the distribution of information loss —
#' to the specified file paths.
#'
#' @examples
#' \donttest{
#' fulltb <- file.path(tempdir(), 'fulltable.rds')
#'
#' savefulltb( census,
#'             hkey = c("CP_CD","CDW_CD","ZONE_CD"),
#'             key = c("RPRSNTV_SEXDSTN", "AGE_FACTOR", "BR_ACT", "BR_JOJIK", "ORG_FORM_CD"),
#'             B=3,
#'             output.path = fulltb)
#'
#' saveaggtb( hkey.level = 2,
#'            key        = c("RPRSNTV_SEXDSTN","AGE_FACTOR","BR_ACT"),
#'            input.path      = fulltb)
#' }
#'
#' @importFrom magrittr %>%
#'
#'
#' @export
saveaggtb <- function(hkey.level, key, input.path = "fulltable.rds", output.table.path = "aggtable.csv", output.infoloss.path = "infoloss.csv") {

  # Check if the input file exists
  if (!file.exists(input.path)) {
    stop(paste0("The specified input file does not exist: ", shQuote(input.path)))
  }

  # Load saved full table
  fulltable <- readRDS(input.path)
  fulltb <- fulltable$fulltb
  B <- fulltable$B
  hkey_full <- fulltable$hkey
  key_full <- fulltable$key

  # Check if hkey.level is valid
  if (hkey.level > length(hkey_full)) {
    stop(paste0("Invalid hkey.level: provided level (", hkey.level,
                ") exceeds the number of hierarchical key variables in the full table (", length(hkey_full), ")."))
  }


  # Check if key is valid
  invalid_key <- setdiff(key, key_full)
  if (length(invalid_key) > 0) {
    stop("The following key variables do not exist in the full table: ", paste(invalid_key, collapse = ", "))
  }

  # Aggregate and summarize counts
  target <- hkey_full[seq_len(hkey.level)]
  summaryData <- fulltb %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(target, key)))) %>%
    dplyr::summarise(
      SumGrOrigin = sum(N[N > B]),
      nLessEqOrigin = sum(N <= B),
      nEqSCA = sum(N_SCA == B),
      SumOrigin = sum(N[N <= B]),
      .groups = "drop"
    )

  # Apply iLBA masking
  iLBAdata <- as.data.frame(summaryData)
  nCol <- ncol(iLBAdata)
  iLBACore <- iLBAdata[, (nCol - 2):nCol]

  Masked <- t(apply(iLBACore, 1, iLBA, B = B))
  colnames(Masked) <- c("Masked", "Type1", "Type2")

  FinalResult <- cbind(iLBAdata, Masked)
  FinalResult$N.Masked <- FinalResult$SumGrOrigin + FinalResult$Masked
  Loss <- FinalResult$N.Masked - FinalResult$SumGrOrigin - FinalResult$SumOrigin

  # Summarize information loss
  InfoLoss <- data.frame(Loss) %>%
    dplyr::group_by(Loss) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    dplyr::mutate(perc = round(n / sum(n) * 100, 2))
  total_row <- data.frame(Loss = "Total", n = sum(InfoLoss$n), perc = 100)
  InfoLoss <- rbind(InfoLoss, total_row)

  # Retain only necessary columns and sort
  FinalResult <- FinalResult[, c(target, key, "N.Masked", "Type1", "Type2")]
  FinalResult <- data.table::as.data.table(FinalResult)
  data.table::setorderv(FinalResult, c(target, key))

  cat("Header of aggregated masked table via iLBA\n\n")
  print(utils::head(FinalResult), row.names = FALSE)
  cat("\nDistribution of Information Loss\n")
  print(as.data.frame(InfoLoss), row.names = FALSE)

  # Save results
  utils::write.csv(FinalResult, file = output.table.path, row.names = FALSE)
  cat("\nAggregated table saved to \"", output.table.path, "\"\n", sep = "")
  utils::write.csv(InfoLoss, file = output.infoloss.path, row.names = FALSE)
  cat("Information Loss saved to \"", output.infoloss.path, "\"\n", sep = "")
}
