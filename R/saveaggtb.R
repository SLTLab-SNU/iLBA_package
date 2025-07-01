#' Generate aggregated masked table and distribution of information loss
#'
#' This function aggregates and masks counts from a previously saved full frequency table,
#' then computes the distribution of information loss.
#'
#' @param hkeyLevel Integer scalar indicating which hierarchy level of `hkey` to aggregate (integer: 1 indicates the top‑level hierarchical variable).
#' @param key Character vector of key variable names to group by (must exist in the saved full table).
#' @param input.path String path to the RDS file produced by `savefulltb()`. Default is "fulltable.rds".
#' @param output.table.path String path to write the aggregated masked table CSV. Default is "aggtable.csv".
#' @param output.infoloss.path String path to write the information loss distribution CSV. Default is "infoloss.csv".
#'
#' @import tictoc
#' @importFrom dplyr group_by summarise mutate
#' @importFrom magrittr %>%

#' @export
#'

saveaggtb <- function(hkeyLevel, key, input.path = "fulltable.rds", output.table.path = "aggtable.csv", output.infoloss.path = "infoloss.csv") {

  tic()

  # 입력 경로 존재 여부 확인
  if (!file.exists(input.path)) {
    stop(paste0("The specified input file does not exist: ", input.path))
  }

  # Input fulltable 불러오기 (.rds 파일)

  fulltable <- readRDS(input.path)
  fulltb <- fulltable$fulltb
  B <- fulltable$B
  hkey_full <- fulltable$hkey
  key_full <- fulltable$key


  # hkeyLevel 제대로 입력했는지 확인
  if(hkeyLevel > length(hkey_full)){
    stop(paste0("Invalid hkeyLevel: provided level (", hkeyLevel,
                ") exceeds the number of available hierarchy levels (", length(hkey_full), ")."))
  }
  target <- hkey_full[hkeyLevel]

  # key가 fulltable에 존재하는 키변수인지 확인

  valid_key = !(key %in% key_full)
  if (any(valid_key)) {
    no_var = key[valid_key]
    stop(paste("The following key variables do not exist in fulltable:", paste(no_var, collapse = ', ')),'\n' )
  }


  # 그룹핑 및 요약 (병렬처리 제거, 일반 dplyr 사용)

  summaryData <- fulltb %>%
    group_by(across(all_of(c(target, key)))) %>%
    summarise(
      SumGrOrigin = sum(N[N > B]),
      nLessEqOrigin = sum(N <= B),
      nEqSCA = sum(N_SCA == B),
      SumOrigin = sum(N[N <= B]),
      .groups = "drop"
    )


  # iLBA 마스킹 계산
  iLBAdata <- as.data.frame(summaryData)
  nCol <- ncol(iLBAdata)
  iLBACore <- iLBAdata[, (nCol - 2):nCol]

  Masked <- t(apply(iLBACore, 1, LBAvec, B = B))
  colnames(Masked) <- c("Masked", "Type1", "Type2")

  FinalResult <- cbind(iLBAdata, Masked)
  FinalResult$N.Masked <- FinalResult$SumGrOrigin + FinalResult$Masked  # 최종 Masked count
  Loss <- FinalResult$N.Masked - FinalResult$SumGrOrigin - FinalResult$SumOrigin
  InfoLoss <- data.frame(Loss) %>% group_by(Loss) %>% summarise(n = n()) %>% mutate(perc = round(n/sum(n)*100,2))
  total_row <- data.frame(Loss = "Total", n = sum(InfoLoss$n), perc = sum(InfoLoss$perc))
  InfoLoss <- rbind(InfoLoss, total_row)


  # 필요한 열만 남기고 정렬
  FinalResult <- FinalResult[, c(target, key, "N.Masked","Type1","Type2")]
  FinalResult <- FinalResult[order(do.call(order, FinalResult[, c(target, key), drop = FALSE])), ]


  # 미리보기
  cat("Header of aggregated masked table via iLBA\n\n")
  print(head(FinalResult),row.names = FALSE)
  cat("\n")
  cat("Distribution of Information Loss\n")
  print(as.data.frame(InfoLoss),row.names = FALSE)


  # 결과 저장
  write.csv(FinalResult, file = output.table.path)
  cat("\nAggregated table saved to", output.table.path, "\n")
  write.csv(InfoLoss, file = output.infoloss.path)
  cat("Information Loss saved to", output.infoloss.path, "\n")
  toc()
}
