#' Generate a masked frequency for a specified cell
#'
#'
#' Calculate the masked frequency for a single cell defined by specified hkey value and key value.
#' @param hkeyLevel Integer scalar indicating which hierarchy level to select (integer: 1 indicates the top‑level hierarchical variable).
#' @param key Character vector of key variable names used for grouping.
#' @param hkey_value A value matching one of the hierarchical key's levels at the chosen level.
#' @param key_value A vector of values (in the same order as `key`) specifying the cell.
#' @param input.path String path to the RDS file produced by `savefulltb()`. Defaults to "fulltable.rds".
#'
#'
#' @import tictoc
#' @importFrom magrittr %>%
#' @export
returnaggfreq <- function(hkeyLevel, key, hkey_value, key_value, input.path = "fulltable.rds"){

  # 입력 경로 존재 여부 확인
  if (!file.exists(input.path)) {
    stop(paste0("The specified input file does not exist: ", input.path))
  }

  # Input fulltable 불러오기 (.rds 파일)

  fulltable <- readRDS(input.path)
  fulltb <- fulltable$fulltb
  B <- fulltable$B
  hkey_full <- fulltable$hkey
  hkey_values_full <- fulltable$hkey_values
  key_full <- fulltable$key
  key_values_full <- fulltable$key_values


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


  # hkey_value 유효성 확인

  var <- target
  value <- hkey_value
  if (!(value %in% hkey_values_full[[var]])) {
    stop(paste0("Invalid hkey value: '", value, "' not found in hkey variable '", var, "'.\n"))
  }


  # key_value 유효성 확인
  for (i in seq_along(key)) {
    var <- key[i]
    value <- key_value[i]
    if (!(value %in% key_values_full[[var]])) {
      stop(paste0("Invalid key value: '", value, "' not found in key variable '", var, "'.\n"))
    }
  }


  condition <- fulltb[[target]] == hkey_value
  for(k in 1:length(key)){
    condition <- condition & (fulltb[[key[k]]] == key_value[k])
  }
  subtb <- fulltb[condition,]

  if(nrow(subtb) == 0){
    message('No matching cell found. Check if the hkey_value and key_value are valid.\n')
    return(0)
  }

  SumGrOrigin   <- sum(subtb$N[subtb$N > B])
  nLessEqOrigin <- sum(subtb$N <= B)
  nEqSCA        <- sum(subtb$N_SCA == B)
  SumOrigin     <- sum(subtb$N[subtb$N <= B])

  input_vec <- c(nLessEqOrigin, nEqSCA, SumOrigin)
  result <- LBAvec(input_vec, B = B)
  names(result) <- c('Masked','Type1','Type2')

  N.Masked <- SumGrOrigin + result[['Masked']]

  # hkey 출력 문자열 만들기
  hkey_str <- paste0(target, " = ", hkey_value)

  # key 출력 문자열 만들기
  key_str <- paste0(key, " = ", key_value, collapse = ", ")

  #결과표
  x <- c(target,key,"N.Masked")
  y <- c(hkey_value,key_value,N.Masked)
  table <- rbind(y) %>% as.data.frame()
  colnames(table) <- x

  # 메시지 출력
  cat("N.Masked for specified cell\n",
      "Hierarchical key:", hkey_str, "\n",
      "Key:", key_str, "\n",
      "N.Masked:", N.Masked, "\n\n")
  print(table, row.names = FALSE)
  return(invisible(N.Masked))


}
