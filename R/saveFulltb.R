#' Generate Full and Masked Frequency Tables
#'
#' Computes a full frequency table from the input data, applies SCA masking,
#' and saves the result—including metadata—as an RDS file.
#' @param data A data.frame or data.table containing the raw data.
#' @param hkey A character vector of column names to use as hierarchical keys.
#' @param key A character vector of column names to use as key variables for cross-classification. Defaults to all columns not in hkey.
#' @param B An integer scalar specifying the masking threshold B. Default is 3.
#' @param rank An optional numeric vector specifying the ranking of hkey variables. Must be same length as hkey.
#' @param key.threshold An integer specifying the maximum number of unique values allowed for key variables. Variables with more unique values will be removed. Default is 100.
#' @param output.path A character string giving the file path to save the resulting RDS object. Default is 'fulltable.rds'.
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
#' @importFrom data.table as.data.table
#' @import tictoc
#' @importFrom magrittr %>%
#' @export
savefulltb <- function(data, hkey, key = NULL, B = 3 ,rank = NULL, key.threshold = 100, output.path = 'fulltable.rds') {

  tic()

  dt <- as.data.table(data)

  # NA 처리
  na_rows <-sum(!complete.cases(dt)) # NA가 포함되어 있는 row의 갯수
  if (na_rows > 0) {
    dt <- dt %>% na.omit() %>% as.data.table()
    cat(paste(na_rows, "rows with missing values have been removed.\n"))
  }

  # hkey가 데이터에 존재하는 변수인지 확인
  valid_hkey = !(hkey %in% colnames(dt))
  if (any(valid_hkey)) {
    no_var = hkey[valid_hkey]
    stop(paste("The following hkey variables do not exist in data:", paste(no_var, collapse = ', ')),'\n' )
  }

  # key를 입력받지 않으면, data의 열들 중 hkey를 제외한 것을 key로 둠
  col_names <- names(dt)

  if (is.null(key)) {
    key <- setdiff(col_names, hkey)
  } else {

    # key가 데이터에 존재하는 변수인지 확인
    valid_key <- !(key %in% colnames(dt))
    if (any(valid_key)) {
      no_var <- key[valid_key]
      stop(paste("The following key variables do not exist in data:", paste(no_var, collapse = ', ')), '\n')
    }
  }

  # 아직 안쓰는듯..?
  cate_var <- names(dt)[sapply(dt, function(x) is.factor(x) || is.character(x))] # 범주형 변수들 이름
  num_var <- names(df)[sapply(df, is.numeric)] # numeric 변수들 이름

  # 위계 변수 순서 설정

  if (!is.null(rank)) {
    if (length(rank) != length(hkey)) {
      stop("Length of 'rank' must be equal to length of 'hkey'\n")
    }
    hkey_ordered <- hkey[order(rank)]
  } else {
    hkey_ordered <- hkey
  }


  # key 변수들 중 unique한 값이 threshold(100)개가 넘는, numeric으로 간주할 수 있는 변수는 제거

  unique_counts <- dt[, lapply(.SD, uniqueN), .SDcols = key]
  unique_counts <- unlist(unique_counts)

  to_remove <- names(unique_counts[unique_counts > key.threshold])

  if (length(to_remove) > 0) {
    dt[, (to_remove) := NULL]  # 열 제거 (data.table 방식)
    key <- setdiff(key, to_remove)  # key 목록 갱신
    cat("Removed 'key' variables with >100 unique values: ", paste(to_remove, collapse = ", "), "\n")
  }


  # rank 올바른지 확인: 이건 이렇게 하는게 맞는지 잘 모르겠음
  for (i in 1:(length(hkey_ordered) - 1)) {
    tmp1<-hkey_ordered[i]
    tmp2<-hkey_ordered[i+1]
    n_unique_prev <- nrow(unique(dt[,..tmp1], drop = FALSE))
    n_unique_next <- nrow(unique(dt[,..tmp2], drop = FALSE))

    if (n_unique_prev > n_unique_next) {
      stop("'rank' is wrong: number of unique values decreased at hkey position ", i+1,'\n')
    }
  }


  # 기저 빈도표 계산
  fulltb <- dt[, .N, by = c(hkey_ordered, key)]
  fulltb$N_SCA <- SCA(fulltb$N, B)


  # hkey별 unique한 값들을 오름차순 정렬해서 저장
  hkey_values <- lapply(dt[, ..hkey_ordered], function(col) sort(unique(col)))
  names(hkey_values) <- hkey_ordered

  # key별 unique한 값들을 오름차순 정렬해서 저장
  key_values <- lapply(dt[, ..key], function(col) sort(unique(col)))
  names(key_values) <- key


  # list로 metadata까지 합침
  result <- list(
    fulltb = fulltb,
    B = B,
    hkey = hkey_ordered,
    hkey_values = hkey_values,
    key = key,
    key_values = key_values
  )

  #saveaggtb에서 사용자가 참고하기 위한 hkey level 테이블
  level<- data.frame(hkey = hkey_ordered,level=1:length(hkey_ordered))

  # list를 rds로 저장
  saveRDS(result, file = output.path)
  cat("Full table saved to", output.path, "\n")
  cat("Hierachical key variables:\n")
  print(level)
  #cat("Hierachical key variables:",paste(hkey,collapse = ', '),"\n")
  cat("\nKey variables:",paste(key, collapse = ', '),"\n")
  cat("B:",paste(B, collapse = ', '),"\n")

  toc()
}
