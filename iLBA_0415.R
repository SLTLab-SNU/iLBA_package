
#---------------------------------
library(dplyr)
library(data.table)
library(tictoc)

#데이터 만들기------------------------------------------------------------------
MDDF<- fread("SBR_SDfull.csv")
Age <- MDDF$BASE_YEAR - MDDF$YBTH
MDDF$AGE_FACTOR <- cut(Age,breaks=c(0,30,40,50,60,max(Age,na.rm=T)),right=F,include.lowest=T)


hkey<-c('CP_CD','CDW_CD', 'ZONE_CD')
key<-c('RPRSNTV_SEXDSTN','AGE_FACTOR','BR_ACT','BR_JOJIK','ORG_FORM_CD')
vars<-c(hkey,key)
data<- MDDF %>% dplyr::select(all_of(vars))

#SCA, iLBA alogrithm----------------------------------------------------------------------------------------

SCA <- function(freq, B=3) {
  if(!is.vector(freq)) freq <- as.vector(freq)
  if(!sum(freq < B)) return(freq)
  n <- length(freq)
  changer <- ifelse(runif(n) < freq/B, B, 0)
  return(ifelse(freq < B, changer,freq))
}

LBAvec <- function(x , B=3) {
  nLessEqOrigin <- x[1] # B이하 cell의 갯수, K
  nEqSCA <- x[2] # 기저빈도표에 SCA를 적용해 B가 된 cell의 갯수, k
  SumSCA <- x[3] # 상위빈도표의 cell 참값, f_i
  Masked <- nEqSCA*B
  type1 <- type2 <- 0
  if(nLessEqOrigin > 1) {
    lbdSCA <- nEqSCA # v
    ubdSCA <- nEqSCA + nLessEqOrigin*(B-1) # w
    a <- (SumSCA-1) %/% B
    Masked <- a*B + trunc(B/2) + 1
    lbdSystem <- a*B + 1 # r
    ubdSystem <- (a + 1)*B # s
    if (a < 0) {
      Masked <- 0
    } else if (lbdSystem < lbdSCA) {
      Masked <- Masked+B
      type1 <- 1
    }else if(ubdSystem > ubdSCA){
      Masked <- Masked - B
      type2 <- 1
    }
    if(Masked > 0 & Masked < B) Masked <- B
  }
  return(c(Masked, type1, type2))
}


# 기저빈도표(true+masked)와 metadata를 list의 형태로 저장, 사실 이거 속도는 안중요할듯  -------------------------------------------------------------------

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


# 상위빈도표, output은 csv형태

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


# 상위빈도표의 특정 hkey, key 값에 대한 masked freq
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


#실행

savefulltb(data,hkey,B = 5)

saveaggtb(2, key[c(1,2,5)])

hkey_value = '11010'
key_value = c('1','[0,30)','1')
returnaggfreq(2,key[c(1,2,5)],hkey_value,key_value)


#오류있는 경우 문구 잘뜨는지 확인

savefulltb(data,c(hkey,'A'))
savefulltb(data,hkey, rank= c(1,2,3,4))
savefulltb(data,hkey, rank= c(2,1,3))
saveaggtb('hello',key[c(1,2,5)])
returnaggfreq('hello',key[c(1,2,5)],hkey_value,key_value)

saveaggtb(2,key[c(1,2,5)],input = '')
