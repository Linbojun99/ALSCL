# 实测数据导入与验证 / Import and validate your own survey observations.
# source("scripts/real_data_workflow.R", encoding="UTF-8")
# 这些是指南辅助函数，不属于 ACL 的公开 API / Guide helpers, not ACL exports.

.numeric_survey_columns <- function(d, name) {
  if (ncol(d) < 3L || nrow(d) < 2L)
    stop(name, ": need a label column, two time columns and two length bins.")
  d <- as.data.frame(d, stringsAsFactors=FALSE, check.names=FALSE)
  if (names(d)[1L] != "LengthBin") stop(name, ": first header must be LengthBin.")
  d[[1L]] <- as.character(d[[1L]])
  for (i in 2:ncol(d)) {
    old <- trimws(as.character(d[[i]]))
    old[is.na(old) | old %in% c("", "NA")] <- NA_character_
    value <- suppressWarnings(as.numeric(old))
    if (any(is.na(value) & !is.na(old)))
      stop(name, ": nonnumeric cell in column ", names(d)[i])
    d[[i]] <- value
  }
  d
}

read_survey_excel <- function(file) {
  if (!requireNamespace("readxl", quietly=TRUE))
    stop("Install the optional reader with install.packages('readxl').")
  sheets <- c("data.CatL", "data.wgt", "data.mat")
  if (!all(sheets %in% readxl::excel_sheets(file)))
    stop("Workbook requires data.CatL, data.wgt and data.mat sheets.")
  setNames(lapply(sheets, function(s) .numeric_survey_columns(
    readxl::read_excel(file, sheet=s, col_types="text", na=c("", "NA"),
      .name_repair="minimal"), s)), sheets)
}

read_survey_csv <- function(folder) {
  sheets <- c("data.CatL", "data.wgt", "data.mat")
  setNames(lapply(sheets, function(s) .numeric_survey_columns(
    utils::read.csv(file.path(folder,paste0(s,".csv")),check.names=FALSE,
      colClasses="character",na.strings=c("", "NA")), s)), sheets)
}

validate_survey_tables <- function(inputs, growth_step=1, zero_action=c("error","missing")) {
  zero_action <- match.arg(zero_action)
  required <- c("data.CatL","data.wgt","data.mat")
  if (!all(required %in% names(inputs))) stop("Three input tables are required.")
  d <- inputs$data.CatL
  if (anyNA(d[[1L]]) || any(!nzchar(d[[1L]])) || anyDuplicated(d[[1L]]))
    stop("LengthBin labels must be present and unique.")
  for (n in required[-1L]) if (!identical(dim(inputs[[n]]),dim(d)) ||
      !identical(names(inputs[[n]]),names(d)) ||
      !identical(as.character(inputs[[n]][[1L]]),as.character(d[[1L]])))
    stop(n, ": labels, periods and their order must match data.CatL.")
  years <- suppressWarnings(as.numeric(names(d)[-1L]))
  if (length(years)<2L || any(!is.finite(years)) || any(diff(years)<=0) ||
      anyDuplicated(years)) stop("Use increasing numeric time headers.")
  if (length(growth_step)!=1L || !is.finite(growth_step) || growth_step<=0 ||
      any(abs(diff(years)-growth_step)>1e-8))
    stop("Time headers must follow growth_step. Insert missing periods explicitly.")
  obs <- as.matrix(d[-1L]); w <- as.matrix(inputs$data.wgt[-1L]); m <- as.matrix(inputs$data.mat[-1L])
  if (!is.numeric(obs) || any(is.infinite(obs)) || any(obs<0,na.rm=TRUE))
    stop("Survey indices must be numeric, nonnegative, finite, or NA.")
  if (zero_action=="error" && any(obs==0,na.rm=TRUE))
    stop("Zero indices found: establish their meaning before excluding them.")
  if (!any(is.finite(obs) & obs>0)) stop("No positive observations.")
  if (!is.numeric(w) || any(!is.finite(w)) || any(w<0))
    stop("Weight must be complete, finite and nonnegative.")
  if (!is.numeric(m) || any(!is.finite(m)) || any(m<0 | m>1))
    stop("Maturity must be complete and in [0,1].")
  # 此处不推断体长边界和生物参数 / Bins and biological parameters are not inferred.
  invisible(data.frame(LengthBins=nrow(d),TimeSteps=ncol(d)-1L,
    Positive=sum(obs>0,na.rm=TRUE),Missing=sum(is.na(obs)),
    Zeros=sum(obs==0,na.rm=TRUE)))
}

fit_survey_tables <- function(inputs, config, model_type=c("acl","alscl")) {
  model_type <- match.arg(model_type)
  required <- c("rec.age","nage","M","sel_L50","sel_L95","growth_step","zero_action")
  if (!all(required %in% names(config)))
    stop("config must supply: ",paste(required,collapse=", "))
  validate_survey_tables(inputs,config$growth_step,config$zero_action)
  if (any(names(config) %in% names(inputs))) stop("Keep input tables out of config.")
  fun <- if (model_type=="acl") ALSCL::run_acl else ALSCL::run_alscl
  do.call(fun,c(inputs,config))
}
