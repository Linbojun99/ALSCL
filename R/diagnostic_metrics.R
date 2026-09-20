#' Compute survey-fit and information-criterion diagnostics
#' @param data.CatL Length labels followed by time columns.
#' @param model_result A fitted ACL or ALSCL result.
#' @return A data frame with Metric and Value. Undefined metrics are NA.
#' @details Only finite positive observations are included. MASE divides MAE by
#'   the pooled absolute one-step differences within each length-bin time series;
#'   pairs containing missing or zero observations are excluded. MAD is the mean
#'   absolute deviation of observations about their mean.
#' @export
diagnostic_metrics <- function(data.CatL, model_result) {
  observed <- as.matrix(data.CatL[, -1L, drop = FALSE]); storage.mode(observed) <- "double"
  predicted <- exp(model_result$report$Elog_index)
  if (!identical(dim(observed), dim(predicted))) stop("Observed and predicted dimensions differ.")
  idx <- is.finite(observed) & observed > 0
  if (!any(idx)) stop("No positive finite observations are available.")
  if (any(!is.finite(predicted[idx]))) stop("Predictions are non-finite at observed cells.")
  obs <- observed[idx]; pred <- predicted[idx]; errors <- obs - pred
  mse <- mean(errors^2); mae <- mean(abs(errors)); sst <- sum((obs-mean(obs))^2)
  scale <- NA_real_
  if (ncol(observed) > 1L) {
    pairs <- idx[, -1L, drop=FALSE] & idx[, -ncol(idx), drop=FALSE]
    differences <- abs(observed[, -1L, drop=FALSE] - observed[, -ncol(observed), drop=FALSE])
    if (any(pairs)) scale <- mean(differences[pairs])
  }
  npar <- length(model_result$opt$par)
  if (!npar) npar <- length(model_result$obj$par)
  nll <- model_result$opt$objective
  values <- c(MSE=mse, MAE=mae, MAD=mean(abs(obs-mean(obs))),
    MASE=if(is.finite(scale) && scale>0) mae/scale else NA_real_, RMSE=sqrt(mse),
    Rsquared=if(sst>0) 1-sum(errors^2)/sst else NA_real_, MAPE=100*mean(abs(errors/obs)),
    SMAPE=100*mean(2*abs(errors)/(abs(obs)+abs(pred))),
    'Explained Variance Score'=if(length(obs)>1L && stats::var(obs)>0) 1-stats::var(errors)/stats::var(obs) else NA_real_,
    'Max Error'=max(abs(errors)), AIC=2*npar+2*nll, BIC=npar*log(length(obs))+2*nll)
  data.frame(Metric=names(values),Value=unname(values),row.names=NULL)
}
.acl_max_gradient <- function(m) {
  x <- m$gradient
  if (is.null(x)) x <- m$max_abs_gradient
  if (is.null(x)) x <- m$final_outer_mgc
  if (is.null(x) || !length(x) || any(!is.finite(x))) return(NA_real_)
  max(abs(x))
}
