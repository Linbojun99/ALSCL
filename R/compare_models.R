utils::globalVariables(c("Model", "Metric", "lower", "upper"))

#' Compare Two Stock Assessment Models (ACL vs ALSCL)
#'
#' Computes a comprehensive set of comparison metrics between two model results,
#' including convergence diagnostics, parameter estimates, goodness-of-fit statistics,
#' time-series correlations, and information criteria.
#'
#' @param model1 A list from \code{run_acl} or \code{run_alscl}. Treated as the reference model.
#' @param model2 A list from \code{run_acl} or \code{run_alscl}. Treated as the alternative model.
#' @param data.CatL A data frame with survey catch-at-length data (used for goodness-of-fit).
#' @param model1_name Character. Label for model1. Default is auto-detected from model_type.
#' @param model2_name Character. Label for model2. Default is auto-detected from model_type.
#'
#' @return A list with:
#' \describe{
#'   \item{summary}{Data frame with convergence, objective, AIC, BIC, number of parameters.}
#'   \item{growth}{Data frame comparing VB growth parameters.}
#'   \item{fit_metrics}{Data frame with MSE, RMSE, R-squared, MAPE, etc. for each model.}
#'   \item{correlation}{Data frame with Pearson correlations of SSB, Rec, B, N, CN between models.}
#'   \item{trend_ratio}{Data frame with mean ratio (model2/model1) for each quantity.}
#' }
#'
#' @export
compare_models <- function(model1, model2, data.CatL,
                           model1_name = NULL, model2_name = NULL) {

  # Auto-detect names
  detect_type <- function(m, fallback) {
    if (!is.null(m$model_type)) return(m$model_type)
    # Infer from report structure
    if (!is.null(m$report$FL)) return("ALSCL")
    if (!is.null(m$report$F))  return("ACL")
    return(fallback)
  }
  if (is.null(model1_name)) model1_name <- detect_type(model1, "Model 1")
  if (is.null(model2_name)) model2_name <- detect_type(model2, "Model 2")

  # ===== 1. Convergence & model summary =====
  get_npar <- function(m) length(m$opt$par)
  get_nobs <- function(m, d) sum(d[, 2:ncol(d)] != 0)

  n_obs <- get_nobs(model1, data.CatL)

  obj1 <- as.numeric(model1$opt$objective)[1]
  obj2 <- as.numeric(model2$opt$objective)[1]
  npar1 <- length(model1$opt$par)
  npar2 <- length(model2$opt$par)

  aic1 <- 2 * npar1 + 2 * obj1
  aic2 <- 2 * npar2 + 2 * obj2
  bic1 <- npar1 * log(n_obs) + 2 * obj1
  bic2 <- npar2 * log(n_obs) + 2 * obj2

  summary_df <- data.frame(
    Metric = c("Model Type", "Convergence", "Boundary Hit",
               "Objective (NLL)", "N Parameters", "N Observations",
               "AIC", "BIC", "Delta AIC", "Delta BIC",
               "Max |Gradient|", "Run Time Ratio"),
    Model1 = as.character(c(
      model1_name,
      as.character(model1$converge)[1],
      as.character(model1$bound_hit)[1],
      round(obj1[1], 2),
      npar1,
      n_obs,
      round(aic1, 2),
      round(bic1, 2),
      "0 (ref)",
      "0 (ref)",
      round(as.numeric(model1$final_outer_mgc)[1], 6),
      "1.0 (ref)"
    )),
    Model2 = as.character(c(
      model2_name,
      as.character(model2$converge)[1],
      as.character(model2$bound_hit)[1],
      round(obj2[1], 2),
      npar2,
      n_obs,
      round(aic2, 2),
      round(bic2, 2),
      round(aic2 - aic1, 2),
      round(bic2 - bic1, 2),
      round(as.numeric(model2$final_outer_mgc)[1], 6),
      ""
    )),
    stringsAsFactors = FALSE
  )
  colnames(summary_df) <- c("Metric", model1_name, model2_name)

  # ===== 2. Growth parameters =====
  safe_round <- function(x, digits = 4) {
    if (is.null(x) || length(x) == 0) return(NA)
    round(as.numeric(x)[1], digits)
  }

  growth_df <- data.frame(
    Parameter = c("Linf", "k (vbk)", "t0", "cv_len", "cv_grow", "sigma_index"),
    Model1 = c(
      safe_round(model1$report$Linf, 3),
      safe_round(model1$report$vbk),
      safe_round(model1$report$t0),
      safe_round(model1$report$cv_len),
      safe_round(model1$report$cv_grow),
      safe_round(model1$report$sigma_index)
    ),
    Model2 = c(
      safe_round(model2$report$Linf, 3),
      safe_round(model2$report$vbk),
      safe_round(model2$report$t0),
      safe_round(model2$report$cv_len),
      safe_round(model2$report$cv_grow),
      safe_round(model2$report$sigma_index)
    ),
    stringsAsFactors = FALSE
  )
  colnames(growth_df) <- c("Parameter", model1_name, model2_name)

  # ===== 3. Goodness-of-fit metrics =====
  calc_fit <- function(model, data.CatL) {
    observed  <- as.matrix(data.CatL[, 2:ncol(data.CatL)])
    predicted <- exp(model$report$Elog_index)
    idx <- observed != 0
    obs <- observed[idx]
    pred <- predicted[idx]
    errors <- obs - pred

    MSE    <- mean(errors^2)
    RMSE   <- sqrt(MSE)
    MAE    <- mean(abs(errors))
    sse    <- sum(errors^2)
    sst    <- sum((obs - mean(obs))^2)
    R2     <- 1 - sse / sst
    MAPE   <- mean(abs(errors / obs)) * 100
    SMAPE  <- mean(abs(errors) / ((abs(obs) + abs(pred)) / 2)) * 100

    # Log-scale residuals
    log_resid <- log(obs) - log(pred)
    SDLR   <- stats::sd(log_resid)

    c(MSE = MSE, RMSE = RMSE, MAE = MAE, R2 = R2,
      MAPE = MAPE, SMAPE = SMAPE, SD_log_resid = SDLR)
  }

  fit1 <- calc_fit(model1, data.CatL)
  fit2 <- calc_fit(model2, data.CatL)

  fit_df <- data.frame(
    Metric = names(fit1),
    Model1 = round(fit1, 4),
    Model2 = round(fit2, 4),
    stringsAsFactors = FALSE
  )
  colnames(fit_df) <- c("Metric", model1_name, model2_name)

  # ===== 4. Time-series correlations =====
  quantities <- c("SSB", "Rec", "B", "N", "CN", "CB")
  cors <- numeric(length(quantities))
  ratios <- numeric(length(quantities))

  for (i in seq_along(quantities)) {
    q <- quantities[i]
    v1 <- model1$report[[q]]
    v2 <- model2$report[[q]]
    if (!is.null(v1) && !is.null(v2) && length(v1) == length(v2)) {
      cors[i]   <- round(stats::cor(v1, v2), 4)
      ratios[i] <- round(mean(v2) / mean(v1), 4)
    } else {
      cors[i]   <- NA
      ratios[i] <- NA
    }
  }

  corr_df <- data.frame(
    Quantity    = quantities,
    Correlation = cors,
    Ratio_Mean  = ratios,
    stringsAsFactors = FALSE
  )

  # ===== Print summary =====
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("  Model Comparison:", model1_name, "vs", model2_name, "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n\n")

  cat("--- Convergence & Information Criteria ---\n")
  print(summary_df, row.names = FALSE)

  cat("\n--- Growth Parameters ---\n")
  print(growth_df, row.names = FALSE)

  cat("\n--- Goodness-of-Fit ---\n")
  print(fit_df, row.names = FALSE)

  cat("\n--- Time-Series Agreement ---\n")
  print(corr_df, row.names = FALSE)

  # AIC interpretation
  delta_aic <- aic2 - aic1
  cat(sprintf("\nDelta AIC = %.2f ", delta_aic))
  if (abs(delta_aic) < 2) {
    cat("(models essentially equivalent)\n")
  } else if (delta_aic < 0) {
    cat(sprintf("(%s preferred)\n", model2_name))
  } else {
    cat(sprintf("(%s preferred)\n", model1_name))
  }

  invisible(list(
    summary     = summary_df,
    growth      = growth_df,
    fit_metrics = fit_df,
    correlation = corr_df,
    model1_name = model1_name,
    model2_name = model2_name
  ))
}
