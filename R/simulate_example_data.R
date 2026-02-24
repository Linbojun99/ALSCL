#' Simulate Example Data for ALSCL Package
#'
#' Generates realistic catch-at-length, weight-at-length, and maturity-at-length
#' data for a generic short-lived species (e.g., krill-like). Supports unequal
#' length bin widths. Useful for package documentation, testing, and tutorials.
#'
#' @param years Numeric vector of years. Default 2000:2025.
#' @param bin_breaks Numeric vector of bin boundary points. Defines n+1
#'   boundaries for n bins. Default \code{c(0, 20, 25, 30, 35, 40, 45, 50, 55, 60)}
#'   giving 9 bins: 0-20, 20-25, 25-30, ..., 55-60. Works with any bin width
#'   (1mm, 2mm, 3mm, 5mm, etc.) and unequal widths.
#' @param Linf Numeric. VB asymptotic length. Default 55.
#' @param vbk Numeric. VB growth coefficient. Default 0.45.
#' @param t0 Numeric. VB theoretical age at length 0. Default -0.1.
#' @param M Numeric. Natural mortality. Default 0.8.
#' @param nage Integer. Maximum age (plus group). Default 7.
#' @param L50_sel Numeric. Length at 50 percent selectivity. Default 28.
#' @param L95_sel Numeric. Length at 95 percent selectivity. Default 36.
#' @param L50_mat Numeric. Length at 50 percent maturity. Default 33.
#' @param L95_mat Numeric. Length at 95 percent maturity. Default 40.
#' @param wgt_a Numeric. Weight-length parameter a in W = a * L^b. Default 2.2e-6.
#' @param wgt_b Numeric. Weight-length parameter b in W = a * L^b. Default 3.2.
#' @param mean_F Numeric. Mean fishing mortality. Default 0.3.
#' @param cv_catch Numeric. CV of observation noise on catches. Default 0.3.
#' @param rec_sigma Numeric. SD of log-recruitment variation. Default 0.5.
#' @param seed Integer or NULL. Random seed for reproducibility. Default NULL
#'   (different result each run). Set e.g. \code{seed = 42} for reproducible output.
#' @param save_csv Logical. Save CSV files to working directory. Default FALSE.
#' @param output_dir Character. Directory for CSV output. Default ".".
#'
#' @return A list with components:
#'   \describe{
#'     \item{data.CatL}{Catch-at-length data.frame (LengthBin x Years)}
#'     \item{data.wgt}{Weight-at-length data.frame (same structure)}
#'     \item{data.mat}{Maturity-at-length data.frame (same structure)}
#'     \item{true_params}{List of true parameter values used in simulation}
#'   }
#'
#' @export
#'
#' @examples
#' # Quick usage with defaults (9 bins like real krill data)
#' sim <- simulate_example_data(seed = 42)
#' sim$data.CatL[, 1]
#' # [1] "0-20" "20-25" "25-30" "30-35" "35-40" "40-45" "45-50" "50-55" "55-60"
#'
#' # Equal-width 5mm bins
#' sim2 <- simulate_example_data(bin_breaks = seq(10, 60, by = 5))
#'
#' # Custom unequal bins
#' sim3 <- simulate_example_data(bin_breaks = c(0, 15, 25, 35, 45, 60))
#'
#' \dontrun{
#' # Run ACL with simulated data
#' sim <- simulate_example_data(seed = 123)
#' result <- run_acl(
#'   data.CatL = sim$data.CatL,
#'   data.wgt  = sim$data.wgt,
#'   data.mat  = sim$data.mat,
#'   rec.age = 1, nage = 7, M = 0.8,
#'   sel_L50 = 28, sel_L95 = 36
#' )
#'
#' # Save CSVs for sharing
#' sim <- simulate_example_data(save_csv = TRUE, output_dir = "example_data")
#' }
simulate_example_data <- function(years      = 2000:2025,
                                  bin_breaks = c(0, 20, 25, 30, 35, 40, 45, 50, 55, 60),
                                  Linf       = 55,
                                  vbk        = 0.45,
                                  t0         = -0.1,
                                  M          = 0.8,
                                  nage       = 7,
                                  L50_sel    = 28,
                                  L95_sel    = 36,
                                  L50_mat    = 33,
                                  L95_mat    = 40,
                                  wgt_a      = 2.2e-6,
                                  wgt_b      = 3.2,
                                  mean_F     = 0.3,
                                  cv_catch   = 0.3,
                                  rec_sigma  = 0.5,
                                  seed       = NULL,
                                  save_csv   = FALSE,
                                  output_dir = ".") {

  if (!is.null(seed)) set.seed(seed)

  nY  <- length(years)
  bin_breaks <- sort(unique(bin_breaks))
  nL  <- length(bin_breaks) - 1

  if (nL < 2) stop("Need at least 3 bin_breaks to define 2+ bins.")

  ages <- 1:nage

  # --- Bin boundaries and midpoints ---
  bin_lower <- bin_breaks[1:nL]
  bin_upper <- bin_breaks[2:(nL + 1)]
  len_mid   <- (bin_lower + bin_upper) / 2

  # --- Build length bin labels ---
  # Always use raw break values: "0-20", "20-25", "25-30", etc.
  # This is mathematically consistent for any bin width (1mm, 2mm, 5mm, ...)
  # and ensures run_acl/run_alscl parsers compute correct midpoints.
  len_labels <- paste0(bin_lower, "-", bin_upper)

  # For pnorm integration, last bin upper is effectively Inf
  bin_upper_eff <- bin_upper
  bin_upper_eff[nL] <- bin_upper[nL] + 50  # large enough to capture tail

  # --- VB: mean length-at-age ---
  mean_La <- Linf * (1 - exp(-vbk * (ages - t0)))

  # --- Selectivity at length (logistic) ---
  sel_L <- 1 / (1 + exp(-log(19) * (len_mid - L50_sel) / (L95_sel - L50_sel)))

  # --- Maturity at length (logistic) ---
  mat_L <- 1 / (1 + exp(-log(19) * (len_mid - L50_mat) / (L95_mat - L50_mat)))

  # --- Weight at length (power law: W = a * L^b) ---
  wgt_L <- wgt_a * len_mid^wgt_b

  # --- Age-length key (proportion of age a in length bin l) ---
  cv_len <- 0.1
  pla <- matrix(0, nrow = nL, ncol = nage)
  for (a in seq_along(ages)) {
    mu_a <- mean_La[a]
    sd_a <- cv_len * mu_a
    for (l in seq_len(nL)) {
      pla[l, a] <- stats::pnorm(bin_upper_eff[l], mu_a, sd_a) -
        stats::pnorm(bin_lower[l], mu_a, sd_a)
    }
    col_sum <- sum(pla[, a])
    if (col_sum > 0) pla[, a] <- pla[, a] / col_sum
  }

  # --- Time-varying fishing mortality ---
  F_trend <- mean_F * (0.5 + 0.5 * sin(seq(-pi/2, 3*pi/2, length.out = nY)))
  F_trend <- pmax(F_trend, 0.05)

  # --- Simulate population dynamics ---
  N <- matrix(0, nrow = nage, ncol = nY)
  log_R0 <- log(1e5)
  log_rec <- log_R0 + stats::rnorm(nY, 0, rec_sigma)
  N[1, ] <- exp(log_rec)

  for (a in 2:nage) {
    N[a, 1] <- N[1, 1] * exp(-M * (a - 1))
  }
  for (y in 2:nY) {
    for (a in 2:nage) {
      N[a, y] <- N[a-1, y-1] * exp(-M - F_trend[y-1])
    }
    N[nage, y] <- N[nage, y] + N[nage, y-1] * exp(-M - F_trend[y-1])
  }

  # --- Convert to catch-at-length ---
  CatL_mat <- matrix(0, nrow = nL, ncol = nY)
  for (y in seq_len(nY)) {
    Fy <- F_trend[y]
    Z  <- M + Fy
    C_a <- N[, y] * (Fy / Z) * (1 - exp(-Z))
    C_L <- as.numeric(pla %*% C_a) * sel_L
    noise <- exp(stats::rnorm(nL, -0.5 * cv_catch^2, cv_catch))
    C_L <- C_L * noise
    CatL_mat[, y] <- pmax(round(C_L, 2), 0)
  }

  # --- Assemble data.frames ---
  year_char <- as.character(years)

  data.CatL <- data.frame(LengthBin = len_labels, CatL_mat, check.names = FALSE)
  colnames(data.CatL) <- c("LengthBin", year_char)

  wgt_mat <- matrix(rep(round(wgt_L, 6), nY), nrow = nL, ncol = nY)
  data.wgt <- data.frame(LengthBin = len_labels, wgt_mat, check.names = FALSE)
  colnames(data.wgt) <- c("LengthBin", year_char)

  mat_mat <- matrix(rep(round(mat_L, 4), nY), nrow = nL, ncol = nY)
  data.mat <- data.frame(LengthBin = len_labels, mat_mat, check.names = FALSE)
  colnames(data.mat) <- c("LengthBin", year_char)

  # --- Save CSV if requested ---
  if (save_csv) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    utils::write.csv(data.CatL, file.path(output_dir, "sim_CatL.csv"), row.names = FALSE)
    utils::write.csv(data.wgt,  file.path(output_dir, "sim_wgt.csv"),  row.names = FALSE)
    utils::write.csv(data.mat,  file.path(output_dir, "sim_mat.csv"),  row.names = FALSE)
    cat("Saved to:", output_dir, "\n")
    cat("  sim_CatL.csv\n  sim_wgt.csv\n  sim_mat.csv\n")
  }

  # --- Return ---
  true_params <- list(
    Linf = Linf, vbk = vbk, t0 = t0, M = M, nage = nage,
    L50_sel = L50_sel, L95_sel = L95_sel,
    L50_mat = L50_mat, L95_mat = L95_mat,
    wgt_a = wgt_a, wgt_b = wgt_b,
    mean_F = mean_F, F_trend = F_trend,
    log_rec = log_rec, cv_catch = cv_catch,
    rec_sigma = rec_sigma, seed = seed,
    bin_breaks = bin_breaks
  )

  result <- list(
    data.CatL   = data.CatL,
    data.wgt    = data.wgt,
    data.mat    = data.mat,
    true_params = true_params
  )

  cat(sprintf("Simulated data: %d length bins x %d years\n", nL, nY))
  cat(sprintf("  Bins: %s\n", paste(len_labels, collapse = ", ")))
  cat(sprintf("  Years: %d-%d\n", min(years), max(years)))
  cat(sprintf("  True VB: Linf=%.1f, k=%.3f, t0=%.2f\n", Linf, vbk, t0))
  cat(sprintf("  True M=%.2f, mean F=%.2f\n", M, mean_F))

  invisible(result)
}
