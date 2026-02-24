#' @title Run Age- and Length-Structured Statistical Catch-at-Length Model (ALSCL)
#'
#' @description Implements the ALSCL model from Zhang and Cadigan (2022) to assess
#' fish population dynamics from survey catch-at-length data. Unlike the ACL model
#' which tracks F at the age level, ALSCL simultaneously tracks 3D population
#' dynamics across time, age, and length (NLA array), with fishing mortality
#' operating at the length level and a growth transition matrix (G) applied
#' at each time step.
#'
#' @param data.CatL A data frame with length bins in column 1 and survey catch-at-length
#'   in subsequent columns (one per year).
#' @param data.wgt A data frame with length bins in column 1 and weight-at-length
#'   in subsequent columns.
#' @param data.mat A data frame with length bins in column 1 and maturity-at-length
#'   (proportion 0-1) in subsequent columns.
#' @param rec.age Numeric, the recruitment age (e.g. 1).
#' @param nage Numeric, the number of age classes.
#' @param M Numeric, the natural mortality rate (e.g. 0.2).
#' @param sel_L50 Numeric, the length at 50 percent survey selectivity.
#' @param sel_L95 Numeric, the length at 95 percent survey selectivity.
#' @param growth_step Numeric, the growth time step. Default is 1 (annual).
#'   Use smaller values for sub-annual data (e.g. 0.25 for quarterly).
#' @param parameters A list of custom initial parameter values (default NULL = use defaults).
#' @param parameters.L A list of custom lower bounds (default NULL = use defaults).
#' @param parameters.U A list of custom upper bounds (default NULL = use defaults).
#' @param map A list for TMB map to fix parameters. Default NULL uses:
#'   \code{log_sigma_log_F = factor(NA)} and \code{log_t0 = factor(NA)}.
#' @param len_mid Numeric vector, user-specified length bin midpoints (default NULL = auto).
#' @param len_border Numeric vector, user-specified length bin borders (default NULL = auto).
#' @param len_lower Numeric vector, lower bounds of each length bin (default NULL = auto).
#'   First element can be -Inf. Must be same length as number of bins.
#' @param len_upper Numeric vector, upper bounds of each length bin (default NULL = auto).
#'   Last element MUST be finite (not Inf) for the growth transition matrix.
#' @param output Logical, if TRUE auto-save all plots and tables. Default FALSE.
#' @param train_times Numeric, number of nlminb restarts for refinement. Default 1.
#' @param ncores Integer, number of CPU cores. Default 1 (single start).
#' @param silent Logical. If TRUE, suppress all progress messages. Default FALSE.
#'   When ncores is greater than 1, runs parallel multi-start optimization with jittered
#'   initial values for robustness against local optima.
#'
#' @return A list containing:
#' \describe{
#'   \item{report}{TMB report with all estimated quantities (NLA, NL, NA, FL, FA, G, pla,
#'     SSB, B, Rec, CN, CNA, etc.)}
#'   \item{est_std}{Parameter estimates with standard errors from sdreport}
#'   \item{converge}{Convergence message from nlminb}
#'   \item{bound_hit}{Logical, whether any parameter hit its boundary}
#'   \item{year}{Year labels}
#'   \item{len_mid}{Length bin midpoints}
#'   \item{len_label}{Length bin labels from input data}
#'   \item{model_type}{\code{"ALSCL"} identifier for distinguishing from ACL results}
#' }
#'
#' @details
#' \strong{Key differences from \code{run_acl()}:}
#' \itemize{
#'   \item F is modelled at \strong{length} level (\code{dev_log_F} is L x Y),
#'     not at age level (A x Y as in ACL).
#'   \item A \strong{growth transition matrix G} is applied at each time step,
#'     allowing length-dependent processes within cohorts.
#'   \item \strong{F-at-age (FA)} is an emergent property derived from NLA to NA,
#'     not a direct parameter.
#'   \item Additional parameter: \code{log_cv_grow} (CV of growth increment).
#'   \item \code{t0} is log-transformed (\code{log_t0}) rather than raw.
#' }
#'
#' @references Zhang, F., and Cadigan, N. G. (2022). An age- and length-structured
#' statistical catch-at-length model for hard-to-age fisheries stocks.
#' \emph{Fish and Fisheries}, 23(5), 1121-1135.
#'
#' @export
run_alscl <- function(data.CatL, data.wgt, data.mat, rec.age, nage, M,
                      sel_L50, sel_L95, growth_step = 1,
                      parameters = NULL, parameters.L = NULL, parameters.U = NULL,
                      map = NULL, len_mid = NULL, len_border = NULL,
                      len_lower = NULL, len_upper = NULL,
                      output = FALSE, train_times = 1, ncores = 1, silent = FALSE)
{

  # Suppress all cat/message output when silent = TRUE
  if (silent) {
    sink(tempfile())
    on.exit(sink(), add = TRUE)
  }

  if (!requireNamespace("TMB", quietly = TRUE)) {
    stop("Package 'TMB' is required. Install with: install.packages('TMB')")
  }

  total_start <- proc.time()

  nyear <- ncol(data.CatL) - 1
  nlen  <- nrow(data.CatL)
  ages  <- c(rec.age:(rec.age + nage - 1))

  # ===== NA matrix: 1 = observed, 0 = zero/missing =====
  na_matrix <- matrix(1, nrow = nlen, ncol = nyear)
  na_matrix[which(data.CatL[, 2:ncol(data.CatL)] == 0)] <- 0

  # ===== Weight and maturity (match CatL year columns) =====
  weight <- data.wgt[, 2:ncol(data.CatL)]
  mat    <- data.mat[, 2:ncol(data.CatL)]

  # ===== Length bin processing =====
  # Helper: check if length labels are single numbers (e.g. "12") vs ranges (e.g. "5-7")
  contains_only_one_number <- function(s) { return(!grepl("-|<|>", s)) }
  contains_special_char <- function(s) { return(grepl("<|>", s)) }

  if (contains_only_one_number(data.CatL[, 1][1])) {
    # --- Single-value length bins (e.g. 12, 13, 14, ...) ---
    len_values <- as.numeric(data.CatL[, 1])

    if (is.null(len_mid)) {
      len_mid <- len_values
    }

    # Compute bin width and boundaries
    bin_width <- if (nlen > 1) len_values[2] - len_values[1] else 1

    if (is.null(len_lower)) {
      len_lower <- len_values - bin_width/2
    }
    if (is.null(len_upper)) {
      len_upper <- len_values + bin_width/2
    }

    if (is.null(len_border)) {
      len_border <- len_values[-nlen] + bin_width/2
    }

    log_q <- log(mat_func(sel_L50, sel_L95, len_mid))

  } else {
    # --- Range-based length bins (e.g. "5-7", "7-9", "<5", ">49") ---

    if (is.null(len_mid)) {
      range_to_bin <- function(range_str) {
        parts <- strsplit(range_str, "-")[[1]]
        return((as.numeric(parts[2]) - as.numeric(parts[1])) / 2)
      }
      bin <- range_to_bin(data.CatL[, 1][2])

      range_to_median <- function(range_str, isFirst = FALSE, isLast = FALSE,
                                  prev_range_str = "", next_range_str = "") {
        if (contains_special_char(range_str)) {
          parts <- strsplit(range_str, "-")[[1]]
          single_number <- as.numeric(gsub("[^0-9.]", "", parts[1]))
          if (isFirst) {
            next_parts <- strsplit(next_range_str, "-")[[1]]
            bin_w <- (as.numeric(next_parts[2]) - as.numeric(next_parts[1])) / 2
            return(round(single_number - bin_w))
          } else {
            prev_parts <- strsplit(prev_range_str, "-")[[1]]
            bin_w <- (as.numeric(prev_parts[2]) - as.numeric(prev_parts[1])) / 2
            return(round(single_number + bin_w))
          }
        } else {
          parts <- strsplit(range_str, "-")[[1]]
          return(round(as.numeric(parts[2]) - bin))
        }
      }

      len_mid <- sapply(seq_along(data.CatL[, 1]), function(i) {
        prev_str <- if (i > 1) data.CatL[, 1][i - 1] else ""
        next_str <- if (i < nlen) data.CatL[, 1][i + 1] else ""
        range_to_median(data.CatL[, 1][i], isFirst = i == 1, isLast = i == nlen,
                        prev_range_str = prev_str, next_range_str = next_str)
      })
    }

    if (is.null(len_border)) {
      extract_last_number <- function(range_str) {
        if (grepl("<", range_str)) {
          return(as.numeric(gsub("[^0-9.]", "", range_str)))
        } else if (grepl(">", range_str)) {
          return(as.numeric(gsub("[^0-9.]", "", range_str)))
        } else {
          parts <- strsplit(range_str, "-")[[1]]
          return(as.numeric(parts[2]))
        }
      }
      len_border <- sapply(data.CatL[-nlen, 1], extract_last_number)
    }

    # Compute len_lower and len_upper from len_border
    # IMPORTANT: len_upper must be FINITE for the growth transition matrix G.
    # Extract the actual upper bound of the last length bin from data labels.
    if (is.null(len_lower) || is.null(len_upper)) {
      last_label <- as.character(data.CatL[nlen, 1])
      last_parts <- strsplit(last_label, "-")[[1]]
      last_upper <- as.numeric(gsub("[^0-9.]", "", last_parts[length(last_parts)]))

      first_label <- as.character(data.CatL[1, 1])
      first_parts <- strsplit(first_label, "-")[[1]]
      first_lower <- as.numeric(gsub("[^0-9.]", "", first_parts[1]))

      if (is.null(len_lower)) len_lower <- c(first_lower, len_border)
      if (is.null(len_upper)) len_upper <- c(len_border, last_upper)
    }

    log_q <- log(mat_func(sel_L50, sel_L95, len_mid))
  }

  # ===== Compile ALSCL.cpp =====
  alscl_cpp_path <- system.file("extdata", "ALSCL.cpp", package = "ACL")

  if (alscl_cpp_path == "") {
    stop("ALSCL.cpp not found in the package. Please ensure it is in inst/extdata/")
  }

  TMB::compile(file = alscl_cpp_path, "&> /tmp/logfile.log")

  # ===== Build TMB data list =====
  logN_at_len <- as.matrix(log(data.CatL[, 2:ncol(data.CatL)] + 1e-5))

  tmb.data <- list(
    logN_at_len = logN_at_len,
    na_matrix   = na_matrix,
    log_q       = log_q,
    len_mid     = len_mid,
    len_lower   = len_lower,
    len_upper   = len_upper,
    len_border  = len_border,
    age         = ages,
    Y           = nyear,
    A           = nage,
    L           = nlen,
    weight      = as.matrix(weight),
    mat         = as.matrix(mat),
    M           = M,
    growth_step = growth_step
  )

  # ===== Parameters with defaults =====
  custom_bounds <- create_parameters("alscl", parameters, parameters.L, parameters.U)
  params   <- custom_bounds$parameters
  params.L <- custom_bounds$parameters.L
  params.U <- custom_bounds$parameters.U

  # ===== Data-adaptive Linf bounds =====
  # If user did NOT provide custom Linf bounds, auto-adjust based on observed
  # length data. This prevents the Linf-k compensation ridge problem where
  # the optimizer drifts to unrealistic Linf values far from the data range.
  max_len <- max(len_upper)  # largest observed length bin upper bound

  # Initial value: ~10% above max observed length (reasonable VB asymptote)
  if (is.null(parameters) || is.null(parameters$log_Linf)) {
    params$log_Linf <- log(max_len * 1.1)
  }
  # Lower bound: Linf must be at least 70% of max observed length
  if (is.null(parameters.L) || is.null(parameters.L$log_Linf)) {
    params.L$log_Linf <- log(max_len * 0.7)
  }
  # Upper bound: Linf up to ~2x max observed length (generous but prevents runaway)
  if (is.null(parameters.U) || is.null(parameters.U$log_Linf)) {
    params.U$log_Linf <- log(max_len * 2.0)
  }

  # --- Data-adaptive starting values for mean_log_R and log_init_Z ---
  # mean_log_R: estimate from total survey catch per year
  if (is.null(parameters) || is.null(parameters$mean_log_R)) {
    catl_numeric <- as.matrix(data.CatL[, 2:ncol(data.CatL)])
    catl_numeric <- apply(catl_numeric, 2, as.numeric)
    total_per_year <- colSums(catl_numeric, na.rm = TRUE)
    total_per_year <- total_per_year[total_per_year > 0]
    if (length(total_per_year) > 0) {
      params$mean_log_R <- log(mean(total_per_year))
    }
  }
  # log_init_Z: use M + 0.3 as a reasonable approximation for total Z
  if (is.null(parameters) || is.null(parameters$log_init_Z)) {
    params$log_init_Z <- log(M + 0.3)
  }

  # Random effects: dev_log_F is (L x Y) for ALSCL, not (A x Y) as in ACL
  params$dev_log_R  <- rep(0, nyear)
  params$dev_log_F  <- array(0, c(nlen, nyear))  # KEY DIFFERENCE: nlen not nage
  params$dev_log_N0 <- rep(0, nage - 1)

  lower <- unlist(params.L)
  upper <- unlist(params.U)

  # ===== Map =====
  default_map <- list(
    log_sigma_log_F = factor(NA),
    log_t0          = factor(NA)
  )
  if (!is.null(map)) {
    default_map <- modifyList(default_map, map)
  }
  # Clean NULL entries (user can "unfix" by setting map = list(log_t0 = NULL))
  default_map <- default_map[!sapply(default_map, is.null)]
  # Convert to factor(NA) format for TMB
  for (nm in names(default_map)) {
    if (!is.factor(default_map[[nm]])) {
      default_map[[nm]] <- factor(NA)
    }
  }

  rnames <- c("dev_log_R", "dev_log_F", "dev_log_N0")

  # ===== Load DLL =====
  alscl_dll_path <- file.path(dirname(alscl_cpp_path), "ALSCL.dll")
  if (!file.exists(alscl_dll_path)) {
    # Try .so for Linux/Mac
    alscl_dll_path <- file.path(dirname(alscl_cpp_path), paste0("ALSCL", .Platform$dynlib.ext))
  }
  dyn.load(alscl_dll_path)

  # ===== Optimization =====
  run_single_opt <- function(par_init) {
    obj <- TMB::MakeADFun(tmb.data, par_init, random = rnames, map = default_map,
                          DLL = "ALSCL", inner.control = list(trace = FALSE, maxit = 500))

    opt <- nlminb(obj$par, obj$fn, obj$gr, lower = lower, upper = upper,
                  control = list(trace = 0, iter.max = 2000, eval.max = 10000))

    for (tt in seq_len(train_times - 1)) {
      opt <- nlminb(opt$par, obj$fn, obj$gr, lower = lower, upper = upper,
                    control = list(trace = 0, iter.max = 2000, eval.max = 10000))
    }

    return(list(obj = obj, opt = opt))
  }

  if (ncores > 1) {
    # ===== Multi-start parallel optimization =====
    cat(sprintf("ALSCL multi-start: %d parallel starts on %d cores\n", ncores, ncores))

    jitter_params <- function(base, sd = 0.1) {
      p <- base
      for (nm in names(p)) {
        if (nm %in% rnames) next
        if (is.numeric(p[[nm]]) && length(p[[nm]]) == 1) {
          p[[nm]] <- p[[nm]] + rnorm(1, 0, abs(p[[nm]]) * sd + 0.01)
        }
      }
      return(p)
    }

    param_list <- c(list(params), lapply(seq_len(ncores - 1), function(i) jitter_params(params)))

    # Decide parallel backend
    use_mclapply <- .Platform$OS.type == "unix" && Sys.info()["sysname"] != "Darwin"

    opt_start <- proc.time()

    if (use_mclapply) {
      results <- parallel::mclapply(param_list, function(p) {
        tryCatch({
          res <- run_single_opt(p)
          list(obj = res$obj, opt = res$opt, obj_val = res$opt$objective)
        }, error = function(e) list(obj = NULL, opt = NULL, obj_val = Inf))
      }, mc.cores = ncores)
    } else {
      cl <- parallel::makeCluster(ncores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterExport(cl, c("tmb.data", "rnames", "default_map", "lower", "upper",
                                    "train_times", "alscl_dll_path"), envir = environment())
      parallel::clusterEvalQ(cl, { requireNamespace("TMB"); dyn.load(alscl_dll_path) })

      results <- parallel::parLapply(cl, param_list, function(p) {
        tryCatch({
          obj <- TMB::MakeADFun(tmb.data, p, random = rnames, map = default_map,
                                DLL = "ALSCL", inner.control = list(trace = FALSE, maxit = 500))
          opt <- nlminb(obj$par, obj$fn, obj$gr, lower = lower, upper = upper,
                        control = list(trace = 0, iter.max = 2000, eval.max = 10000))
          for (tt in seq_len(train_times - 1)) {
            opt <- nlminb(opt$par, obj$fn, obj$gr, lower = lower, upper = upper,
                          control = list(trace = 0, iter.max = 2000, eval.max = 10000))
          }
          list(obj = obj, opt = opt, obj_val = opt$objective)
        }, error = function(e) list(obj = NULL, opt = NULL, obj_val = Inf))
      })
    }

    opt_time <- (proc.time() - opt_start)[3]

    obj_vals <- sapply(results, function(r) r$obj_val)
    best_idx <- which.min(obj_vals)
    cat(sprintf("  Best start: #%d (obj = %.4f) | Total: %.1f sec\n",
                best_idx, obj_vals[best_idx], opt_time))

    obj <- results[[best_idx]]$obj
    opt <- results[[best_idx]]$opt

  } else {
    # ===== Single-start optimization =====
    cat("\nRunning ALSCL optimization with nlminb...\n")
    opt_start <- proc.time()

    single <- run_single_opt(params)
    obj <- single$obj
    opt <- single$opt

    opt_time <- (proc.time() - opt_start)[3]
    cat(sprintf("  Optimization done: %.1f sec (obj = %.4f)\n", opt_time, opt$objective))
  }

  # ===== sdreport =====
  cat("Running sdreport...\n")
  sd_start <- proc.time()

  sdresult <- tryCatch({
    TMB::sdreport(obj)
  }, error = function(e) {
    cat("  sdreport warning: ", e$message, "\n")
    tryCatch(TMB::sdreport(obj, getReportCovariance = FALSE),
             error = function(e2) { cat("  sdreport failed.\n"); NULL })
  })

  sd_time <- (proc.time() - sd_start)[3]
  cat(sprintf("  sdreport done: %.1f sec\n", sd_time))

  # ===== Extract results =====
  report <- obj$report()
  est_std <- if (!is.null(sdresult)) summary(sdresult) else NULL

  # Boundary check
  bound_check <- c((as.vector(opt$par) - as.vector(lower)),
                   (as.vector(upper) - as.vector(opt$par)))
  bound_hit <- min(bound_check) == 0
  final_outer_mgc <- max(abs(obj$gr(opt$par)))
  par_low_up <- cbind(opt$par, lower, upper)

  # Year labels
  cl_l <- tidyr::gather(data.CatL, key = "Year", value = "length", 2:ncol(data.CatL))
  year <- cl_l %>%
    dplyr::mutate(Year = as.numeric(gsub("X", "", Year))) %>%
    dplyr::distinct(Year) %>%
    dplyr::pull(Year)

  len_label <- data.CatL[, 1]

  total_time <- (proc.time() - total_start)[3]
  cat(sprintf("\n=== run_alscl total time: %.1f sec (%.1f min) ===\n", total_time, total_time / 60))

  # ===== Unload DLL =====
  dyn.unload(alscl_dll_path)

  # ===== Build result list =====
  result <- list(
    obj             = obj,
    opt             = opt,
    report          = report,
    est_std         = est_std,
    year            = year,
    len_mid         = len_mid,
    len_label       = len_label,
    len_border      = len_border,
    bound_hit       = bound_hit,
    bound_check     = bound_check,
    converge        = opt$message,
    final_outer_mgc = final_outer_mgc,
    par_low_up      = par_low_up,
    model_type      = "ALSCL",
    growth_step     = growth_step
  )

  # ===== Auto output =====
  if (output == TRUE) {
    dir.create("output/figures/result", recursive = TRUE, showWarnings = FALSE)
    dir.create("output/figures/diagnostic", recursive = TRUE, showWarnings = FALSE)
    dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

    save_fig <- function(fname, p) {
      tryCatch({
        print(p)
        ggsave(filename = fname, width = 16, height = 9, units = "in", dpi = 600)
      }, error = function(e) cat(sprintf("  Skipped %s: %s\n", fname, e$message)))
    }

    save_fig("output/figures/result/plot_recruitment.png",            plot_recruitment(result, se = TRUE))
    save_fig("output/figures/result/plot_SSB.png",                    plot_SSB(result, type = "SSB", se = TRUE))
    save_fig("output/figures/result/plot_SBL.png",                    plot_SSB(result, type = "SBL"))
    save_fig("output/figures/result/plot_biomass_B.png",              plot_biomass(result, type = "B", se = TRUE))
    save_fig("output/figures/result/plot_biomass_BL.png",             plot_biomass(result, type = "BL"))
    save_fig("output/figures/result/plot_abundance_N.png",            plot_abundance(result, type = "N", se = TRUE))
    save_fig("output/figures/result/plot_abundance_NA.png",           plot_abundance(result, type = "NA", se = TRUE))
    save_fig("output/figures/result/plot_abundance_NL.png",           plot_abundance(result, type = "NL"))
    save_fig("output/figures/result/plot_catch_CN.png",               plot_catch(result, type = "CN", se = TRUE))
    save_fig("output/figures/result/plot_catch_CNA.png",              plot_catch(result, type = "CNA"))
    save_fig("output/figures/result/plot_CatL_length.png",            plot_CatL(result, type = "length"))
    save_fig("output/figures/result/plot_CatL_Year(exp=T).png",       plot_CatL(result, type = "year", exp_transform = TRUE))
    save_fig("output/figures/result/plot_fishing_mortality_year.png",  plot_fishing_mortality(result, type = "year", se = TRUE))
    save_fig("output/figures/result/plot_fishing_mortality_age.png",   plot_fishing_mortality(result, type = "age", se = TRUE))
    save_fig("output/figures/result/plot_ridges.png",                 plot_ridges(result))
    save_fig("output/figures/result/plot_SSB_Rec.png",                plot_SSB_Rec(result))
    save_fig("output/figures/result/plot_VB.png",                     plot_VB(result, se = TRUE))
    save_fig("output/figures/result/plot_pla.png",                    plot_pla(result))
    save_fig("output/figures/diagnostic/plot_residuals_length.png",   plot_residuals(result, type = "length"))
    save_fig("output/figures/diagnostic/plot_residuals_year.png",     plot_residuals(result, type = "year"))

    # Tables
    tryCatch({
      diagnostics <- diagnose_model(data.CatL = data.CatL, model_result = result)
      write.csv(diagnostics, file = "output/tables/diagnostics.csv", row.names = FALSE)
    }, error = function(e) cat("  Skipped diagnostics table\n"))

    cat("Output saved to output/figures/ and output/tables/\n")
  }

  return(result)
}
