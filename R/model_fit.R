.acl_map <- function(parameters, map, model) {
  if (!is.null(map) && (!is.list(map) || (length(map) &&
      (is.null(names(map)) || any(!nzchar(names(map))) || anyDuplicated(names(map))))))
    stop("map must be a named list.")
  defaults <- if (model == "ACL") generate_map(map) else {
    out <- list(log_sigma_log_F = factor(NA), log_t0 = factor(NA))
    if (!is.null(map)) for (nm in names(map)) out[nm] <- map[nm]
    out[!vapply(out, is.null, logical(1))]
  }
  if (length(setdiff(names(defaults), names(parameters)))) stop("Unknown map parameter: ", paste(setdiff(names(defaults), names(parameters)), collapse = ", "))
  for (nm in names(defaults)) {
    value <- defaults[[nm]]
    if (length(value) == 1L && is.na(value)) defaults[[nm]] <- value <- factor(NA)
    if (!is.factor(value)) stop("map entries must be factors (factor(NA) fixes a parameter), or NULL to unfix.")
    if (length(value) != length(parameters[[nm]])) stop("Map length does not match parameter ", nm)
  }
  defaults
}
.acl_bounds <- function(par, config) {
  lower <- upper <- par
  for (nm in names(par)) {
    lo <- config$parameters.L[[nm]]; hi <- config$parameters.U[[nm]]
    if (is.null(lo)) lo <- -Inf
    if (is.null(hi)) hi <- Inf
    if (length(lo) != 1L || length(hi) != 1L || is.na(lo) || is.na(hi) || lo >= hi) stop("Invalid bounds for ", nm)
    lower[nm] <- lo; upper[nm] <- hi
  }
  if (any(!is.finite(par))) stop("Initial free parameters must be finite.")
  list(lower = lower, upper = upper, start = pmax(lower, pmin(upper, par)))
}
.acl_optimize <- function(obj, start, lower, upper, train_times, control) {
  for (i in seq_len(train_times)) {
    opt <- stats::nlminb(start, obj$fn, obj$gr, lower = lower, upper = upper, control = control)
    start <- opt$par
  }
  opt
}
.acl_fit_start <- function(i, data, params, mapping, random, dll, bounds,
                           train_times, control, obj = NULL) {
  started_at <- as.numeric(Sys.time())
  clock <- proc.time()
  # Identical starting points regardless of worker count, without changing the
  # caller's random-number stream during sequential multistart fits.
  if (i > 1L) {
    previous_kind <- RNGkind()
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) previous_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      do.call(RNGkind, as.list(previous_kind))
      if (had_seed) assign(".Random.seed", previous_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
  }
  start <- bounds$start
  if (i > 1L) {
    set.seed(137L * i, kind = "Mersenne-Twister", normal.kind = "Inversion")
    start <- start + stats::rnorm(length(start), 0, .15)
  }
  start <- pmax(bounds$lower, pmin(bounds$upper, start))
  opt <- tryCatch({
    if (is.null(obj)) obj <- TMB::MakeADFun(data, params, random = random,
      map = mapping, DLL = dll, inner.control = list(trace = FALSE, maxit = 500), silent = TRUE)
    .acl_optimize(obj, start, bounds$lower, bounds$upper, train_times, control)
  }, error = function(e) list(objective = Inf, message = conditionMessage(e), convergence = 1L))
  elapsed <- proc.time() - clock
  opt$start_id <- i
  opt$initial <- start
  opt$pid <- Sys.getpid()
  opt$started_at <- started_at
  opt$finished_at <- as.numeric(Sys.time())
  opt$elapsed <- unname(elapsed[["elapsed"]])
  opt$cpu_seconds <- unname(sum(elapsed[c("user.self", "sys.self")]))
  opt
}
.acl_fit <- function(model, prepared, parameters, parameters.L, parameters.U, map,
                     train_times, ncores, silent, output, control, nstarts = ncores) {
  train_times <- .acl_positive_integer(train_times, "train_times")
  ncores <- .acl_positive_integer(ncores, "ncores")
  nstarts <- .acl_positive_integer(nstarts, "nstarts")
  workers <- min(ncores, nstarts)
  started <- proc.time()[["elapsed"]]
  config <- create_parameters(model_type = tolower(model), parameters = parameters,
                              parameters.L = parameters.L, parameters.U = parameters.U)
  p <- config$parameters; d <- prepared$data
  if (model == "ALSCL") {
    d$len_mid <- prepared$bins$mid; d$len_lower <- prepared$bins$lower
    d$len_upper <- prepared$bins$upper; d$growth_step <- prepared$growth_step
    max_len <- max(prepared$bins$mid)
    if (is.null(parameters$log_Linf)) p$log_Linf <- log(1.1 * max_len)
    if (is.null(parameters.L$log_Linf)) config$parameters.L$log_Linf <- log(.7 * max_len)
    if (is.null(parameters.U$log_Linf)) config$parameters.U$log_Linf <- log(2 * max_len)
  }
  if (is.null(parameters$mean_log_R)) {
    totals <- colSums(prepared$observations, na.rm = TRUE)
    p$mean_log_R <- log(mean(totals[totals > 0]))
  }
  if (is.null(parameters$log_init_Z)) p$log_init_Z <- log(d$M + .3)
  # Length-at-age standard deviations must stay positive during optimization.
  t0_name <- if (model == "ACL") "t0" else "log_t0"
  t0_upper <- if (model == "ACL") min(d$age)*(1-1e-6) else log(min(d$age)*(1-1e-6))
  config$parameters.U[[t0_name]] <- min(config$parameters.U[[t0_name]], t0_upper)
  if ((if (model == "ACL") p$t0 else exp(p$log_t0)) >= min(d$age)) stop("t0 must be below the recruitment age.")
  p$dev_log_R <- rep(0, d$Y)
  p$dev_log_F <- array(0, c(if (model == "ACL") d$A else d$L, d$Y))
  p$dev_log_N0 <- rep(0, d$A-1L)
  mapping <- .acl_map(p, map, model)
  random <- c("dev_log_R", "dev_log_F", "dev_log_N0")
  # Clamp all starting values to the valid named bounds before building the tape.
  for (nm in names(config$parameters)) {
    if (is.null(mapping[[nm]]) || any(!is.na(mapping[[nm]])))
      p[[nm]] <- max(config$parameters.L[[nm]], min(config$parameters.U[[nm]], p[[nm]]))
  }
  if (silent) {
    invisible(utils::capture.output(info <- .acl_compile(model)))
  } else info <- .acl_compile(model)
  make_obj <- function() TMB::MakeADFun(d, p, random = random, map = mapping,
    DLL = info$dll_name, inner.control = list(trace = FALSE, maxit = 500), silent = TRUE)
  obj <- make_obj(); bounds <- .acl_bounds(obj$par, config)
  if (!length(obj$par)) stop("At least one fixed-effect parameter must remain free.")
  control <- utils::modifyList(list(iter.max = 2000L, eval.max = 10000L), control)
  if (!silent) message("Fitting ", model, " with ", nstarts, " start(s) on ", workers, " worker(s).")
  if (workers == 1L) {
    starts <- lapply(seq_len(nstarts), function(i) .acl_fit_start(i, d, p,
      mapping, random, info$dll_name, bounds, train_times, control,
      obj = if (nstarts == 1L) obj else NULL))
  } else {
    cl <- parallel::makeCluster(workers)
    on.exit(if (!is.null(cl)) parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, function(lib, path) {
      .libPaths(lib); loadNamespace("TMB"); loadNamespace("ALSCL"); dyn.load(path); NULL
    }, .libPaths(), info$dll_path)
    starts <- parallel::parLapply(cl, seq_len(nstarts), .acl_fit_start,
       data = d, params = p, mapping = mapping, random = random, dll = info$dll_name,
       bounds = bounds, train_times = train_times, control = control)
    parallel::stopCluster(cl)
    cl <- NULL
  }
  values <- vapply(starts, function(x) if (is.finite(x$objective)) x$objective else Inf, numeric(1))
  if (all(!is.finite(values))) stop("All optimization starts failed: ", paste(vapply(starts, `[[`, character(1), "message"), collapse = "; "))
  opt <- starts[[which.min(values)]]
  start_diagnostics <- do.call(rbind, lapply(starts, function(s) data.frame(
    start_id = s$start_id, pid = s$pid, started_at = s$started_at, finished_at = s$finished_at,
    elapsed = s$elapsed, cpu_seconds = s$cpu_seconds, objective = s$objective,
    convergence = s$convergence)))
  if (!is.finite(opt$objective)) stop("Optimization produced a non-finite objective.")
  # Re-evaluate in this process to synchronize random effects, reports and gradients.
  opt$objective <- obj$fn(opt$par)
  gradient <- as.numeric(obj$gr(opt$par)); names(gradient) <- names(opt$par)
  report <- obj$report()
  sdresult <- tryCatch(TMB::sdreport(obj, par.fixed = opt$par, getReportCovariance = FALSE),
    error = function(e) {warning("Standard errors unavailable: ", conditionMessage(e)); NULL})
  est_std <- if (is.null(sdresult)) NULL else summary(sdresult)
  covariance <- if (is.null(sdresult)) NULL else sdresult$cov.fixed
  if (!is.null(covariance)) dimnames(covariance) <- list(names(opt$par), names(opt$par))
  distances <- c(opt$par - bounds$lower, bounds$upper - opt$par)
  mgc <- max(abs(gradient))
  result <- list(obj = obj, opt = opt, report = report, est_std = est_std,
    vcov = covariance, pdHess = if (is.null(sdresult)) NA else isTRUE(sdresult$pdHess),
    gradient = gradient, max_abs_gradient = mgc, final_outer_mgc = mgc,
    bound_hit = any(distances <= 1e-7), bound_check = distances,
    converge = opt$message, convergence_code = opt$convergence,
    par_low_up = cbind(estimate = opt$par, lower = bounds$lower, upper = bounds$upper),
    model_type = model, year = prepared$year, age = d$age,
    len_mid = prepared$bins$mid, len_label = prepared$bins$labels,
    len_border = prepared$bins$border, len_lower = prepared$bins$lower, len_upper = prepared$bins$upper,
    growth_step = prepared$growth_step, elapsed = proc.time()[["elapsed"]]-started,
    starts = starts, start_diagnostics = start_diagnostics,
    nstarts = nstarts, workers = workers, observations = prepared$observations)
  if (output) .acl_save_output(result, prepared$data.CatL)
  result
}
.acl_save_output <- function(result, data.CatL) {
  dir.create("output/figures/result", recursive = TRUE, showWarnings = FALSE)
  dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
  plots <- list(recruitment = function() plot_recruitment(result), SSB = function() plot_SSB(result),
                biomass = function() plot_biomass(result), abundance = function() plot_abundance(result),
                catch = function() plot_catch(result), CatL = function() plot_CatL(result),
                fishing_mortality = function() plot_fishing_mortality(result),
                residuals = function() plot_residuals(result), VB = function() plot_VB(result))
  for (nm in names(plots)) tryCatch({
    ggplot2::ggsave(file.path("output/figures/result", paste0("plot_", nm, ".png")),
                   plot = plots[[nm]](), width = 12, height = 7, dpi = 300)
  }, error = function(e) warning("Could not save ", nm, ": ", conditionMessage(e)))
  utils::write.csv(diagnose_model(data.CatL, result), "output/tables/diagnostics.csv", row.names = FALSE)
  invisible(NULL)
}
