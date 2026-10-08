# Equal-work benchmark of the PUBLIC model-fitting functions.
# From the repository root, after R CMD INSTALL .:
# Rscript scripts/benchmark_parallel.R [output_dir] [repeats] [years]
# No optimizer limits, dummy work, sleeps, or discarded slow starts are used.
suppressPackageStartupMessages(library(ALSCL))
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1] else 'docs/results/parallel'
repeats <- if (length(args) > 1) as.integer(args[2]) else 3L
nyears <- if (length(args) > 2) as.integer(args[3]) else 10L
stopifnot(repeats >= 1L, nyears >= 2L, nyears <= 20L,
          'nstarts' %in% names(formals(run_acl)))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
cache <- Sys.getenv('ALSCL_BENCHMARK_CACHE', unset = file.path(tempdir(), 'benchmark-tmb'))
options(ALSCL.tmb.cache = cache)
x <- YTF_example
inputs <- lapply(x[c('data.CatL','data.wgt','data.mat')], function(d) d[, seq_len(nyears + 1L)])
workers <- c(1L, 2L, 4L)
nstarts <- 4L
train_times <- 2L
sources <- c('R/model_fit.R','R/run_acl.R','R/run_alscl.R',
             'inst/extdata/ACL.cpp','inst/extdata/ALSCL.cpp','inst/extdata/model_math.hpp',
             'data/YTF_example.rda','scripts/benchmark_parallel.R')
write.csv(data.frame(path = sources, md5 = unname(tools::md5sum(sources))),
          file.path(out, 'source_hashes.csv'), row.names = FALSE)
metadata <- c(paste('UTC:', format(Sys.time(), tz = 'UTC', usetz = TRUE)),
 paste('R:', R.version.string), paste('Platform:', R.version$platform),
 paste('ALSCL:', packageVersion('ALSCL')), paste('TMB:', packageVersion('TMB')),
 paste('RcppEigen:', packageVersion('RcppEigen')),
 paste('Logical CPUs:', parallel::detectCores()),
 paste('Data: synthetic YTF_example; seed 4;', nyears, 'years; 23 length bins; 15 ages'),
 paste('Work:', nstarts, 'fixed starts;', train_times, 'successive optimizer passes each'),
 paste('Repeats:', repeats), 'Warm compiled templates; compilation and warm-up excluded.',
 'Total fit timing includes preparation, worker startup/shutdown, fitting, final report and sdreport.',
 'Each start builds its own real TMB objective and calls nlminb with analytic gradients.',
 'Default optimizer controls and original conditional teaching maps; no early truncation.',
 'Other user analyses may share this computer. No unrelated processes are stopped.')
if (Sys.info()[['sysname']] == 'Darwin') metadata <- c(metadata,
 paste('CPU:', system2('sysctl', c('-n', 'machdep.cpu.brand_string'), stdout = TRUE)),
 paste('RAM bytes:', system2('sysctl', c('-n', 'hw.memsize'), stdout = TRUE)))
writeLines(metadata, file.path(out, 'environment.txt'))
# Avoid personal library paths in published metadata.
versions <- installed.packages()[, c('Package','Version')]
write.csv(versions, file.path(out, 'package_versions.csv'), row.names = FALSE)

peak_overlap <- function(d) {
  events <- rbind(data.frame(t = d$started_at, change = 1L),
                  data.frame(t = d$finished_at, change = -1L))
  max(cumsum(events$change[order(events$t, events$change)]))
}
all_runs <- all_starts <- all_parameters <- list()
references <- list()
write_results <- function() {
  write.csv(do.call(rbind, all_runs), file.path(out, 'runs.csv'), row.names = FALSE)
  write.csv(do.call(rbind, all_starts), file.path(out, 'starts.csv'), row.names = FALSE)
  write.csv(do.call(rbind, all_parameters), file.path(out, 'parameters.csv'), row.names = FALSE)
}
fit_model <- function(model, workers, starts) do.call(
  get(paste0('run_', tolower(model))),
  c(inputs, x$fit_args, x$fit_config[[tolower(model)]],
    list(ncores = workers, nstarts = starts, train_times = train_times, silent = TRUE)))

# Compile and execute one excluded warm-up for BOTH models before timed trials.
for (model in c('ACL','ALSCL')) {
  message('Warm-up ', model)
  warm <- fit_model(model, 1L, 1L)
  stopifnot(warm$convergence_code == 0L, isTRUE(warm$pdHess),
            warm$max_abs_gradient < .001)
  rm(warm); gc()
}
# Rotate the order to reduce systematic first/last-run and thermal effects.
for (rep in seq_len(repeats)) {
  order <- workers[((seq_along(workers) + rep - 2L) %% length(workers)) + 1L]
  for (model in if (rep %% 2L) c('ACL','ALSCL') else c('ALSCL','ACL')) {
    for (nc in order) {
      gc()
      run_id <- paste(model, rep, nc, sep = '-')
      message('START ', run_id, ' ', format(Sys.time()))
      before <- as.numeric(Sys.time())
      timed <- system.time(fit <- fit_model(model, nc, nstarts))
      d <- fit$start_diagnostics
      if (is.null(references[[model]])) references[[model]] <- list(
        objectives = d$objective, initial = lapply(fit$starts, `[[`, 'initial'),
        coefficients = lapply(fit$starts, `[[`, 'par'), B = fit$report$B)
      ref <- references[[model]]
      objective_diff <- max(abs(d$objective - ref$objectives))
      biomass_relative_diff <- max(abs(fit$report$B - ref$B) / pmax(1, abs(ref$B)))
      same_initial <- isTRUE(all.equal(lapply(fit$starts, `[[`, 'initial'), ref$initial, tolerance = 0))
      z <- data.frame(run_id, model, repetition = rep, ncores = nc, nstarts,
        elapsed_seconds = timed[['elapsed']], worker_cpu_seconds = sum(d$cpu_seconds),
        optimization_span_seconds = max(d$finished_at) - min(d$started_at),
        unique_pids = length(unique(d$pid)), peak_concurrent_starts = peak_overlap(d),
        objective = fit$opt$objective, convergence = fit$convergence_code,
        gradient = fit$max_abs_gradient, pdHess = fit$pdHess, bound_hit = fit$bound_hit,
        failed_starts = sum(!is.finite(d$objective) | d$convergence != 0L),
        max_objective_diff = objective_diff, max_biomass_relative_diff = biomass_relative_diff,
        same_initial, started_at = before)
      all_runs[[run_id]] <- z
      all_starts[[run_id]] <- cbind(run_id = run_id, model = model, ncores = nc, d)
      all_parameters[[run_id]] <- do.call(rbind, lapply(fit$starts, function(s)
        data.frame(run_id, start_id = s$start_id, parameter = names(s$initial),
                   initial = s$initial, estimate = s$par, row.names = NULL)))
      write_results()
      print(z)
      stopifnot(nrow(d) == nstarts, same_initial, objective_diff < 1e-6,
        biomass_relative_diff < 1e-5, all(is.finite(d$objective)), all(d$convergence == 0L),
        fit$convergence_code == 0L, isTRUE(fit$pdHess), fit$max_abs_gradient < .001,
        length(unique(d$pid)) == nc, all(d$cpu_seconds > 0),
        if (nc > 1L) peak_overlap(d) > 1L else all(d$pid == Sys.getpid()))
      rm(fit); gc()
    }
  }
}
runs <- do.call(rbind, all_runs)
summary <- do.call(rbind, lapply(split(runs, list(runs$model, runs$ncores), drop = TRUE), function(d)
  data.frame(model = d$model[1], ncores = d$ncores[1], nstarts = nstarts,
    repetitions = nrow(d), median_seconds = median(d$elapsed_seconds),
    min_seconds = min(d$elapsed_seconds), max_seconds = max(d$elapsed_seconds),
    median_cpu_seconds = median(d$worker_cpu_seconds),
    peak_concurrent_starts = max(d$peak_concurrent_starts))))
summary$speedup <- vapply(seq_len(nrow(summary)), function(i)
  summary$median_seconds[summary$model == summary$model[i] & summary$ncores == 1L] /
    summary$median_seconds[i], numeric(1))
summary <- summary[order(summary$model, summary$ncores), ]
write.csv(summary, file.path(out, 'summary.csv'), row.names = FALSE)
print(summary)
writeLines('PASS: all timed runs passed equality, convergence, PID and CPU checks.',
           file.path(out, 'validation.txt'))
