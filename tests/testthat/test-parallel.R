test_that("fixed multistart work is reproducible across socket worker counts", {
  for (model in c("ACL", "ALSCL")) {
    serial <- fit_fixture(model, ncores = 1, nstarts = 4)$result
    concurrent <- fit_fixture(model, ncores = 2, nstarts = 4)$result
    expect_length(serial$starts, 4)
    expect_length(concurrent$starts, 4)
    expect_equal(lapply(serial$starts, `[[`, "initial"),
                 lapply(concurrent$starts, `[[`, "initial"))
    expect_equal(serial$start_diagnostics$objective,
                 concurrent$start_diagnostics$objective, tolerance = 1e-8)
    expect_equal(serial$report$B, concurrent$report$B, tolerance = 1e-8)
    expect_equal(serial$opt$par, concurrent$opt$par, tolerance = 1e-8)
    expect_true(all(serial$start_diagnostics$convergence == 0))
    expect_true(all(concurrent$start_diagnostics$convergence == 0))
    expect_equal(unique(serial$start_diagnostics$pid), Sys.getpid())
    expect_length(unique(concurrent$start_diagnostics$pid), 2)
    expect_false(Sys.getpid() %in% concurrent$start_diagnostics$pid)
    expect_true(all(concurrent$start_diagnostics$cpu_seconds > 0))
    expect_true(all(concurrent$start_diagnostics$finished_at >=
                      concurrent$start_diagnostics$started_at))
    expect_equal(concurrent$workers, 2)
  }
})

test_that("worker cap and original ncores-only calls remain valid", {
  capped <- fit_fixture("ACL", ncores = 4, nstarts = 1)$result
  expect_equal(capped$workers, 1)
  expect_equal(capped$start_diagnostics$pid, Sys.getpid())
  x <- fit_fixture("ACL")
  legacy <- do.call(run_acl, c(x$input, list(parameters = x$parameters,
    map = x$map, ncores = 2, silent = TRUE)))
  expect_equal(legacy$nstarts, 2)
  expect_equal(legacy$workers, 2)
  expect_error(fit_fixture("ACL", nstarts = 0), "nstarts")
  expect_error(fit_fixture("ACL", nstarts = 1.5), "nstarts")
})

test_that("sequential jitter preserves the caller's RNG state and RNG kind", {
  original_kind <- RNGkind()
  on.exit(do.call(RNGkind, as.list(original_kind)), add = TRUE)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(912)
  before <- .Random.seed
  serial <- fit_fixture("ACL", ncores = 1, nstarts = 2)$result
  expect_identical(.Random.seed, before)
  expect_identical(RNGkind()[1], "L'Ecuyer-CMRG")
  concurrent <- fit_fixture("ACL", ncores = 2, nstarts = 2)$result
  expect_equal(lapply(serial$starts, `[[`, "initial"),
               lapply(concurrent$starts, `[[`, "initial"))
})

test_that("a failed start retains useful diagnostics", {
  bad <- ALSCL:::.acl_fit_start(1, list(), list(), list(), character(),
    "nonexistent_test_model", list(start = c(x = 1), lower = c(x = 0), upper = c(x = 2)),
    1, list())
  expect_identical(bad$objective, Inf)
  expect_equal(bad$convergence, 1)
  expect_equal(bad$pid, Sys.getpid())
  expect_true(nzchar(bad$message))
})

test_that("both retrospective models agree with 1, 2 and 4 requested workers", {
  for (model in c("ACL", "ALSCL")) {
    x <- fit_fixture(model)
    args <- c(list(nyear = 3, model_type = tolower(model)), x$input,
              list(parameters = x$parameters, map = x$map, silent = TRUE))
    serial <- do.call(retro_model, c(args, list(ncores = 1)))
    for (workers in c(2, 4)) {
      concurrent <- do.call(retro_model, c(args, list(ncores = workers)))
      expect_equal(concurrent$results, serial$results, tolerance = 1e-8)
      expect_equal(concurrent$rho_text, serial$rho_text, tolerance = 1e-8)
    }
  }
})

test_that("sim_acl fits actual replicates on the requested socket workers", {
  x <- fit_fixture("ACL"); f <- x$input; m <- x$result
  path <- tempfile("parallel-simulation-"); dir.create(path)
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  sim.data <- list(SN_at_len = t(as.matrix(f$data.CatL[-1])), nyear = length(m$year),
    nage = length(m$age), ages = m$age, growth_step = 1, len_mid = m$len_mid,
    len_border = m$len_border, weight = as.matrix(f$data.wgt[-1]),
    mat = as.matrix(f$data.mat[-1]), q_surv_L50 = f$sel_L50, q_surv_L95 = f$sel_L95)
  for (i in 1:4) save(sim.data, file = file.path(path, paste0("sim_rep", i)))
  for (workers in c(1, 2, 4)) {
    fits <- sim_acl(1:4, path, path, parameters = x$parameters, map = x$map,
                    M = f$M, ncores = workers)
    expect_true(all(vapply(fits, function(z) is.null(z$error), logical(1))))
    expect_equal(vapply(fits, function(z) z$opt$objective, numeric(1)),
                 setNames(rep(m$opt$objective, 4), names(fits)), tolerance = 1e-8)
    pids <- vapply(fits, function(z) z$start_diagnostics$pid, integer(1))
    expect_length(unique(pids), workers)
    if (workers > 1) expect_false(Sys.getpid() %in% pids)
  }
})
