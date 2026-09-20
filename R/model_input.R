# Shared input validation and bin parsing.
.acl_positive_integer <- function(x, name, minimum = 1L) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x != floor(x) || x < minimum || x > .Machine$integer.max)
    stop(name, " must be an integer >= ", minimum, call. = FALSE)
  as.integer(x)
}
.acl_scalar <- function(x, name, lower = -Inf, strict = FALSE) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
      (if (strict) x <= lower else x < lower))
    stop(name, " must be finite and ", if (strict) "> " else ">= ", lower, call. = FALSE)
  x
}
.acl_length_bins <- function(labels, len_mid = NULL, len_border = NULL,
                             len_lower = NULL, len_upper = NULL) {
  labels <- trimws(as.character(labels)); n <- length(labels)
  if (n < 2L || anyNA(labels) || anyDuplicated(labels)) stop("At least two distinct length bins are required.")
  number <- "[+-]?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?"
  single <- grepl(paste0("^", number, "$"), labels, perl = TRUE)
  lower <- upper <- mid <- rep(NA_real_, n)
  if (all(single)) {
    mid <- as.numeric(labels)
    edges <- (mid[-n] + mid[-1L]) / 2
    lower <- c(-Inf, edges); upper <- c(edges, Inf)
  } else {
    for (i in seq_len(n)) {
      s <- labels[i]
      if (grepl(paste0("^<\\s*", number, "$"), s, perl = TRUE) && i == 1L) {
        lower[i] <- -Inf; upper[i] <- as.numeric(sub("^<\\s*", "", s))
      } else if (grepl(paste0("^>\\s*", number, "$"), s, perl = TRUE) && i == n) {
        lower[i] <- as.numeric(sub("^>\\s*", "", s)); upper[i] <- Inf
      } else {
        pat <- paste0("^(", number, ")\\s*-\\s*(", number, ")$")
        parts <- regmatches(s, regexec(pat, s, perl = TRUE))[[1L]]
        if (length(parts) != 3L) {
          if (is.null(len_mid) || is.null(len_border)) stop("Cannot parse length bin '", s, "'; supply len_mid and len_border.")
        } else {
          lower[i] <- as.numeric(parts[2L]); upper[i] <- as.numeric(parts[3L])
          mid[i] <- (lower[i] + upper[i]) / 2
        }
      }
    }
    if (all(!is.na(lower)) && all(!is.na(upper))) {
      if (any(upper <= lower) || any(abs(upper[-n] - lower[-1L]) > 1e-8)) stop("Length intervals must have positive widths and be contiguous.")
      if (!is.finite(mid[1L]) && n > 2L) mid[1L] <- upper[1L] - (upper[2L] - lower[2L]) / 2
      if (!is.finite(mid[n]) && n > 2L) mid[n] <- lower[n] + (upper[n-1L] - lower[n-1L]) / 2
    }
    edges <- upper[-n]
  }
  if (!is.null(len_mid)) mid <- len_mid
  if (!is.null(len_border)) edges <- len_border
  if (!is.numeric(mid) || length(mid) != n || any(!is.finite(mid)) || any(diff(mid) <= 0)) stop("len_mid must contain one finite increasing midpoint per bin.")
  if (!is.numeric(edges) || length(edges) != n-1L || any(!is.finite(edges)) || any(diff(edges) <= 0)) stop("len_border must contain nbin - 1 finite increasing boundaries.")
  lower <- if (is.null(len_lower)) c(-Inf, edges) else len_lower
  upper <- if (is.null(len_upper)) c(edges, Inf) else len_upper
  if (!is.numeric(lower) || !is.numeric(upper) || length(lower) != n || length(upper) != n ||
      anyNA(lower) || anyNA(upper) || any(lower >= upper) || any(mid <= lower | mid >= upper) ||
      any(!is.finite(lower[-1L])) || any(!is.finite(upper[-n])) ||
      any(abs(lower[-1L] - edges) > 1e-8) || any(abs(upper[-n] - edges) > 1e-8))
    stop("Length bounds must contain the midpoints and agree with len_border.")
  list(mid = as.numeric(mid), border = as.numeric(edges), lower = c(-Inf, edges),
       upper = c(edges, Inf), labels = labels)
}
.acl_prepare_data <- function(data.CatL, data.wgt, data.mat, rec.age, nage, M,
                              sel_L50, sel_L95, growth_step = NULL,
                              len_mid = NULL, len_border = NULL, len_lower = NULL,
                              len_upper = NULL, zero_action = c("missing", "error")) {
  zero_action <- match.arg(zero_action)
  nage <- .acl_positive_integer(nage, "nage", 2L)
  rec.age <- .acl_scalar(rec.age, "rec.age", 0, TRUE); M <- .acl_scalar(M, "M", 0)
  if (is.null(growth_step)) growth_step <- if (rec.age < 1) rec.age else 1
  growth_step <- .acl_scalar(growth_step, "growth_step", 0, TRUE)
  .acl_scalar(sel_L50, "sel_L50"); .acl_scalar(sel_L95, "sel_L95")
  if (sel_L95 <= sel_L50) stop("sel_L95 must be greater than sel_L50.")
  frames <- list(data.CatL = data.CatL, data.wgt = data.wgt, data.mat = data.mat)
  if (any(vapply(frames, function(x) length(dim(x)) != 2L || ncol(x) < 3L || nrow(x) < 2L, logical(1)))) stop("Inputs require a length-label column, at least two time columns and two length bins.")
  frames <- lapply(frames, as.data.frame, stringsAsFactors = FALSE); catl <- frames$data.CatL
  if (anyDuplicated(names(catl)[-1L])) stop("Time column names must be unique.")
  for (nm in names(frames)[-1L]) {
    x <- frames[[nm]]
    if (!identical(dim(x), dim(catl)) || !identical(as.character(x[[1L]]), as.character(catl[[1L]])) || !identical(names(x)[-1L], names(catl)[-1L])) stop(nm, " must match data.CatL length labels and time columns, including order.")
  }
  numeric_matrix <- function(x, name) {
    original <- as.matrix(x[-1L])
    z <- suppressWarnings(matrix(as.numeric(original), nrow = nrow(x), dimnames = list(NULL, names(x)[-1L])))
    if (any(is.na(z) & !is.na(original))) stop(name, " contains nonnumeric observations.")
    z
  }
  obs <- numeric_matrix(catl, "data.CatL"); weight <- numeric_matrix(frames$data.wgt, "data.wgt"); mat <- numeric_matrix(frames$data.mat, "data.mat")
  if (any(is.infinite(obs)) || any(obs < 0, na.rm = TRUE)) stop("Catch must be nonnegative and finite, or NA for missing observations.")
  if (zero_action == "error" && any(obs == 0, na.rm = TRUE)) stop("Zero observations are incompatible with the lognormal observation model.")
  if (any(!is.finite(weight)) || any(weight < 0)) stop("Weights must be finite and nonnegative.")
  if (any(!is.finite(mat)) || any(mat < 0 | mat > 1)) stop("Maturity must be finite and in [0, 1].")
  observed <- !is.na(obs) & obs > 0
  if (!any(observed)) stop("At least one positive survey observation is required.")
  log_obs <- matrix(0, nrow(obs), ncol(obs), dimnames = dimnames(obs)); log_obs[observed] <- log(obs[observed])
  bins <- .acl_length_bins(catl[[1L]], len_mid, len_border, len_lower, len_upper)
  years <- suppressWarnings(as.numeric(sub("^X", "", names(catl)[-1L])))
  if (any(!is.finite(years)) || any(diff(years) <= 0)) stop("Time columns must have increasing numeric labels (decimal years are supported).")
  log_q <- stats::plogis(log(19) * (bins$mid - sel_L50) / (sel_L95 - sel_L50), log.p = TRUE)
  list(data = list(logN_at_len = log_obs, na_matrix = observed * 1, log_q = log_q,
                  len_border = bins$border, age = rec.age + seq.int(0L, nage-1L) * growth_step,
                  Y = ncol(obs), A = nage, L = nrow(obs), weight = weight, mat = mat, M = M),
       bins = bins, year = years, growth_step = growth_step, observations = obs, data.CatL = catl)
}
.acl_observed_log <- function(model_result) {
  d <- model_result$obj$env$.data; z <- d$logN_at_len
  if (!is.null(d$na_matrix)) z[d$na_matrix == 0] <- NA_real_
  z
}
