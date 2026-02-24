#' Simulate Fishery Dynamics
#'
#' This function simulates fish population dynamics using either an age-based
#' (flatfish/krill) or length-based (tuna) operating model. The model type is
#' determined automatically from the \code{bio_vars$model_type} field set by
#' \code{sim_cal()}.
#'
#' \strong{Age-based pathway} (M7 flatfish, krill): Standard forward-projection
#' of \code{N_at_age} with age-length transition matrix (pla). F is age-dependent.
#'
#' \strong{Length-based pathway} (M6 tuna): Three-dimensional tracking of
#' \code{nla[length, age]} with growth transition matrix (Gij). F is
#' length-dependent. Supports quarterly (sub-annual) time steps.
#'
#' @param bio_vars A list from \code{sim_cal()} containing biological variables.
#' @param params A list from \code{initialize_params()} containing parameters.
#' @param sim_year Integer, total number of simulation years (default: 100).
#' @param output_dir Character, directory for saving results (default: tempdir()).
#' @param iter_range Integer vector, which iterations (seeds) to run (default: 4:100).
#' @param return_iter Integer or NULL, which iteration to return in memory.
#'
#' @return A list containing simulated fishery data for the observation window
#'   (post-burn-in), including SN_at_len, N_at_len, N_at_age, SSB, etc.
#'
#' @export
sim_data <- function(bio_vars, params, sim_year = 100,
                     output_dir = ".", iter_range = 4:100, return_iter = NULL) {

  # Determine model type from bio_vars (set by sim_cal)
  model_type <- if (!is.null(bio_vars$model_type)) bio_vars$model_type else "age_based"

  if (model_type == "length_based") {
    .sim_data_length_based(bio_vars, params, sim_year, output_dir, iter_range, return_iter)
  } else {
    .sim_data_age_based(bio_vars, params, sim_year, output_dir, iter_range, return_iter)
  }
}


# ===========================================================================
# Age-based simulation engine (M7 flatfish / krill)
# Source: M7.R operating model
# ===========================================================================
.sim_data_age_based <- function(bio_vars, params, sim_year, output_dir,
                                iter_range, return_iter) {

  # --- Extract variables ---
  nlen    <- bio_vars$nlen
  W_at_len <- bio_vars$W_at_len
  mat     <- bio_vars$mat
  pla     <- bio_vars$pla      # nage x nlen
  q_surv  <- bio_vars$q_surv
  sel     <- bio_vars$sel      # at-age (length nage)

  nage    <- params$nage
  M       <- params$M
  init_Z  <- params$init_Z
  len_mid <- params$len_mid
  std_logR  <- params$std_logR
  std_logN0 <- params$std_logN0
  std_SN    <- params$std_SN
  alpha   <- params$alpha
  beta    <- params$beta
  R_init  <- if (!is.null(params$R_init)) params$R_init else 500
  F_mean  <- if (!is.null(params$F_mean)) params$F_mean else 0.3
  F_ar    <- if (!is.null(params$F_ar))   params$F_ar   else 0.75
  F_sd    <- if (!is.null(params$F_sd))   params$F_sd   else 0.2
  R_ar    <- if (!is.null(params$R_ar))   params$R_ar   else 0.1
  burn_in <- if (!is.null(params$burn_in)) params$burn_in else 80

  # Observation window indices
  obs_start <- burn_in + 1
  obs_end   <- sim_year
  nyear_obs <- obs_end - obs_start + 1

  # Determine return iteration
  if (is.null(return_iter)) return_iter <- max(iter_range)
  if (!(return_iter %in% iter_range))
    stop("return_iter = ", return_iter, " is not in iter_range")

  return_data <- NULL

  for (iter in iter_range) {

    # --- Recruitment deviations (AR1) ---
    set.seed(iter)
    dev_logR <- arima.sim(list(order = c(1, 0, 0), ar = R_ar), n = sim_year) * std_logR
    Rec <- exp(log(R_init) + dev_logR)

    # --- Fishing mortality (AR1) ---
    set.seed(iter)
    F_yr <- F_mean * exp(arima.sim(list(order = c(1, 0, 0), ar = F_ar), n = sim_year) * F_sd)

    Z_at_age <- matrix(NA, nrow = sim_year, ncol = nage)
    F_at_age <- matrix(NA, nrow = sim_year, ncol = nage)
    M_at_age <- matrix(M, nrow = sim_year, ncol = nage)
    for (i in 1:sim_year) {
      F_at_age[i, ] <- F_yr[i] * sel
      Z_at_age[i, ] <- F_at_age[i, ] + M_at_age[i, ]
    }

    # --- Initial age structure ---
    set.seed(iter)
    N0_at_age <- rep(NA, nage)
    dev_logN0 <- rnorm(nage - 1, 0, std_logN0)
    N0_at_age[1] <- Rec[1]
    for (i in 2:nage) {
      N0_at_age[i] <- N0_at_age[i - 1] * exp(-init_Z + dev_logN0[i - 1])
    }

    # --- Pre-allocate ---
    N_at_age <- matrix(NA, nrow = sim_year, ncol = nage)
    N_at_len <- matrix(NA, nrow = sim_year, ncol = nlen)
    B_at_len <- matrix(NA, nrow = sim_year, ncol = nlen)
    SB_at_len <- matrix(NA, nrow = sim_year, ncol = nlen)
    CN_at_age <- matrix(NA, nrow = sim_year, ncol = nage)
    CN_at_len <- matrix(NA, nrow = sim_year, ncol = nlen)
    CB_at_len <- matrix(NA, nrow = sim_year, ncol = nlen)
    SSB <- TN <- TB <- CN_tot <- CB_tot <- rep(NA, sim_year)

    # --- Year 1: initialization ---
    N_at_age[1, ] <- N0_at_age
    N_at_len[1, ] <- N_at_age[1, ] %*% pla
    B_at_len[1, ] <- N_at_len[1, ] * W_at_len
    SB_at_len[1, ] <- B_at_len[1, ] * mat
    SSB[1] <- sum(SB_at_len[1, ])

    # --- Forward dynamics ---
    for (i in 2:sim_year) {
      # Beverton-Holt recruitment with AR1 error
      logR <- log(alpha * SSB[i - 1] / (beta + SSB[i - 1]))
      Rec[i] <- exp(logR + dev_logR[i])
      N_at_age[i, 1] <- Rec[i]

      # Survival (ages 2 to nage-1)
      for (j in 2:(nage - 1)) {
        N_at_age[i, j] <- N_at_age[i - 1, j - 1] * exp(-Z_at_age[i - 1, j - 1])
      }
      # Plus group
      N_at_age[i, nage] <- N_at_age[i - 1, nage - 1] * exp(-Z_at_age[i - 1, nage - 1]) +
        N_at_age[i - 1, nage] * exp(-Z_at_age[i - 1, nage])

      # Age -> Length conversion
      N_at_len[i, ] <- N_at_age[i, ] %*% pla
      B_at_len[i, ] <- N_at_len[i, ] * W_at_len
      SB_at_len[i, ] <- B_at_len[i, ] * mat
      SSB[i] <- sum(SB_at_len[i, ])
    }

    # --- Derived quantities ---
    for (i in 1:sim_year) {
      for (j in 1:nage) {
        CN_at_age[i, j] <- N_at_age[i, j] * (1 - exp(-Z_at_age[i, j])) *
          (F_at_age[i, j] / Z_at_age[i, j])
      }
      CN_at_len[i, ] <- CN_at_age[i, ] %*% pla
      CB_at_len[i, ] <- CN_at_len[i, ] * W_at_len
      TN[i] <- sum(N_at_len[i, ])
      TB[i] <- sum(B_at_len[i, ])
      CN_tot[i] <- sum(CN_at_len[i, ])
      CB_tot[i] <- sum(CB_at_len[i, ])
    }

    # --- Survey index ---
    RVN_at_len <- matrix(NA, nrow = sim_year, ncol = nlen)
    RVB_at_len <- matrix(NA, nrow = sim_year, ncol = nlen)
    set.seed(iter)
    surv_error <- rnorm(sim_year, 0, std_SN)
    for (i in 1:sim_year) {
      RVN_at_len[i, ] <- N_at_len[i, ] * q_surv * exp(surv_error[i])
      RVB_at_len[i, ] <- RVN_at_len[i, ] * W_at_len
    }
    # Floor small values
    RVN_at_len[RVN_at_len < 1e-6] <- 1e-6

    # --- Extract observation window ---
    idx <- obs_start:obs_end

    sim.data <- list(
      SN_at_len  = RVN_at_len[idx, ],
      q_surv     = q_surv,
      len_mid    = len_mid,
      nyear      = nyear_obs,
      nage       = nage,
      nlen       = nlen,
      ages       = params$ages,
      weight     = matrix(rep(W_at_len, nyear_obs), nrow = nlen, ncol = nyear_obs),
      mat        = matrix(rep(mat, nyear_obs), nrow = nlen, ncol = nyear_obs),
      sel        = sel,
      F_at_age   = F_at_age[idx, ],
      M_at_age   = M_at_age[idx, ],
      N_at_len   = N_at_len[idx, ],
      N_at_age   = N_at_age[idx, ],
      B_at_len   = B_at_len[idx, ],
      SB_at_len  = SB_at_len[idx, ],
      CN_at_len  = CN_at_len[idx, ],
      CN_at_age  = CN_at_age[idx, ],
      CB_at_len  = CB_at_len[idx, ],
      RVN_at_len = RVN_at_len[idx, ],
      RVB_at_len = RVB_at_len[idx, ],
      Rec        = N_at_age[idx, 1],
      TN         = TN[idx],
      TB         = TB[idx],
      CN         = CN_tot[idx],
      CB         = CB_tot[idx],
      SSB        = SSB[idx]
    )

    save(sim.data, file = file.path(output_dir, paste0("sim_rep", iter)))
    if (iter == return_iter) {
      return_data <- sim.data
      return_data$iter <- iter
    }
  }

  cat("Simulated iterations:", min(iter_range), "to", max(iter_range),
      "(", length(iter_range), "replicates )\n")
  cat("Returned iteration:", return_iter, "\n")
  return(return_data)
}


# ===========================================================================
# Length-based simulation engine (M6 tuna)
# Source: ALSCL_simulation.R (M6 operating model)
# Three-dimensional tracking: nla[nlen, nage] at each time step
# Growth via Gij transition matrix; F is length-dependent
# ===========================================================================
.sim_data_length_based <- function(bio_vars, params, sim_year, output_dir,
                                   iter_range, return_iter) {

  # --- Extract variables ---
  nlen     <- bio_vars$nlen
  W_at_len <- bio_vars$W_at_len
  mat      <- bio_vars$mat
  pla      <- bio_vars$pla      # nlen x nage (transposed orientation for length-based)
  Gij      <- bio_vars$Gij      # nlen x nlen growth transition matrix
  q_surv   <- bio_vars$q_surv
  sel      <- bio_vars$sel      # at-length (length nlen)
  len_lower <- bio_vars$len_lower
  len_upper <- bio_vars$len_upper

  nage      <- params$nage
  M         <- params$M
  init_Z    <- params$init_Z
  len_mid   <- params$len_mid
  std_logR  <- params$std_logR
  std_logN0 <- params$std_logN0
  std_SN    <- params$std_SN
  alpha     <- params$alpha
  beta      <- params$beta
  R_init    <- if (!is.null(params$R_init)) params$R_init else 10000
  F_mean    <- if (!is.null(params$F_mean)) params$F_mean else 0.2
  F_ar      <- if (!is.null(params$F_ar))   params$F_ar   else 0.75
  F_sd      <- if (!is.null(params$F_sd))   params$F_sd   else 0.2
  R_ar      <- if (!is.null(params$R_ar))   params$R_ar   else 0.1
  burn_in   <- if (!is.null(params$burn_in)) params$burn_in else 95
  growth_step <- if (!is.null(params$growth_step)) params$growth_step else 0.25

  # Quarterly time steps: total number of steps
  steps_per_year <- round(1 / growth_step)
  nstep <- sim_year * steps_per_year

  # Observation window (in quarterly steps)
  obs_start <- burn_in * steps_per_year + 1
  obs_end   <- nstep
  nstep_obs <- obs_end - obs_start + 1

  # Determine return iteration
  if (is.null(return_iter)) return_iter <- max(iter_range)
  if (!(return_iter %in% iter_range))
    stop("return_iter = ", return_iter, " is not in iter_range")

  return_data <- NULL

  for (iter in iter_range) {

    # --- Recruitment deviations (AR1, at quarterly resolution) ---
    set.seed(iter)
    dev_logR <- arima.sim(list(order = c(1, 0, 0), ar = R_ar), n = nstep) * std_logR
    Rec <- exp(log(R_init) + dev_logR)

    # --- Length-dependent fishing mortality (AR1) ---
    set.seed(iter)
    F_yr <- F_mean * exp(arima.sim(list(order = c(1, 0, 0), ar = F_ar), n = nstep) * F_sd)

    Z_at_len <- matrix(NA, nrow = nstep, ncol = nlen)
    F_at_len <- matrix(NA, nrow = nstep, ncol = nlen)
    M_at_len <- matrix(M, nrow = nstep, ncol = nlen)
    for (i in 1:nstep) {
      F_at_len[i, ] <- F_yr[i] * sel
      Z_at_len[i, ] <- F_at_len[i, ] + M_at_len[i, ]
    }

    # --- Initial age structure ---
    set.seed(iter)
    N0_at_age <- rep(NA, nage)
    dev_logN0 <- rnorm(nage - 1, 0, std_logN0)
    N0_at_age[1] <- Rec[1]
    for (i in 2:nage) {
      N0_at_age[i] <- N0_at_age[i - 1] * exp(-init_Z + dev_logN0[i - 1])
    }

    # --- 3D arrays: nla[nlen, nage] at each step ---
    nla_list  <- vector("list", nstep)
    bla_list  <- vector("list", nstep)
    sbla_list <- vector("list", nstep)
    cnla_list <- vector("list", nstep)
    cbla_list <- vector("list", nstep)

    nla  <- matrix(NA, nrow = nlen, ncol = nage)
    bla  <- matrix(NA, nrow = nlen, ncol = nage)
    sbla <- matrix(NA, nrow = nlen, ncol = nage)
    cnla <- matrix(NA, nrow = nlen, ncol = nage)
    cbla <- matrix(NA, nrow = nlen, ncol = nage)

    SSB <- TN <- TB <- CN_tot <- CB_tot <- rep(NA, nstep)

    # --- Step 1: Initialize from age structure ---
    for (a in 1:nage) {
      nla[, a]  <- N0_at_age[a] * pla[, a]
      bla[, a]  <- nla[, a] * W_at_len
      sbla[, a] <- bla[, a] * mat
      cnla[, a] <- nla[, a] * (1 - exp(-Z_at_len[1, ])) * (F_at_len[1, ] / Z_at_len[1, ])
      cbla[, a] <- cnla[, a] * W_at_len
    }
    nla_list[[1]]  <- nla
    bla_list[[1]]  <- bla
    sbla_list[[1]] <- sbla
    cnla_list[[1]] <- cnla
    cbla_list[[1]] <- cbla
    TN[1]  <- sum(nla)
    TB[1]  <- sum(bla)
    SSB[1] <- sum(sbla)
    CN_tot[1] <- sum(cnla)
    CB_tot[1] <- sum(cbla)

    # --- Forward dynamics: 3D length-age tracking with growth transition ---
    for (i in 2:nstep) {
      # Recruitment: age-1 group distributed across length bins
      logR <- log(alpha * SSB[i - 1] / (beta + SSB[i - 1]))
      Rec[i] <- exp(logR + dev_logR[i])
      nla[, 1]  <- Rec[i] * pla[, 1]
      bla[, 1]  <- nla[, 1] * W_at_len
      sbla[, 1] <- bla[, 1] * mat
      cnla[, 1] <- nla[, 1] * (1 - exp(-Z_at_len[i, ])) * (F_at_len[i, ] / Z_at_len[i, ])
      cbla[, 1] <- cnla[, 1] * W_at_len

      # Ages 2 to nage-1: survive + grow via Gij
      for (j in 2:(nage - 1)) {
        nla_survival <- nla_list[[i - 1]][, j - 1] * exp(-Z_at_len[i - 1, ])
        nla[, j]  <- Gij %*% nla_survival
        bla[, j]  <- nla[, j] * W_at_len
        sbla[, j] <- bla[, j] * mat
        cnla[, j] <- nla[, j] * (1 - exp(-Z_at_len[i, ])) * (F_at_len[i, ] / Z_at_len[i, ])
        cbla[, j] <- cnla[, j] * W_at_len
      }

      # Plus group: combine last two ages, survive + grow
      nla_survival <- nla_list[[i - 1]][, nage - 1] * exp(-Z_at_len[i - 1, ]) +
        nla_list[[i - 1]][, nage] * exp(-Z_at_len[i - 1, ])
      nla[, nage]  <- Gij %*% nla_survival
      bla[, nage]  <- nla[, nage] * W_at_len
      sbla[, nage] <- bla[, nage] * mat
      cnla[, nage] <- nla[, nage] * (1 - exp(-Z_at_len[i, ])) * (F_at_len[i, ] / Z_at_len[i, ])
      cbla[, nage] <- cnla[, nage] * W_at_len

      nla_list[[i]]  <- nla
      bla_list[[i]]  <- bla
      sbla_list[[i]] <- sbla
      cnla_list[[i]] <- cnla
      cbla_list[[i]] <- cbla

      TN[i]     <- sum(nla)
      TB[i]     <- sum(bla)
      SSB[i]    <- sum(sbla)
      CN_tot[i] <- sum(cnla)
      CB_tot[i] <- sum(cbla)
    }

    # --- Aggregate 3D -> 2D arrays ---
    N_at_age  <- matrix(NA, nrow = nstep, ncol = nage)
    N_at_len  <- matrix(NA, nrow = nstep, ncol = nlen)
    B_at_len  <- matrix(NA, nrow = nstep, ncol = nlen)
    SB_at_len <- matrix(NA, nrow = nstep, ncol = nlen)
    CN_at_len <- matrix(NA, nrow = nstep, ncol = nlen)
    CB_at_len <- matrix(NA, nrow = nstep, ncol = nlen)
    CN_at_age <- matrix(NA, nrow = nstep, ncol = nage)

    for (i in 1:nstep) {
      N_at_age[i, ]  <- colSums(nla_list[[i]])
      N_at_len[i, ]  <- rowSums(nla_list[[i]])
      B_at_len[i, ]  <- rowSums(bla_list[[i]])
      SB_at_len[i, ] <- rowSums(sbla_list[[i]])
      CN_at_len[i, ] <- rowSums(cnla_list[[i]])
      CB_at_len[i, ] <- rowSums(cbla_list[[i]])
      CN_at_age[i, ] <- colSums(cnla_list[[i]])
    }

    # --- Survey index (3D-aware, like M6 original) ---
    RVN_at_len <- matrix(NA, nrow = nstep, ncol = nlen)
    RVB_at_len <- matrix(NA, nrow = nstep, ncol = nlen)
    set.seed(iter)
    surv_error <- rnorm(nstep, 0, std_SN)
    for (i in 1:nstep) {
      # Apply q_surv to each age column, then sum across ages
      rvnla <- nla_list[[i]] * q_surv * exp(surv_error[i])
      RVN_at_len[i, ] <- rowSums(rvnla)
      RVB_at_len[i, ] <- RVN_at_len[i, ] * W_at_len
    }
    # Floor small values
    RVN_at_len[RVN_at_len < 1e-6] <- 1e-6

    # --- Extract observation window ---
    idx <- obs_start:obs_end

    sim.data <- list(
      SN_at_len  = RVN_at_len[idx, ],
      q_surv     = q_surv,
      len_mid    = len_mid,
      nyear      = nstep_obs,              # number of time steps in obs window
      nage       = nage,
      nlen       = nlen,
      ages       = params$ages,
      weight     = matrix(rep(W_at_len, nstep_obs), nrow = nlen, ncol = nstep_obs),
      mat        = matrix(rep(mat, nstep_obs), nrow = nlen, ncol = nstep_obs),
      sel        = sel,
      F_at_len   = F_at_len[idx, ],
      M_at_len   = M_at_len[idx, ],
      N_at_len   = N_at_len[idx, ],
      N_at_age   = N_at_age[idx, ],
      B_at_len   = B_at_len[idx, ],
      SB_at_len  = SB_at_len[idx, ],
      CN_at_len  = CN_at_len[idx, ],
      CN_at_age  = CN_at_age[idx, ],
      CB_at_len  = CB_at_len[idx, ],
      RVN_at_len = RVN_at_len[idx, ],
      RVB_at_len = RVB_at_len[idx, ],
      Rec        = N_at_age[idx, 1],
      TN         = TN[idx],
      TB         = TB[idx],
      CN         = CN_tot[idx],
      CB         = CB_tot[idx],
      SSB        = SSB[idx]
    )

    save(sim.data, file = file.path(output_dir, paste0("sim_rep", iter)))
    if (iter == return_iter) {
      return_data <- sim.data
      return_data$iter <- iter
    }
  }

  cat("Simulated iterations:", min(iter_range), "to", max(iter_range),
      "(", length(iter_range), "replicates )\n")
  cat("Returned iteration:", return_iter, "\n")
  return(return_data)
}
