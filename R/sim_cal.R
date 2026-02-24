#' Calculate Biological and Fishing Variables for Simulation
#'
#' This function calculates derived biological variables (length-at-age,
#' weight-at-length, maturity, selectivity, survey catchability) and the
#' appropriate transition matrix for the simulation model. For age-based
#' species (flatfish, krill) it builds the age-length transition matrix (pla);
#' for length-based species (tuna) it builds the growth transition matrix (Gij).
#'
#' @param params A list of parameters from \code{initialize_params()}.
#'   Must contain \code{model_type} ("age_based" or "length_based") to
#'   determine which transition matrix to build.
#'
#' @return A list containing:
#'   \describe{
#'     \item{model_type}{Character: "age_based" or "length_based"}
#'     \item{LatA}{Length-at-age vector}
#'     \item{W_at_len}{Weight-at-length vector}
#'     \item{mat}{Maturity-at-length vector}
#'     \item{q_surv}{Survey catchability-at-length vector}
#'     \item{sel}{Selectivity vector (at-age for age_based, at-length for length_based)}
#'     \item{pla}{Age-length transition matrix (always built for initial conditions)}
#'     \item{Gij}{Growth transition matrix (length_based only)}
#'     \item{nlen, len_lower, len_upper}{Length bin dimensions}
#'   }
#'
#' @export
sim_cal <- function(params) {

  # =========================================================================
  # Extract common parameters
  # =========================================================================
  nage       <- params$nage
  ages       <- params$ages
  rec.age    <- params$rec.age
  vbk        <- params$vbk
  Linf       <- params$Linf
  t0         <- params$t0
  len_mid    <- params$len_mid
  len_border <- params$len_border
  cv_L       <- params$cv_L
  cv_inc     <- params$cv_inc
  q_surv_L50 <- params$q_surv_L50
  q_surv_L95 <- params$q_surv_L95

  model_type <- if (!is.null(params$model_type)) params$model_type else "age_based"
  growth_step <- if (!is.null(params$growth_step)) params$growth_step else 1

  # =========================================================================
  # Length bin setup
  # =========================================================================
  nlen      <- length(len_mid)
  len_lower <- len_border[1:nlen]
  len_upper <- len_border[2:(nlen + 1)]

  # =========================================================================
  # Common derived variables
  # =========================================================================
  # Length-at-age (VB growth curve)
  LatA <- VB_func(Linf, vbk, t0, ages)

  # Weight-at-length: W = a * L^b
  W_at_len <- params$a * (len_mid ^ params$b)

  # Maturity-at-length: logistic function
  mat <- mat_func(params$mat_L50, params$mat_L95, len_mid)

  # Survey catchability-at-length: logistic function
  q_surv <- mat_func(q_surv_L50, q_surv_L95, len_mid)

  # =========================================================================
  # Age-length transition matrix (pla) -- always built, needed for
  # initial conditions in both age-based and length-based models
  # NOTE: In M7 (age-based), pla is nage x nlen (rows=ages, cols=lengths)
  #       In M6 (length-based), pla is nlen x nage (rows=lengths, cols=ages)
  # =========================================================================
  ml <- VB_func(Linf, vbk, t0, ages)
  sl <- cv_L * ml

  if (model_type == "age_based") {
    # pla[age, length] -- probability of length given age
    pla <- matrix(NA, nrow = nage, ncol = nlen)
    for (i in 1:nlen) {
      pla[, i] <- pnorm(len_border[i + 1], ml, sl) - pnorm(len_border[i], ml, sl)
    }
  } else {
    # pla[length, age] -- probability of length given age (transposed orientation)
    pla <- matrix(NA, nrow = nlen, ncol = nage)
    for (i in 1:nlen) {
      pla[i, ] <- pnorm(len_border[i + 1], ml, sl) - pnorm(len_border[i], ml, sl)
    }
  }

  # =========================================================================
  # Selectivity
  # =========================================================================
  sel_type <- if (!is.null(params$sel_type)) params$sel_type else "flat"

  if (sel_type == "flat") {
    # Flat selectivity at age (all ages fully selected)
    sel <- rep(1, nage)

  } else if (sel_type == "dome") {
    # Dome-shaped double logistic selectivity at length
    # Ascending limb (small fish becoming selected)
    L50_asc  <- params$sel_L50_asc
    L95_asc  <- params$sel_L95_asc
    # Descending limb (large fish becoming unselected; note L50 > L95)
    L50_desc <- params$sel_L50_desc
    L95_desc <- params$sel_L95_desc

    # Find the break point: where ascending > 0.5 and descending starts
    mid_len  <- (L95_asc + L50_desc) / 2
    asc_idx  <- which(len_mid <= mid_len)
    desc_idx <- which(len_mid > mid_len)

    sel <- numeric(nlen)
    sel[asc_idx]  <- mat_func(L50_asc, L95_asc, len_mid[asc_idx])
    sel[desc_idx] <- mat_func(L50_desc, L95_desc, len_mid[desc_idx])

  } else {
    # Default fallback: flat at age
    sel <- rep(1, nage)
  }

  # =========================================================================
  # Growth transition matrix Gij (length-based model only)
  # Gij[i,j] = probability of growing from bin j to bin i in one time step
  # Uses decreasing logistic growth increment model from M6 (tuna)
  # =========================================================================
  Gij <- NULL
  if (model_type == "length_based") {
    Gij <- matrix(NA, nrow = nlen, ncol = nlen)

    for (j in 1:nlen) {
      # Decreasing logistic growth increment
      # Maximum expected size increment (annual)
      delta_max <- (1 - exp(-vbk)) * Linf
      # Size at which growth increment is 50% of maximum
      l50_grow <- 0.5 * Linf
      # Size at which growth increment is 95% of maximum
      l95_grow <- 0.05 * Linf
      # Expected growth increment per time step (scaled by growth_step)
      ml_inc <- growth_step * delta_max /
        (1 + exp(-log(19) * (len_mid[j] - l50_grow) / (l95_grow - l50_grow)))
      # Variability of growth increment
      sl_inc <- ml_inc * cv_inc

      for (i in 1:nlen) {
        if (i < j) {
          # Cannot shrink: no transition to smaller bins
          Gij[i, j] <- 0
        } else if (i == j) {
          # Stay in same bin: CDF up to upper border
          Gij[i, j] <- pnorm(len_upper[i] - len_mid[j], ml_inc, sl_inc)
        } else {
          # Grow to larger bin: CDF difference
          Gij[i, j] <- pnorm(len_upper[i] - len_mid[j], ml_inc, sl_inc) -
            pnorm(len_lower[i] - len_mid[j], ml_inc, sl_inc)
        }
      }
    }
  }

  # =========================================================================
  # Return all calculated variables
  # =========================================================================
  list(
    model_type = model_type,
    LatA       = LatA,
    W_at_len   = W_at_len,
    mat        = mat,
    pla        = pla,
    Gij        = Gij,
    q_surv     = q_surv,
    sel        = sel,
    len_lower  = len_lower,
    len_upper  = len_upper,
    nlen       = nlen
  )
}
