#' Initialize Parameters for the Fish Population Simulation
#'
#' Initializes all biological, fishing, and survey parameters needed for the
#' population simulation via \code{sim_cal()} and \code{sim_data()}. Users can
#' select a built-in species preset (\code{"flatfish"}, \code{"tuna"}, or
#' \code{"krill"}) to get biologically consistent parameter sets, or specify
#' all parameters manually.
#'
#' @param species Character or NULL. Species preset name, one of
#'   \code{"flatfish"}, \code{"tuna"}, \code{"krill"}, or \code{"?"} to list
#'   presets. NULL (default) uses generic defaults equivalent to flatfish.
#' @param nyear Integer. Total simulation years including burn-in (default 100).
#' @param rec.age Numeric. Age of recruitment (default 1; use 0.25 for quarterly).
#' @param first.year Integer. First year of the observation window (default 2020).
#' @param nage Integer. Number of age classes including plus-group (default 15).
#' @param M Numeric. Natural mortality rate (default 0.2).
#' @param init_Z Numeric. Initial total mortality for equilibrium (default 0.5).
#' @param vbk Numeric. Von Bertalanffy k (default 0.2).
#' @param Linf Numeric. Asymptotic length (default 60).
#' @param t0 Numeric. Von Bertalanffy t0 (default 1/60).
#' @param len_mid Numeric vector. Midpoints of length bins.
#' @param len_border Numeric vector. Borders of length bins.
#' @param a Numeric. Length-weight coefficient (default exp(-12)).
#' @param b Numeric. Length-weight exponent (default 3).
#' @param mat_L50 Numeric. Length at 50 pct maturity (default 35).
#' @param mat_L95 Numeric. Length at 95 pct maturity (default 40).
#' @param cv_L Numeric. CV of length-at-age (default 0.2).
#' @param cv_inc Numeric. CV of growth increment (default 0.2).
#' @param std_logR Numeric. SD of log-recruitment deviations (default 0.3).
#' @param std_logN0 Numeric. SD of initial log-N deviations (default 0.2).
#' @param alpha Numeric. Beverton-Holt alpha (default 400).
#' @param beta Numeric. Beverton-Holt beta (default 10).
#' @param std_SN Numeric. Survey observation error SD (default 0.2).
#' @param q_surv_L50 Numeric. Survey selectivity L50 (default 15).
#' @param q_surv_L95 Numeric. Survey selectivity L95 (default 20).
#' @param R_init Numeric. Initial recruitment (default 500).
#' @param F_mean Numeric. Mean fishing mortality (default 0.3).
#' @param F_ar Numeric. AR1 coefficient for F deviations (default 0.75).
#' @param F_sd Numeric. SD of F deviations (default 0.2).
#' @param R_ar Numeric. AR1 coefficient for recruitment deviations (default 0.1).
#'
#' @return A list containing all initialized parameters and derived variables
#'   (ages, years, etc.), ready to pass to \code{sim_cal()}.
#'
#' @details
#' Species presets and their key parameters:
#' \itemize{
#'   \item flatfish: nage=15, Linf=60, vbk=0.2, M=0.2, annual (M7 model)
#'   \item tuna: nage=20, Linf=152, vbk=0.38, M=0.2, quarterly (M6 model)
#'   \item krill: nage=5, Linf=60, vbk=0.3, M=0.3, annual
#' }
#'
#' @examples
#' # List available presets
#' initialize_params(species = "?")
#'
#' # M7 flatfish (same as default)
#' params <- initialize_params(species = "flatfish")
#'
#' # M6 tuna with quarterly time steps
#' params <- initialize_params(species = "tuna")
#'
#' # Krill with custom M
#' params <- initialize_params(species = "krill", M = 0.4)
#'
#' # Fully custom (no preset)
#' params <- initialize_params(nage = 10, Linf = 80, vbk = 0.15)
#'
#' @export
initialize_params <- function(species    = NULL,
                              nyear      = NULL,
                              rec.age    = NULL,
                              first.year = NULL,
                              nage       = NULL,
                              M          = NULL,
                              init_Z     = NULL,
                              vbk        = NULL,
                              Linf       = NULL,
                              t0         = NULL,
                              len_mid    = NULL,
                              len_border = NULL,
                              a          = NULL,
                              b          = NULL,
                              mat_L50    = NULL,
                              mat_L95    = NULL,
                              cv_L       = NULL,
                              cv_inc     = NULL,
                              std_logR   = NULL,
                              std_logN0  = NULL,
                              alpha      = NULL,
                              beta       = NULL,
                              std_SN     = NULL,
                              q_surv_L50 = NULL,
                              q_surv_L95 = NULL,
                              R_init     = NULL,
                              F_mean     = NULL,
                              F_ar       = NULL,
                              F_sd       = NULL,
                              R_ar       = NULL) {

  # =========================================================================
  # Help: list available presets
  # =========================================================================
  if (!is.null(species) && is.character(species) && length(species) == 1 && species == "?") {
    message("Available species presets for initialize_params():\n",
            "  'flatfish' -- M7 yellowtail flounder (annual, Linf=60, nage=15, M=0.2)\n",
            "  'tuna'     -- M6 tuna (quarterly, Linf=152, nage=20, rec.age=0.25)\n",
            "  'krill'    -- Antarctic krill (Linf=60, nage=5, M=0.3, 5mm bins)\n",
            "Use species=NULL for generic defaults (same as 'flatfish').")
    return(invisible(NULL))
  }

  # =========================================================================
  # Load species preset defaults
  # =========================================================================
  if (!is.null(species)) {
    species <- match.arg(species, c("flatfish", "tuna", "krill"))
  }

  # -----------------------------------------------------------------------
  # FLATFISH (M7) -- Annual, age-based, Linf=60
  # Source: M7.R operating model
  # -----------------------------------------------------------------------
  if (is.null(species) || species == "flatfish") {
    sp <- list(
      model_type  = "age_based",           # age-based dynamics with pla matrix
      growth_step = 1,                     # annual time step
      burn_in     = 80,                    # 80-year burn-in -> 20-year obs window
      nyear      = 100,
      rec.age    = 1,
      first.year = 2020,
      nage       = 15,
      M          = 0.2,
      init_Z     = 0.5,
      vbk        = 0.2,
      Linf       = 60,
      t0         = 1 / 60,
      len_mid    = seq(6, 50, 2),          # 23 bins, 2mm width
      len_border = seq(5, 51, 2),          # 24 borders -> -Inf/Inf applied below
      a          = exp(-12),
      b          = 3,
      mat_L50    = 35,
      mat_L95    = 40,
      cv_L       = 0.2,
      cv_inc     = 0.2,
      std_logR   = 0.3,
      std_logN0  = 0.2,
      alpha      = 400,
      beta       = 10,
      std_SN     = 0.2,
      q_surv_L50 = 15,
      q_surv_L95 = 20,
      sel_type   = "flat",                 # flat selectivity at age: sel = rep(1, nage)
      R_init     = 500,
      F_mean     = 0.3,
      F_ar       = 0.75,
      F_sd       = 0.2,
      R_ar       = 0.1
    )

  # -----------------------------------------------------------------------
  # TUNA (M6) -- Quarterly, length-based dynamics, Linf=152
  # Source: ALSCL_simulation.R (M6 operating model)
  # -----------------------------------------------------------------------
  } else if (species == "tuna") {
    sp <- list(
      model_type  = "length_based",        # 3D length-based dynamics with Gij matrix
      growth_step = 0.25,                  # quarterly time step
      burn_in     = 95,                    # 95-year burn-in -> 5-year obs (20 quarterly steps)
      nyear      = 100,
      rec.age    = 0.25,                   # quarterly recruitment
      first.year = 2020,
      nage       = 20,                     # quarterly ages: 0.25, 0.5, ..., 5.0
      M          = 0.2,                    # quarterly natural mortality
      init_Z     = 0.1,
      vbk        = 0.38,
      Linf       = 152,
      t0         = 1 / 152,
      len_mid    = seq(12.5, 122.5, 5),    # 23 bins, 5mm width
      len_border = seq(10, 125, 5),        # 24 borders
      a          = 2e-5,
      b          = 3,
      mat_L50    = 100,
      mat_L95    = 120,
      cv_L       = 0.2,
      cv_inc     = 0.3,
      std_logR   = 0.5,
      std_logN0  = 0.3,
      alpha      = 15000,
      beta       = 150,
      std_SN     = 0.1,
      q_surv_L50 = 30,
      q_surv_L95 = 50,
      sel_type   = "dome",                 # dome-shaped double logistic at length
      sel_L50_asc  = 30,                   # ascending limb L50
      sel_L95_asc  = 40,                   # ascending limb L95
      sel_L50_desc = 120,                  # descending limb L50 (note: L50 > L95)
      sel_L95_desc = 90,                   # descending limb L95
      R_init     = 10000,
      F_mean     = 0.2,                    # quarterly F
      F_ar       = 0.75,
      F_sd       = 0.2,
      R_ar       = 0.1
    )

  # -----------------------------------------------------------------------
  # KRILL -- Annual, short-lived, Linf=60
  # Source: ACL_krill.R biological parameters + operating model inference
  # -----------------------------------------------------------------------
  } else if (species == "krill") {
    sp <- list(
      model_type  = "age_based",           # age-based dynamics (like flatfish)
      growth_step = 1,                     # annual time step
      burn_in     = 80,                    # 80-year burn-in -> 20-year obs window
      nyear      = 100,
      rec.age    = 1,
      first.year = 2020,
      nage       = 5,
      M          = 0.3,
      init_Z     = 0.5,
      vbk        = 0.3,
      Linf       = 60,
      t0         = 1 / 60,
      len_mid    = c(15, 23, 28, 33, 38, 43, 48, 53, 58),  # 9 bins, ~5mm
      len_border = c(10, 20, 25, 30, 35, 40, 45, 50, 55, 60),  # 10 borders
      a          = exp(-12),
      b          = 3,
      mat_L50    = 40,
      mat_L95    = 45,
      cv_L       = 0.2,
      cv_inc     = 0.2,
      std_logR   = 0.3,
      std_logN0  = 0.2,
      alpha      = 400,
      beta       = 10,
      std_SN     = 0.2,
      q_surv_L50 = 40,                    # survey selects larger krill
      q_surv_L95 = 45,
      sel_type   = "flat",                 # flat selectivity at age
      R_init     = 500,
      F_mean     = 0.3,
      F_ar       = 0.75,
      F_sd       = 0.2,
      R_ar       = 0.1
    )
  }

  # =========================================================================
  # Override with user-supplied values (explicit args beat preset defaults)
  # =========================================================================
  if (!is.null(nyear))      sp$nyear      <- nyear
  if (!is.null(rec.age))    sp$rec.age    <- rec.age
  if (!is.null(first.year)) sp$first.year <- first.year
  if (!is.null(nage))       sp$nage       <- nage
  if (!is.null(M))          sp$M          <- M
  if (!is.null(init_Z))     sp$init_Z     <- init_Z
  if (!is.null(vbk))        sp$vbk        <- vbk
  if (!is.null(Linf))       sp$Linf       <- Linf
  if (!is.null(t0))         sp$t0         <- t0
  if (!is.null(len_mid))    sp$len_mid    <- len_mid
  if (!is.null(len_border)) sp$len_border <- len_border
  if (!is.null(a))          sp$a          <- a
  if (!is.null(b))          sp$b          <- b
  if (!is.null(mat_L50))    sp$mat_L50    <- mat_L50
  if (!is.null(mat_L95))    sp$mat_L95    <- mat_L95
  if (!is.null(cv_L))       sp$cv_L       <- cv_L
  if (!is.null(cv_inc))     sp$cv_inc     <- cv_inc
  if (!is.null(std_logR))   sp$std_logR   <- std_logR
  if (!is.null(std_logN0))  sp$std_logN0  <- std_logN0
  if (!is.null(alpha))      sp$alpha      <- alpha
  if (!is.null(beta))       sp$beta       <- beta
  if (!is.null(std_SN))     sp$std_SN     <- std_SN
  if (!is.null(q_surv_L50)) sp$q_surv_L50 <- q_surv_L50
  if (!is.null(q_surv_L95)) sp$q_surv_L95 <- q_surv_L95
  if (!is.null(R_init))     sp$R_init     <- R_init
  if (!is.null(F_mean))     sp$F_mean     <- F_mean
  if (!is.null(F_ar))       sp$F_ar       <- F_ar
  if (!is.null(F_sd))       sp$F_sd       <- F_sd
  if (!is.null(R_ar))       sp$R_ar       <- R_ar

  # Ensure model_type & growth_step exist (for backward compat / custom usage)
  if (is.null(sp$model_type))  sp$model_type  <- "age_based"
  if (is.null(sp$growth_step)) sp$growth_step <- 1
  if (is.null(sp$burn_in))     sp$burn_in     <- 80
  if (is.null(sp$sel_type))    sp$sel_type    <- "flat"

  # =========================================================================
  # Derived variables
  # =========================================================================
  # Apply -Inf / Inf to first/last border
  sp$len_border[1] <- -Inf
  sp$len_border[length(sp$len_border)] <- Inf

  # Generate age and year vectors
  if (sp$rec.age < 1) {
    # Quarterly or sub-annual: ages are seq(rec.age, max_age, rec.age)
    max_age <- sp$rec.age * sp$nage
    sp$ages <- seq(sp$rec.age, max_age, sp$rec.age)
  } else {
    sp$ages <- sp$rec.age:(sp$rec.age + sp$nage - 1)
  }
  sp$years <- sp$first.year:(sp$first.year + sp$nyear - 1)

  return(sp)
}
