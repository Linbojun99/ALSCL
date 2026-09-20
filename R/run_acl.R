#' Fit the ACL survey catch-at-length model
#' @param data.CatL Data frame or matrix with length labels in its first column,
#'   followed by at least two numerically labelled time columns.
#' @param data.wgt Weight-at-length, with identical labels, columns and dimensions.
#' @param data.mat Maturity-at-length between 0 and 1, with the same layout.
#' @param rec.age Recruitment age in years, greater than t0.
#' @param nage Number of age classes including the plus group (at least two).
#' @param M Nonnegative natural mortality per model time step.
#' @param sel_L50,sel_L95 Lengths at 50 and 95 percent survey catchability.
#' @param parameters Named list of initial fixed-effect values.
#' @param parameters.L,parameters.U Named lower and upper bounds. Bounds are aligned
#'   with free parameters after applying map. Initial free values are clamped to bounds.
#' @param map Named TMB factor mappings. factor(NA) fixes a parameter at its supplied
#'   value; NULL releases a default fixed parameter. Other defaults remain in force.
#' @param len_mid Finite length-bin midpoints; overrides automatic parsing.
#' @param len_border Interior boundaries, one fewer than the number of length bins.
#' @param output Save diagnostic tables and plots below output/ when TRUE.
#' @param train_times Positive integer; exact number of optimization passes per start.
#' @param ncores Positive integer; number of independent starts. Values above one
#'   use socket workers and jitter free parameters only.
#' @param silent Suppress fitting progress messages when TRUE.
#' @param growth_step Time between consecutive ages, in years. Use 0.25 for quarters.
#'   ACL infers rec.age when below one, otherwise one, if this argument is NULL.
#' @param zero_action Exclude zero observations as missing (historical behavior),
#'   or reject them with "error". NA is always missing; negative/infinite values fail.
#' @param control A named list of nlminb control settings.
#' @details M and F are per model time step; vbk is per year. The first and last
#'   length bins absorb the distribution tails. Compilation uses a writable session
#'   cache. Optimizer success, gradient size, Hessian status and boundary hits should
#'   be assessed separately before interpreting estimates or uncertainty.
#' @return A list with report, opt, obj, est_std, vcov, pdHess, gradient,
#'   max_abs_gradient (also final_outer_mgc), convergence_code, bound_hit,
#'   year, age, length-bin metadata, growth_step, elapsed and multi-start diagnostics.
#'   The DLL remains loaded so the returned obj can be evaluated in the same session.
#' @export
run_acl <- function(data.CatL, data.wgt, data.mat, rec.age, nage, M, sel_L50, sel_L95,
                    parameters = NULL, parameters.L = NULL, parameters.U = NULL,
                    map = NULL, len_mid = NULL, len_border = NULL, output = FALSE,
                    train_times = 1, ncores = 1, silent = FALSE, growth_step = NULL,
                    zero_action = c("missing", "error"), control = list()) {
  prepared <- .acl_prepare_data(data.CatL, data.wgt, data.mat, rec.age, nage, M,
    sel_L50, sel_L95, growth_step, len_mid, len_border, zero_action = zero_action)
  .acl_fit("ACL", prepared, parameters, parameters.L, parameters.U, map,
           train_times, ncores, silent, output, control)
}
