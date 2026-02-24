#' Plot Deviation of Recruitment (R) and Fishing Mortality (F) Over Years
#'
#' This function calculates the deviations of R and F in the ACL or ALSCL model results.
#' For ACL, F deviations are at the age level (A x Y).
#' For ALSCL, F deviations are at the length level (L x Y).
#'
#' @param model_result A list from \code{run_acl} or \code{run_alscl}.
#' @param se Logical, whether to draw the standard error bar.
#' @param line_size Numeric, the line size.
#' @param line_color Character, the line color.
#' @param line_type Character, the line type.
#' @param point_color Character. The color of the point. Default is "white".
#' @param point_size Numeric. The size of the point. Default is 3.
#' @param point_shape Numeric. The shape of the point. Default is 21.
#' @param se_color Character, the color of the standard error.
#' @param se_width Numeric, the width of the standard error.
#' @param log Logical, whether to keep log scale (TRUE) or apply exp transform (FALSE).
#' @param facet_ncol Number of columns in facet wrap. Default is NULL.
#' @param facet_scales Scales for facet wrap. Default is "free".
#' @param type Character. "R" for recruitment deviation, "F" for fishing mortality deviation.
#'
#' @return A ggplot2 object.
#' @export
plot_deviance <- function(model_result, se = TRUE, point_size = 3, point_color = "white",
                          point_shape = 21, line_size = 1, line_color = "black",
                          line_type = "solid", se_color = "black", se_width = 0.5,
                          facet_ncol = NULL, facet_scales = "free", log = TRUE,
                          type = c("R", "F")) {

  type <- match.arg(type)

  # Detect model type
  is_alscl <- is.null(model_result[["report"]][["F"]]) && !is.null(model_result[["report"]][["FL"]])

  # Get theme
  current_theme <- tryCatch(get("acl_get_theme", envir = asNamespace("ACL"))(), error = function(e) ggplot2::theme_minimal())

  # ==============================================================
  # Type = "R": Recruitment deviations (same for ACL and ALSCL)
  # ==============================================================
  if (type == "R") {
    dev_log_R <- model_result[["est_std"]][grep("^dev_log_R", rownames(model_result[["est_std"]])), ]
    if (!is.data.frame(dev_log_R)) dev_log_R <- as.data.frame(dev_log_R)

    ci <- data.frame(
      Year     = model_result[["year"]],
      estimate = dev_log_R[, "Estimate"],
      lower    = dev_log_R[, "Estimate"] - 1.96 * dev_log_R[, "Std. Error"],
      upper    = dev_log_R[, "Estimate"] + 1.96 * dev_log_R[, "Std. Error"]
    )

    if (!log) {
      ci[, c("estimate", "lower", "upper")] <- exp(ci[, c("estimate", "lower", "upper")])
    }

    plot <- ggplot2::ggplot(ci, ggplot2::aes(x = Year, y = estimate)) +
      ggplot2::geom_line(linewidth = line_size, color = line_color, linetype = line_type) +
      ggplot2::geom_point(size = point_size, fill = point_color, shape = point_shape) +
      ggplot2::labs(y = "Recruitment Deviance", x = "Year", title = "Recruitment Deviance Over Years") +
      current_theme

    if (se) {
      plot <- plot +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), color = se_color, width = se_width) +
        ggplot2::geom_point(size = point_size, fill = point_color, shape = point_shape)
    }
  }

  # ==============================================================
  # Type = "F": Fishing mortality deviations
  # ==============================================================
  if (type == "F") {
    dev_log_F <- model_result[["est_std"]][grep("^dev_log_F", rownames(model_result[["est_std"]])), ]
    dev_log_F <- as.data.frame(dev_log_F)

    if (is_alscl) {
      # ALSCL: dev_log_F is (L x Y), faceted by length bin
      FL <- model_result[["report"]][["FL"]]
      len_label <- model_result[["len_label"]]
      if (is.null(len_label)) len_label <- paste0("L", seq_len(nrow(FL)))

      group_labels <- .acl_fix_len_labels(len_label, prefix = "Length")
      n_groups <- nrow(FL)
      n_years  <- ncol(FL)

      dev_log_F$Group <- factor(rep(group_labels, times = n_years), levels = group_labels)
      dev_log_F$Year  <- rep(model_result[["year"]], each = n_groups)

      facet_title <- "F Deviance at Length (ALSCL)"
    } else {
      # ACL: dev_log_F is (A x Y), faceted by age group
      Fmat <- model_result[["report"]][["F"]]
      n_groups <- nrow(Fmat)
      n_years  <- ncol(Fmat)

      dev_log_F$Group <- rep(paste0("Age group ", seq_len(n_groups)), times = n_years)
      dev_log_F$Year  <- rep(model_result[["year"]], each = n_groups)

      facet_title <- "F Deviance Over Years"
    }

    ci <- data.frame(
      Group    = dev_log_F$Group,
      Year     = dev_log_F$Year,
      estimate = dev_log_F[, "Estimate"],
      lower    = dev_log_F[, "Estimate"] - 1.96 * dev_log_F[, "Std. Error"],
      upper    = dev_log_F[, "Estimate"] + 1.96 * dev_log_F[, "Std. Error"]
    )

    if (!log) {
      ci[, c("estimate", "lower", "upper")] <- exp(ci[, c("estimate", "lower", "upper")])
    }

    plot <- ggplot2::ggplot(ci, ggplot2::aes(x = Year, y = estimate)) +
      ggplot2::geom_line(linewidth = line_size, color = line_color, linetype = line_type) +
      ggplot2::geom_point(size = point_size, fill = point_color, shape = point_shape) +
      ggplot2::labs(y = "F Deviance", x = "Year", title = facet_title) +
      current_theme +
      ggplot2::facet_wrap(~Group, ncol = facet_ncol, scales = facet_scales)

    if (se) {
      plot <- plot +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), color = se_color, width = se_width) +
        ggplot2::geom_point(size = point_size, fill = point_color, shape = point_shape)
    }
  }

  return(plot)
}
