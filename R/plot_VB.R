#' Plot the Von Bertalanffy growth function
#'
#' This function takes Linf, k, and t0 as input, creates a sequence of ages, applies the VB function to each age, and then plots the result using ggplot2.
#' @param model_result A list that contains model output. The list should have a "report" component which contains "Linf", "vbk" and "t0" components.
#' @param age_range Numeric vector of length 2, defining the range of ages to consider.
#' @param line_size Numeric. The thickness of the line in the plot. Default is 1.2.
#' @param line_color Character. The color of the line in the plot. Default is "black".
#' @param line_type Character. The type of the line in the plot. Default is "solid".
#' @param se_color Character. The color of the confidence interval ribbon. Default is "blue".
#' @param se_alpha Numeric. The transparency of the confidence interval ribbon. Default is 0.2.
#' @param se_type Character. Type of CI display: "ribbon" (shaded area) or "errorbar" (error bars). Default is "ribbon".
#' @param se Logical. Whether to calculate and plot standard error as confidence intervals. Default is FALSE.
#' @param text_size Numeric. The thickness of the text in the plot. Default is 5.
#' @param text_color Character. The color of the text in the plot. Default is "black".
#' @param title Character or NULL. Custom plot title. If NULL, uses global theme setting. See \code{acl_theme_set()}.
#' @param xlab Character or NULL. Custom x-axis label. If NULL, uses global theme setting.
#' @param ylab Character or NULL. Custom y-axis label. If NULL, uses global theme setting.
#' @param font_family Character or NULL. Custom font family. If NULL, uses global theme setting (default "Arial").
#' @param title_size Numeric or NULL. Plot title size in pt. If NULL, uses global theme (default 14).
#' @param axis_title_size Numeric or NULL. Axis title size in pt. If NULL, uses global theme (default 12).
#' @param axis_text_size Numeric or NULL. Axis tick label size in pt. If NULL, uses global theme (default 10).
#' @param strip_text_size Numeric or NULL. Facet label size in pt. If NULL, uses global theme (default 10).
#' @param legend_text_size Numeric or NULL. Legend text size in pt. If NULL, uses global theme (default 10).
#' @param x_breaks Numeric vector or NULL. Custom x-axis breaks (e.g. \code{seq(1, 20, by = 2)}). NULL = auto.
#' @param base_theme Character or NULL. Base ggplot2 theme name (e.g. "theme_bw"). NULL = global setting.
#' @param title_hjust Numeric or NULL. Title horizontal alignment: 0 = left, 0.5 = center, 1 = right. NULL = global setting.
#' @return A ggplot object representing the plot.
#' @export
#' @examples
#' \dontrun{
#' plot_VB(model_result, age_range = c(1, 25))
#' }
plot_VB <- function(model_result, age_range = c(1, 25), line_size = 1.2, line_color = "red", line_type = "solid", se = FALSE, se_color = "red", se_alpha = 0.2, se_type = c("ribbon", "errorbar"),text_color="black",text_size=5, title = NULL, xlab = NULL, ylab = NULL, font_family = NULL, title_size = NULL, axis_title_size = NULL, axis_text_size = NULL, strip_text_size = NULL, legend_text_size = NULL, x_breaks = NULL, base_theme = NULL, title_hjust = NULL){
  # Define the VB function
  VB_func <- function(Linf, k, t0, age) {
    Lt = Linf * (1 - exp(-k * (age - t0)))
    return(Lt)
  }

  Linf = model_result[["report"]][["Linf"]]
  k = model_result[["report"]][["vbk"]]
  t0 = model_result[["report"]][["t0"]]

  # Create a data frame with a fine sequence of ages for smooth curve
  data <- data.frame(age = seq(age_range[1], age_range[2], by = 0.1))
  data <- data %>%
    mutate(length = VB_func(Linf, k, t0, age))





  if (!se)
  {
    # Plot the VB function using ggplot2
    p <- ggplot2::ggplot(data, aes(x = age, y = length)) +
      ggplot2::geom_line(linewidth = line_size, color = line_color, linetype = line_type) +
      ggplot2::labs(x = if (!is.null(xlab)) xlab else .acl_lab("x", "age"), y = "Length", title = if (!is.null(title)) title else .acl_title("VB")) +
      ggplot2::annotate("text", x = -Inf, y = Inf,
                        label = paste("Linf =", round(Linf, 2), "\nk =", round(k, 2)),
                        hjust = -0.1, vjust = 1.5, size = text_size, color = text_color) +
      .acl_scale_x(x_breaks, n_breaks = 10) +
      .acl_base_theme(font_family, title_size, axis_title_size, axis_text_size, strip_text_size, legend_text_size, base_theme = base_theme, title_hjust = title_hjust)

  }
  else{

    # --- Extract SE for log_vbk ---
    ss_logk <- model_result[["est_std"]][grep("^log_vbk", rownames(model_result[["est_std"]])), ]
    ss_logk <- as.data.frame(ss_logk)
    se_val  <- ss_logk["Std. Error", ]

    # Check if SE is NaN (common when parameter hits a boundary)
    se_available <- is.finite(se_val) && se_val > 0

    if (!se_available) {
      # Detect if the parameter hit a boundary
      par_lu <- model_result[["par_low_up"]]
      if (!is.null(par_lu) && "log_vbk" %in% rownames(par_lu)) {
        est <- par_lu["log_vbk", 1]
        lo  <- par_lu["log_vbk", 2]
        hi  <- par_lu["log_vbk", 3]
        if (isTRUE(abs(est - lo) < 1e-6)) {
          message("plot_VB: Cannot compute CI -- log_vbk hit its LOWER bound (vbk = ",
                  round(exp(lo), 4), "). The Hessian is singular at the boundary, ",
                  "so TMB returns NaN for standard errors.\n",
                  "  Fix: widen the lower bound for log_vbk in create_parameters(), ",
                  "or use a different starting value.")
        } else if (isTRUE(abs(est - hi) < 1e-6)) {
          message("plot_VB: Cannot compute CI -- log_vbk hit its UPPER bound (vbk = ",
                  round(exp(hi), 4), "). The Hessian is singular at the boundary.")
        } else {
          message("plot_VB: Cannot compute CI -- Std. Error for log_vbk is NaN. ",
                  "This usually means a growth parameter hit a boundary or the ",
                  "Hessian is not positive definite.")
        }
      } else {
        message("plot_VB: Cannot compute CI -- Std. Error for log_vbk is NaN. ",
                "Falling back to line-only plot.")
      }

      # Fall back: draw line only, with annotation noting CI unavailable
      data <- data.frame(age = seq(age_range[1], age_range[2], by = 0.1))
      data$length <- VB_func(Linf, k, t0, data$age)

      p <- ggplot2::ggplot(data, ggplot2::aes(x = age, y = length)) +
        ggplot2::geom_line(linewidth = line_size, color = line_color, linetype = line_type) +
        ggplot2::labs(
          x = if (!is.null(xlab)) xlab else .acl_lab("x", "age"),
          y = "Length",
          title    = if (!is.null(title)) title else .acl_title("VB"),
          subtitle = "CI unavailable (parameter hit boundary; SE = NaN)"
        ) +
        ggplot2::annotate("text", x = -Inf, y = Inf,
                          label = paste("Linf =", round(Linf, 2), "\nk =", round(k, 2)),
                          hjust = -0.1, vjust = 1.5, size = text_size, color = text_color) +
        .acl_scale_x(x_breaks, n_breaks = 10) +
        .acl_base_theme(font_family, title_size, axis_title_size, axis_text_size,
                        strip_text_size, legend_text_size, base_theme = base_theme,
                        title_hjust = title_hjust)

    } else {
      # SE is valid: compute CI via delta method on k
      ss_k <- data.frame(
        estimate = exp(ss_logk["Estimate", ]),
        lower    = exp(ss_logk["Estimate", ] - 1.96 * se_val),
        upper    = exp(ss_logk["Estimate", ] + 1.96 * se_val)
      )

      data <- data.frame(age = seq(age_range[1], age_range[2], by = 0.1))
      data <- data %>%
        dplyr::mutate(
          length       = VB_func(Linf, k, t0, age),
          length_lower = VB_func(Linf, ss_k$lower, t0, age),
          length_upper = VB_func(Linf, ss_k$upper, t0, age)
        )

      p <- ggplot2::ggplot(data, ggplot2::aes(x = age, y = length)) +
        { if (se_type[1] == "ribbon")
          ggplot2::geom_ribbon(ggplot2::aes(ymin = length_lower, ymax = length_upper),
                               alpha = se_alpha, fill = se_color)
          else
            ggplot2::geom_errorbar(ggplot2::aes(ymin = length_lower, ymax = length_upper),
                                   color = se_color, width = 0.3) } +
        ggplot2::geom_line(linewidth = line_size, color = line_color, linetype = line_type) +
        ggplot2::labs(
          x = if (!is.null(xlab)) xlab else .acl_lab("x", "age"),
          y = "Length",
          title = if (!is.null(title)) title else .acl_title("VB")
        ) +
        ggplot2::annotate("text", x = -Inf, y = Inf,
                          label = paste("Linf =", round(Linf, 2), "\nk =", round(k, 2)),
                          hjust = -0.1, vjust = 1.5, size = text_size, color = text_color) +
        .acl_scale_x(x_breaks, n_breaks = 10) +
        .acl_base_theme(font_family, title_size, axis_title_size, axis_text_size,
                        strip_text_size, legend_text_size, base_theme = base_theme,
                        title_hjust = title_hjust)
    }
  }

  return(p)
}
