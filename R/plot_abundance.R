utils::globalVariables(c("AgeGroup", "LengthGroup", "Count", "number.age"))
#' Plot Abundance over the Years in ACL or ALSCL model
#'
#' This function visualizes fish abundance data from ACL or ALSCL assessments.
#' It handles three types: total number ("N"), number at age ("NA"), and number at length ("NL").
#'
#' @param model_result A list from \code{run_acl} or \code{run_alscl}.
#' @param line_size Numeric. Line thickness. Default is 1.2.
#' @param line_color Character. Line color. Default is "red".
#' @param line_type Character. Line type. Default is "solid".
#' @param se Logical. Whether to plot confidence intervals. Default is FALSE.
#' @param se_color Character. CI ribbon color. Default is "red".
#' @param se_alpha Numeric. CI ribbon transparency. Default is 0.2.
#' @param type Character. "N" (total), "NA" (at age), or "NL" (at length). Default is "N".
#' @param facet_ncol Number of columns in facet wrap. Default is NULL.
#' @param facet_scales Scales for facet wrap. Default is "free".
#' @param return_data Logical. Whether to return processed data. Default is FALSE.
#' @return A ggplot object or a list with plot and data.
#' @export
plot_abundance <- function(model_result, line_size = 1.2, line_color = "red", line_type = "solid",
                           se = FALSE, se_color = "red", se_alpha = 0.2,
                           type = c("N", "NA", "NL"), facet_ncol = NULL, facet_scales = "free",
                           return_data = FALSE) {

  type <- match.arg(type)
  len_label <- model_result[["len_label"]]
  Year <- model_result[["year"]]
  current_theme <- tryCatch(get("acl_get_theme", envir = asNamespace("ACL"))(), error = function(e) ggplot2::theme_minimal())

  # ==========================
  # Type = "N": Total abundance
  # ==========================
  if (type == "N") {
    number <- model_result[["report"]][["N"]]
    if (!is.data.frame(number)) number <- as.data.frame(number)
    number$Year <- Year

    if (!se) {
      p <- ggplot2::ggplot(number, ggplot2::aes(x = Year, y = number)) +
        ggplot2::geom_line(linewidth = line_size, color = line_color, linetype = line_type) +
        ggplot2::labs(x = "Year", y = "Relative abundance", title = "Total Abundance Over Years") +
        current_theme
      data_out <- number
    } else {
      # Use exact match ^N$ to avoid matching NA, NL, NLA
      ss_n <- model_result[["est_std"]][grepl("^N$", rownames(model_result[["est_std"]])), ]
      ci <- data.frame(
        Year     = Year,
        estimate = ss_n[, "Estimate"],
        lower    = ss_n[, "Estimate"] - 1.96 * ss_n[, "Std. Error"],
        upper    = ss_n[, "Estimate"] + 1.96 * ss_n[, "Std. Error"]
      )
      p <- ggplot2::ggplot(ci, ggplot2::aes(x = Year, y = estimate)) +
        ggplot2::geom_line(linewidth = line_size, color = line_color, linetype = line_type) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = se_color, alpha = se_alpha) +
        ggplot2::labs(y = "Relative abundance", x = "Year", title = "Total Abundance Over Years") +
        current_theme
      data_out <- ci
    }
  }

  # ==========================
  # Type = "NA": Abundance at age
  # ==========================
  if (type == "NA") {
    na <- model_result[["report"]][["NA"]]

    AgeGroup <- factor(paste("Age bin", seq_len(nrow(na))),
                       levels = paste("Age bin", seq_len(nrow(na))))

    na_long <- reshape2::melt(na)
    colnames(na_long) <- c("AgeGroup", "Year", "number.age")
    na_long$Year <- Year[as.numeric(na_long$Year)]
    na_long$AgeGroup <- AgeGroup[as.numeric(na_long$AgeGroup)]

    if (!se) {
      p <- ggplot2::ggplot(na_long, ggplot2::aes(x = Year, y = number.age)) +
        ggplot2::geom_line(linewidth = line_size, color = line_color, linetype = line_type) +
        ggplot2::facet_wrap(~AgeGroup, ncol = facet_ncol, scales = facet_scales) +
        ggplot2::labs(x = "Year", y = "Relative abundance", title = "Abundance at Age Over Years") +
        current_theme
      data_out <- na_long
    } else {
      # Use exact match ^NA$ to avoid matching NLA
      ss_na <- model_result[["est_std"]][grep("^NA$", rownames(model_result[["est_std"]])), ]
      ss_na <- as.data.frame(ss_na)

      ci <- data.frame(
        AgeGroup = rep(levels(AgeGroup), times = ncol(na)),
        Year     = rep(Year, each = nrow(na)),
        estimate = ss_na[, "Estimate"],
        lower    = ss_na[, "Estimate"] - 1.96 * ss_na[, "Std. Error"],
        upper    = ss_na[, "Estimate"] + 1.96 * ss_na[, "Std. Error"]
      )

      p <- ggplot2::ggplot(ci, ggplot2::aes(x = Year, y = estimate)) +
        ggplot2::geom_line(linewidth = line_size, color = line_color, linetype = line_type) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = se_color, alpha = se_alpha) +
        ggplot2::facet_wrap(~AgeGroup, ncol = facet_ncol, scales = facet_scales) +
        ggplot2::labs(x = "Year", y = "Relative abundance", title = "Abundance at Age Over Years") +
        current_theme
      data_out <- ci
    }
  }

  # ==========================
  # Type = "NL": Abundance at length
  # ==========================
  if (type == "NL") {
    NL <- model_result[["report"]][["NL"]]
    if (!is.matrix(NL)) NL <- as.matrix(NL)

    LengthLabels <- .acl_fix_len_labels(len_label)

    NL_long <- reshape2::melt(NL)
    colnames(NL_long) <- c("LengthGroup", "Year", "Count")
    NL_long$Year <- Year[as.numeric(NL_long$Year)]
    NL_long$LengthGroup <- factor(LengthLabels[as.numeric(NL_long$LengthGroup)],
                                  levels = LengthLabels)

    if (!se) {
      p <- ggplot2::ggplot(NL_long, ggplot2::aes(x = Year, y = Count)) +
        ggplot2::geom_line(linewidth = line_size, color = line_color, linetype = line_type) +
        ggplot2::facet_wrap(~LengthGroup, ncol = facet_ncol, scales = facet_scales) +
        ggplot2::labs(x = "Year", y = "Relative abundance", title = "Abundance at Length Over Years") +
        current_theme
      data_out <- NL_long
    } else {
      # Use exact match ^NL$ to avoid matching NLA
      ss_NL <- model_result[["est_std"]][grep("^NL$", rownames(model_result[["est_std"]])), ]
      ss_NL <- as.data.frame(ss_NL)

      ci <- data.frame(
        LengthGroup = factor(rep(LengthLabels, times = ncol(NL)), levels = LengthLabels),
        Year        = rep(Year, each = nrow(NL)),
        estimate    = ss_NL[, "Estimate"],
        lower       = ss_NL[, "Estimate"] - 1.96 * ss_NL[, "Std. Error"],
        upper       = ss_NL[, "Estimate"] + 1.96 * ss_NL[, "Std. Error"]
      )

      p <- ggplot2::ggplot(ci, ggplot2::aes(x = Year, y = estimate)) +
        ggplot2::geom_line(linewidth = line_size, color = line_color, linetype = line_type) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = se_color, alpha = se_alpha) +
        ggplot2::facet_wrap(~LengthGroup, ncol = facet_ncol, scales = facet_scales) +
        ggplot2::labs(x = "Year", y = "Relative abundance", title = "Abundance at Length Over Years") +
        current_theme
      data_out <- ci
    }
  }

  if (return_data) {
    return(list(plot = p, data = data_out))
  } else {
    return(p)
  }
}
