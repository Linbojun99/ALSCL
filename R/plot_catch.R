utils::globalVariables(c("AgeGroup", "LengthGroup", "CNL", "CNA"))
#' Plot Catch Number (CN), Catch by Age (CNA), or Catch by Length (CNL) Over Years
#'
#' Visualizes catch data from ACL or ALSCL model outputs.
#' Supports total catch ("CN"), catch at age ("CNA"), and catch at length ("CNL").
#'
#' @param model_result A list from \code{run_acl} or \code{run_alscl}.
#' @param line_size Numeric. Line thickness. Default is 1.2.
#' @param line_color Character. Line color. Default is "red".
#' @param line_type Character. Line type. Default is "solid".
#' @param se Logical. Whether to plot confidence intervals. Default is FALSE.
#' @param se_color Character. CI ribbon color. Default is "red".
#' @param se_alpha Numeric. CI ribbon transparency. Default is 0.2.
#' @param facet_ncol Numeric. Columns in facet wrap. Default is NULL.
#' @param facet_scales Character. Scales for facet wrap. Default is "free".
#' @param type Character. "CN" (total), "CNA" (at age), or "CNL" (at length). Default is "CN".
#' @param return_data Logical. Whether to return processed data. Default is FALSE.
#'
#' @return A ggplot object or a list with plot and data.
#' @export
plot_catch <- function(model_result, line_size = 1.2, line_color = "red", line_type = "solid",
                       se = FALSE, se_color = "red", se_alpha = 0.2,
                       facet_ncol = NULL, facet_scales = "free",
                       type = c("CN", "CNA", "CNL"), return_data = FALSE) {

  type <- match.arg(type)
  Year <- model_result[["year"]]
  len_label <- model_result[["len_label"]]
  current_theme <- tryCatch(get("acl_get_theme", envir = asNamespace("ACL"))(), error = function(e) ggplot2::theme_minimal())

  # ==========================
  # Type = "CN": Total catch
  # ==========================
  if (type == "CN") {
    CN <- model_result[["report"]][["CN"]]
    if (!is.data.frame(CN)) CN <- as.data.frame(CN)
    CN$Year <- Year

    if (!se) {
      p <- ggplot2::ggplot(CN, ggplot2::aes(x = Year, y = CN)) +
        ggplot2::geom_line(linewidth = line_size, color = line_color, linetype = line_type) +
        ggplot2::labs(x = "Year", y = "Total Catch Numbers", title = "Total Catch Over Years") +
        current_theme
      data_out <- CN
    } else {
      # Use exact match ^CN$ to avoid matching CNA, CNL
      ss_CN <- model_result[["est_std"]][grep("^CN$", rownames(model_result[["est_std"]])), ]
      ci <- data.frame(
        Year     = Year,
        estimate = ss_CN[, "Estimate"],
        lower    = ss_CN[, "Estimate"] - 1.96 * ss_CN[, "Std. Error"],
        upper    = ss_CN[, "Estimate"] + 1.96 * ss_CN[, "Std. Error"]
      )
      p <- ggplot2::ggplot(ci, ggplot2::aes(x = Year, y = estimate)) +
        ggplot2::geom_line(linewidth = line_size, color = line_color, linetype = line_type) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = se_color, alpha = se_alpha) +
        ggplot2::labs(y = "Total Catch Numbers", x = "Year", title = "Total Catch Over Years") +
        current_theme
      data_out <- ci
    }
  }

  # ==========================
  # Type = "CNA": Catch at age
  # ==========================
  if (type == "CNA") {
    CNA <- model_result[["report"]][["CNA"]]

    AgeLabels <- factor(paste("Age bin", seq_len(nrow(CNA))),
                        levels = paste("Age bin", seq_len(nrow(CNA))))

    CNA_long <- reshape2::melt(CNA)
    colnames(CNA_long) <- c("AgeGroup", "Year", "CNA")
    CNA_long$Year <- Year[as.numeric(CNA_long$Year)]
    CNA_long$AgeGroup <- AgeLabels[as.numeric(CNA_long$AgeGroup)]

    if (!se) {
      p <- ggplot2::ggplot(CNA_long, ggplot2::aes(x = Year, y = CNA)) +
        ggplot2::geom_line(linewidth = line_size, color = line_color, linetype = line_type) +
        ggplot2::facet_wrap(~AgeGroup, ncol = facet_ncol, scales = facet_scales) +
        ggplot2::labs(x = "Year", y = "Catch Numbers", title = "Catch at Age Over Years") +
        current_theme
      data_out <- CNA_long
    } else {
      ss_cna <- model_result[["est_std"]][grep("^CNA$", rownames(model_result[["est_std"]])), ]
      ss_cna <- as.data.frame(ss_cna)
      ci <- data.frame(
        AgeGroup = rep(levels(AgeLabels), times = ncol(CNA)),
        Year     = rep(Year, each = nrow(CNA)),
        estimate = ss_cna[, "Estimate"],
        lower    = ss_cna[, "Estimate"] - 1.96 * ss_cna[, "Std. Error"],
        upper    = ss_cna[, "Estimate"] + 1.96 * ss_cna[, "Std. Error"]
      )
      p <- ggplot2::ggplot(ci, ggplot2::aes(x = Year, y = estimate)) +
        ggplot2::geom_line(linewidth = line_size, color = line_color, linetype = line_type) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = se_color, alpha = se_alpha) +
        ggplot2::facet_wrap(~AgeGroup, ncol = facet_ncol, scales = facet_scales) +
        ggplot2::labs(x = "Year", y = "Catch Numbers", title = "Catch at Age Over Years") +
        current_theme
      data_out <- ci
    }
  }

  # ==========================
  # Type = "CNL": Catch at length
  # ==========================
  if (type == "CNL") {
    CNL <- model_result[["report"]][["CNL"]]
    if (!is.matrix(CNL)) CNL <- as.matrix(CNL)

    LengthLabels <- .acl_fix_len_labels(len_label)

    CNL_long <- reshape2::melt(CNL)
    colnames(CNL_long) <- c("LengthGroup", "Year", "CNL")
    CNL_long$Year <- Year[as.numeric(CNL_long$Year)]
    CNL_long$LengthGroup <- factor(LengthLabels[as.numeric(CNL_long$LengthGroup)],
                                   levels = LengthLabels)

    if (!se) {
      p <- ggplot2::ggplot(CNL_long, ggplot2::aes(x = Year, y = CNL)) +
        ggplot2::geom_line(linewidth = line_size, color = line_color, linetype = line_type) +
        ggplot2::facet_wrap(~LengthGroup, ncol = facet_ncol, scales = facet_scales) +
        ggplot2::labs(x = "Year", y = "Catch Numbers", title = "Catch at Length Over Years") +
        current_theme
      data_out <- CNL_long
    } else {
      ss_cnl <- model_result[["est_std"]][grep("^CNL$", rownames(model_result[["est_std"]])), ]
      ss_cnl <- as.data.frame(ss_cnl)
      ci <- data.frame(
        LengthGroup = factor(rep(LengthLabels, times = ncol(CNL)), levels = LengthLabels),
        Year        = rep(Year, each = nrow(CNL)),
        estimate    = ss_cnl[, "Estimate"],
        lower       = ss_cnl[, "Estimate"] - 1.96 * ss_cnl[, "Std. Error"],
        upper       = ss_cnl[, "Estimate"] + 1.96 * ss_cnl[, "Std. Error"]
      )
      p <- ggplot2::ggplot(ci, ggplot2::aes(x = Year, y = estimate)) +
        ggplot2::geom_line(linewidth = line_size, color = line_color, linetype = line_type) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = se_color, alpha = se_alpha) +
        ggplot2::facet_wrap(~LengthGroup, ncol = facet_ncol, scales = facet_scales) +
        ggplot2::labs(x = "Year", y = "Catch Numbers", title = "Catch at Length Over Years") +
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
