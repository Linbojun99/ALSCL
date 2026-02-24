#' ACL Plot Theme Configuration
#'
#' Global settings for all ACL plot functions, including model comparison plots.
#' Change titles, axis labels, font family, colors, line sizes, and
#' comparison plot styles in one place.
#'
#' @section Usage:
#' \preformatted{
#' # View current settings
#' acl_theme()
#'
#' # Change just one thing
#' acl_theme_set(font_family = "Times New Roman")
#' acl_theme_set(base_theme = "theme_classic")
#' acl_theme_set(compare_linewidth = 2)
#'
#' # ACL vs ALSCL comparison style
#' acl_theme_set(
#'   compare_colors = c("steelblue", "tomato"),
#'   compare_linetypes = c("solid", "dashed")
#' )
#'
#' # Same model, different parameters
#' acl_theme_set(
#'   compare_colors = c("M=0" = "blue", "M=0.3" = "red")
#' )
#' plot_compare_ts(result1, result2,
#'   model1_name = "M=0", model2_name = "M=0.3")
#'
#' # Or just override per-call (no global setting needed)
#' plot_compare_ts(result1, result2,
#'   model1_name = "Run A", model2_name = "Run B",
#'   colors = c("darkgreen", "purple"))
#'
#' # All plots return ggplot objects, so you can layer:
#' plot_compare_ts(result1, result2) +
#'   ggplot2::theme(legend.position = "top") +
#'   ggplot2::labs(title = "My Title")
#'
#' # Reset to defaults
#' acl_theme_reset()
#' }
#'
#' @name acl_theme
NULL

# ---------------------------------------------------------------------------
# Default configuration (internal)
# ---------------------------------------------------------------------------
.acl_defaults <- list(

  # --- Base ggplot2 theme --------------------------------------------------
  base_theme = "theme_bw",

  # --- Font ----------------------------------------------------------------
  font_family = "sans",

  # --- Text sizes (in pt) -------------------------------------------------
  title_size      = 14,
  title_hjust     = 0.5,
  axis_title_size = 12,
  axis_text_size  = 10,
  strip_text_size = 10,
  legend_text_size = 10,

  # --- Axis settings -------------------------------------------------------
  x_breaks = NULL,
  x_expand = c(0.01, 0.01),

  # --- Default colors and line sizes (single-model plots) ------------------
  line_color  = "#D32F2F",
  line_size   = 1.5,
  se_color    = "#D32F2F",
  se_alpha    = 0.2,

  # --- Comparison plot settings (plot_compare_* functions) ------------------
  # Colors/linetypes: positional (1st = model1, 2nd = model2).
  # Names are optional; if provided they override positional matching.
  # Examples:
  #   c("#2166AC", "#B2182B")                       positional
  #   c(ACL = "#2166AC", ALSCL = "#B2182B")         named
  #   c("M=0" = "steelblue", "M=0.5" = "tomato")   custom names
  compare_colors     = c("#2166AC", "#B2182B"),
  compare_linetypes  = c("solid", "dashed"),
  compare_linewidth  = 1.2,
  compare_point_size = 2,
  compare_se_alpha   = 0.15,
  compare_legend_pos = "bottom",
  compare_facet_ncol = 2,
  compare_facet_scales = "free_y",

  # --- Titles --------------------------------------------------------------
  titles = list(
    # plot_abundance
    N       = "Estimated Total Abundance (N)",
    N_se    = "Estimated Total Abundance (N) with 95% CI",
    NA_     = "Estimated Abundance-at-Age",
    NA_se   = "Estimated Abundance-at-Age with 95% CI",
    NL      = "Estimated Abundance-at-Length",
    NL_se   = "Estimated Abundance-at-Length with 95% CI",

    # plot_biomass
    B       = "Estimated Total Biomass (B)",
    B_se    = "Estimated Total Biomass (B) with 95% CI",
    BL      = "Estimated Biomass-at-Length",

    # plot_SSB
    SSB     = "Estimated Spawning Stock Biomass (SSB)",
    SSB_se  = "Estimated Spawning Stock Biomass (SSB) with 95% CI",
    SBL     = "Estimated SSB-at-Length",

    # plot_recruitment
    Rec     = "Estimated Recruitment (R)",
    Rec_se  = "Estimated Recruitment (R) with 95% CI",

    # plot_SSB_Rec
    SSB_Rec = "Spawning Stock Biomass vs. Recruitment",

    # plot_catch
    CN      = "Estimated Catch Abundance (C)",
    CN_se   = "Estimated Catch Abundance (C) with 95% CI",
    CNA     = "Estimated Catch-at-Age",

    # plot_fishing_mortality
    F_year     = "Estimated Fishing Mortality by Year (F)",
    F_year_se  = "Estimated Fishing Mortality by Year (F) with 95% CI",
    F_age      = "Estimated Fishing Mortality-at-Age (F)",
    F_age_se   = "Estimated Fishing Mortality-at-Age (F) with 95% CI",

    # plot_deviance
    dev_R   = "Recruitment Deviance (dev log R)",
    dev_F   = "Fishing Mortality Deviance (dev log F)",

    # plot_VB
    VB      = "Von Bertalanffy Growth Curve",

    # plot_CatL
    CatL_obs_length = "Observed Catch-at-Length Over Years",
    CatL_est_length = "Estimated Catch-at-Length Over Years",
    CatL_year       = "Estimated(Line) and Observed(Point) Catch-at-Length Over Years",
    CatL_year_dist  = "Estimated(Red) and Observed(Blue) Catch-at-Length Distribution Yearly",

    # plot_residuals
    resid_length = "Residuals by Length Bins",
    resid_year   = "Residuals by Year",

    # plot_ridges
    ridges   = "Observed vs. Estimated Catch-at-Length",

    # retro_acl
    retro    = "Retrospective Analysis",

    # plot_pla
    pla      = "Age-Length Transition Probability Matrix",

    # Comparison plot titles
    compare_ts      = "{m1} vs {m2}",
    compare_F       = "Fishing Mortality Comparison",
    compare_resid   = "Residual Diagnostics: {m1} vs {m2}",
    compare_growth  = "Growth Comparison",
    compare_CatL    = "Observed vs Predicted Catch-at-Length",
    compare_metrics = "Model Fit Comparison: {m1} vs {m2}",
    compare_sel     = "Survey Selectivity at Length",
    compare_annualF = "Mean Annual Fishing Mortality"
  ),

  # --- Axis labels ---------------------------------------------------------
  xlab = list(
    year   = "Year",
    age    = "Age",
    length = "Length",
    ssb    = "SSB"
  ),

  ylab = list(
    abundance = "Relative abundance",
    biomass   = "Biomass",
    ssb       = "SSB",
    rec       = "Recruitment",
    catch     = "Catch abundance",
    F         = "Fishing mortality",
    deviance  = "Deviance",
    residual  = "Residual",
    growth    = "Length"
  )
)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

#' Get Current ACL Plot Theme
#'
#' Returns the current global theme settings used by all ACL plot functions.
#'
#' @param what Optional character. Retrieve a specific element, e.g.
#'   \code{"titles"}, \code{"font_family"}, \code{"compare_colors"}.
#' @return A list of current settings, or a single element if \code{what}
#'   is specified.
#' @export
#' @examples
#' acl_theme()
#' acl_theme("font_family")
#' acl_theme("compare_colors")
acl_theme <- function(what = NULL) {
  current <- getOption("acl.theme", .acl_defaults)
  if (is.null(what)) return(current)
  current[[what]]
}


#' Set ACL Plot Theme Options
#'
#' Modify one or more global theme settings. Partial updates are supported:
#' only the fields you specify will be changed; all others keep their
#' current values.
#'
#' @param base_theme Character. Base ggplot2 theme name.
#' @param font_family Character. Font family for all text.
#' @param title_size Numeric. Plot title size in pt.
#' @param title_hjust Numeric. Title alignment (0=left, 0.5=center, 1=right).
#' @param axis_title_size Numeric. Axis title size.
#' @param axis_text_size Numeric. Axis tick label size.
#' @param strip_text_size Numeric. Facet label size.
#' @param legend_text_size Numeric. Legend text size.
#' @param x_breaks Numeric vector or NULL. Custom x-axis breaks.
#' @param x_expand Numeric vector of length 2. X-axis expansion.
#' @param line_color Character. Default line color for single-model plots.
#' @param line_size Numeric. Default line width for single-model plots.
#' @param se_color Character. Default CI color.
#' @param se_alpha Numeric. Default CI transparency.
#' @param compare_colors Character vector of 2 colors. Can be unnamed (positional)
#'   or named. The first color is always model1, second is model2.
#'   Examples: \code{c("blue", "red")}, \code{c(ACL = "blue", ALSCL = "red")},
#'   \code{c("M=0" = "steelblue", "M=0.5" = "tomato")}.
#' @param compare_linetypes Character vector of 2 linetypes. Same rules as colors.
#' @param compare_linewidth Numeric. Line width for comparison plots.
#' @param compare_point_size Numeric. Point size for comparison plots.
#' @param compare_se_alpha Numeric. CI ribbon alpha for comparison plots.
#' @param compare_legend_pos Character. Legend position for comparison plots.
#' @param compare_facet_ncol Integer. Default facet columns in comparison plots.
#' @param compare_facet_scales Character. Facet scales for comparison plots.
#' @param titles Named list. Override specific titles.
#' @param xlab Named list. Override specific x-axis labels.
#' @param ylab Named list. Override specific y-axis labels.
#' @return Invisible previous settings.
#' @export
#' @examples
#' # --- Set just one thing ---
#' acl_theme_set(base_theme = "theme_classic")
#' acl_theme_set(font_family = "Times New Roman")
#' acl_theme_set(compare_linewidth = 2)
#'
#' # --- Compare ACL vs ALSCL (default) ---
#' acl_theme_set(
#'   compare_colors = c("steelblue", "tomato"),
#'   compare_linetypes = c("solid", "longdash")
#' )
#'
#' # --- Compare same model, different parameters ---
#' acl_theme_set(
#'   compare_colors    = c("M=0" = "#2166AC", "M=0.3" = "#B2182B"),
#'   compare_linetypes = c("M=0" = "solid",   "M=0.3" = "dotted")
#' )
#' # plot_compare_ts(result_M0, result_M03,
#' #   model1_name = "M=0", model2_name = "M=0.3")
#'
#' # --- Or just override per-call, no global setting needed ---
#' # plot_compare_ts(result1, result2,
#' #   model1_name = "Run A", model2_name = "Run B",
#' #   colors = c("darkgreen", "purple"))
acl_theme_set <- function(base_theme = NULL, font_family = NULL, title_size = NULL,
                          title_hjust = NULL,
                          axis_title_size = NULL, axis_text_size = NULL,
                          strip_text_size = NULL, legend_text_size = NULL,
                          x_breaks = NULL, x_expand = NULL,
                          line_color = NULL, line_size = NULL,
                          se_color = NULL, se_alpha = NULL,
                          compare_colors = NULL, compare_linetypes = NULL,
                          compare_linewidth = NULL, compare_point_size = NULL,
                          compare_se_alpha = NULL, compare_legend_pos = NULL,
                          compare_facet_ncol = NULL, compare_facet_scales = NULL,
                          titles = NULL, xlab = NULL, ylab = NULL) {

  current <- getOption("acl.theme", .acl_defaults)
  old <- current

  if (!is.null(base_theme))           current$base_theme           <- base_theme
  if (!is.null(font_family))          current$font_family          <- font_family
  if (!is.null(title_size))           current$title_size           <- title_size
  if (!is.null(title_hjust))          current$title_hjust          <- title_hjust
  if (!is.null(axis_title_size))      current$axis_title_size      <- axis_title_size
  if (!is.null(axis_text_size))       current$axis_text_size       <- axis_text_size
  if (!is.null(strip_text_size))      current$strip_text_size      <- strip_text_size
  if (!is.null(legend_text_size))     current$legend_text_size     <- legend_text_size
  if (!is.null(x_breaks))            current$x_breaks             <- x_breaks
  if (!is.null(x_expand))            current$x_expand             <- x_expand
  if (!is.null(line_color))           current$line_color           <- line_color
  if (!is.null(line_size))            current$line_size            <- line_size
  if (!is.null(se_color))             current$se_color             <- se_color
  if (!is.null(se_alpha))             current$se_alpha             <- se_alpha
  if (!is.null(compare_colors))       current$compare_colors       <- compare_colors
  if (!is.null(compare_linetypes))    current$compare_linetypes    <- compare_linetypes
  if (!is.null(compare_linewidth))    current$compare_linewidth    <- compare_linewidth
  if (!is.null(compare_point_size))   current$compare_point_size   <- compare_point_size
  if (!is.null(compare_se_alpha))     current$compare_se_alpha     <- compare_se_alpha
  if (!is.null(compare_legend_pos))   current$compare_legend_pos   <- compare_legend_pos
  if (!is.null(compare_facet_ncol))   current$compare_facet_ncol   <- compare_facet_ncol
  if (!is.null(compare_facet_scales)) current$compare_facet_scales <- compare_facet_scales
  if (!is.null(titles))               current$titles               <- modifyList(current$titles, titles)
  if (!is.null(xlab))                 current$xlab                 <- modifyList(current$xlab, xlab)
  if (!is.null(ylab))                 current$ylab                 <- modifyList(current$ylab, ylab)

  options(acl.theme = current)
  invisible(old)
}


#' Reset ACL Plot Theme to Defaults
#'
#' @export
#' @examples
#' acl_theme_reset()
acl_theme_reset <- function() {
  options(acl.theme = .acl_defaults)
  cat("ACL plot theme reset to defaults.\n")
}


# ---------------------------------------------------------------------------
# Internal helpers (used by plot functions)
# ---------------------------------------------------------------------------

#' Get a title from the theme, with model name substitution
#' @param key Character. Title key.
#' @param m1 Optional model 1 name for substitution.
#' @param m2 Optional model 2 name for substitution.
#' @return Character string.
#' @keywords internal
.acl_title <- function(key, m1 = NULL, m2 = NULL) {
  result <- acl_theme("titles")[[key]]
  if (is.null(result)) return(key)
  if (!is.null(m1)) result <- gsub("\\{m1\\}", m1, result)
  if (!is.null(m2)) result <- gsub("\\{m2\\}", m2, result)
  result
}

#' @keywords internal
.acl_lab <- function(axis = "x", key) {
  labs <- if (axis == "x") acl_theme("xlab") else acl_theme("ylab")
  result <- labs[[key]]
  if (!is.null(result)) result else key
}

#' Build the standard ACL ggplot theme layer
#' @keywords internal
.acl_base_theme <- function(font_family = NULL, title_size = NULL,
                            axis_title_size = NULL, axis_text_size = NULL,
                            strip_text_size = NULL, legend_text_size = NULL,
                            base_theme = NULL, title_hjust = NULL) {
  bt  <- if (!is.null(base_theme))        base_theme       else acl_theme("base_theme")
  ff  <- if (!is.null(font_family))       font_family      else acl_theme("font_family")
  ts  <- if (!is.null(title_size))        title_size       else acl_theme("title_size")
  th  <- if (!is.null(title_hjust))       title_hjust      else acl_theme("title_hjust")
  ats <- if (!is.null(axis_title_size))   axis_title_size  else acl_theme("axis_title_size")
  atx <- if (!is.null(axis_text_size))    axis_text_size   else acl_theme("axis_text_size")
  sts <- if (!is.null(strip_text_size))   strip_text_size  else acl_theme("strip_text_size")
  lts <- if (!is.null(legend_text_size))  legend_text_size else acl_theme("legend_text_size")

  theme_fn <- switch(bt,
                     "theme_bw"        = ggplot2::theme_bw,
                     "theme_minimal"   = ggplot2::theme_minimal,
                     "theme_classic"   = ggplot2::theme_classic,
                     "theme_gray"      = ggplot2::theme_gray,
                     "theme_grey"      = ggplot2::theme_grey,
                     "theme_light"     = ggplot2::theme_light,
                     "theme_linedraw"  = ggplot2::theme_linedraw,
                     "theme_dark"      = ggplot2::theme_dark,
                     "theme_void"      = ggplot2::theme_void,
                     ggplot2::theme_bw
  )

  list(
    theme_fn(),
    ggplot2::theme(
      text             = ggplot2::element_text(family = ff),
      plot.title       = ggplot2::element_text(size = ts, hjust = th),
      axis.title       = ggplot2::element_text(size = ats),
      axis.text        = ggplot2::element_text(size = atx),
      strip.text       = ggplot2::element_text(size = sts),
      legend.text      = ggplot2::element_text(size = lts),
      legend.title     = ggplot2::element_text(size = lts)
    )
  )
}


#' Build x-axis scale with breaks and expand
#' @keywords internal
.acl_scale_x <- function(x_breaks = NULL, n_breaks = 10) {
  brks <- if (!is.null(x_breaks)) x_breaks else acl_theme("x_breaks")
  expd <- acl_theme("x_expand")
  if (is.null(expd)) expd <- c(0.01, 0.01)

  if (is.null(brks)) {
    ggplot2::scale_x_continuous(
      breaks = scales::pretty_breaks(n = n_breaks),
      expand = ggplot2::expansion(mult = expd)
    )
  } else {
    ggplot2::scale_x_continuous(
      breaks = brks,
      expand = ggplot2::expansion(mult = expd)
    )
  }
}


#' Get comparison plot defaults, merging global settings with per-call overrides
#'
#' Resolves colors/linetypes with this priority:
#' \enumerate{
#'   \item Per-call argument (e.g. \code{colors = c("blue","red")})
#'   \item Global theme with matching names (e.g. \code{compare_colors = c(ACL = "blue", ...)})
#'   \item Global theme positional (first value -> model1, second -> model2)
#'   \item Hard-coded fallback
#' }
#'
#' @param colors Per-call override. Named or unnamed vector of 2.
#' @param linetypes Per-call override.
#' @param linewidth Per-call override.
#' @param point_size Per-call override.
#' @param se_alpha Per-call override.
#' @param legend_pos Per-call override.
#' @param ncol Per-call override.
#' @param scales Per-call override.
#' @param m1 Model 1 name.
#' @param m2 Model 2 name.
#' @return A named list of resolved settings.
#' @keywords internal
.acl_compare_defaults <- function(colors = NULL, linetypes = NULL,
                                  linewidth = NULL, point_size = NULL,
                                  se_alpha = NULL, legend_pos = NULL,
                                  ncol = NULL, scales = NULL,
                                  m1 = "Model 1", m2 = "Model 2") {

  # Helper: resolve a 2-element vector (color or linetype) to named c(m1=..., m2=...)
  resolve_pair <- function(user_val, theme_key, fallback) {
    # Priority 1: per-call argument
    val <- user_val
    # Priority 2: global theme
    if (is.null(val)) val <- acl_theme(theme_key)
    # Priority 3: fallback
    if (is.null(val) || length(val) < 2) val <- fallback

    # If named and both model names match, use by name
    if (!is.null(names(val)) && m1 %in% names(val) && m2 %in% names(val)) {
      return(stats::setNames(val[c(m1, m2)], c(m1, m2)))
    }
    # Otherwise use positionally: first = model1, second = model2
    stats::setNames(val[1:2], c(m1, m2))
  }

  list(
    colors     = resolve_pair(colors,    "compare_colors",    c("#2166AC", "#B2182B")),
    linetypes  = resolve_pair(linetypes, "compare_linetypes", c("solid", "dashed")),
    linewidth  = if (!is.null(linewidth))  linewidth  else (acl_theme("compare_linewidth")  %||% 1.2),
    point_size = if (!is.null(point_size)) point_size else (acl_theme("compare_point_size") %||% 2),
    se_alpha   = if (!is.null(se_alpha))   se_alpha   else (acl_theme("compare_se_alpha")   %||% 0.15),
    legend_pos = if (!is.null(legend_pos)) legend_pos else (acl_theme("compare_legend_pos") %||% "bottom"),
    ncol       = if (!is.null(ncol))       ncol       else (acl_theme("compare_facet_ncol") %||% 2),
    scales     = if (!is.null(scales))     scales     else (acl_theme("compare_facet_scales") %||% "free_y")
  )
}

# base R fallback for %||% (available in R >= 4.0 but define for safety)
`%||%` <- function(x, y) if (is.null(x)) y else x


#' Fix length bin labels
#' @keywords internal
.acl_fix_len_labels <- function(len_label, prefix = "Length bin") {
  LengthGroup <- len_label
  n_bins <- length(LengthGroup)
  if (n_bins >= 2) {
    first_parts <- strsplit(as.character(LengthGroup[1]), "-")[[1]]
    last_parts  <- strsplit(as.character(LengthGroup[n_bins]), "-")[[1]]
    if (length(first_parts) == 2) LengthGroup[1] <- paste0("<", first_parts[2])
    if (length(last_parts)  == 2) LengthGroup[n_bins] <- paste0(">", last_parts[1])
  }
  paste(prefix, LengthGroup)
}


#' Internal: Build ggplot2 theme object from acl_theme()
#'
#' This is the function called by plot functions via
#' \code{get("acl_get_theme", envir = asNamespace("ACL"))()}.
#'
#' @return A list of ggplot2 theme layers.
#' @keywords internal
acl_get_theme <- function() {
  .acl_base_theme()
}
