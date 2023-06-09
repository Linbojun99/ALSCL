% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_VB.R
\name{plot_VB}
\alias{plot_VB}
\title{Plot the Von Bertalanffy growth function}
\usage{
plot_VB(
  model_result,
  age_range = c(1, 25),
  line_size = 1.2,
  line_color = "red",
  line_type = "solid",
  se = FALSE,
  se_color = "red",
  se_alpha = 0.2,
  text_color = "black",
  text_size = 5
)
}
\arguments{
\item{model_result}{A list that contains model output. The list should have a "report" component which contains "Linf", "vbk" and "t0" components.}

\item{age_range}{Numeric vector of length 2, defining the range of ages to consider.}

\item{line_size}{Numeric. The thickness of the line in the plot. Default is 1.2.}

\item{line_color}{Character. The color of the line in the plot. Default is "black".}

\item{line_type}{Character. The type of the line in the plot. Default is "solid".}

\item{se}{Logical. Whether to calculate and plot standard error as confidence intervals. Default is FALSE.}

\item{se_color}{Character. The color of the confidence interval ribbon. Default is "blue".}

\item{se_alpha}{Numeric. The transparency of the confidence interval ribbon. Default is 0.2.}

\item{text_color}{Character. The color of the text in the plot. Default is "black".}

\item{text_size}{Numeric. The thickness of the text in the plot. Default is 5.}
}
\value{
A ggplot object representing the plot.
}
\description{
This function takes Linf, k, and t0 as input, creates a sequence of ages, applies the VB function to each age, and then plots the result using ggplot2.
}
\examples{
\dontrun{
plot_VB(model_result, age_range = c(1, 25))
}
}
