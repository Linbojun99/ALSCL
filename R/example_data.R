#' Example Catch-at-Length Survey Data
#'
#' A list containing three data frames for demonstrating ACL and ALSCL model
#' fitting. The data represent anonymized survey catch-at-length observations
#' across 9 length bins and 21 years.
#'
#' @format A list with three elements:
#' \describe{
#'   \item{data.CatL}{Survey catch-at-length. A data frame with 9 rows
#'     (length bins) and 22 columns (LengthBin + 21 years, 2001--2021).}
#'   \item{data.wgt}{Weight-at-length. Same dimensions as \code{data.CatL}.
#'     Values are constant across years.}
#'   \item{data.mat}{Maturity-at-length (proportion mature). Same dimensions
#'     as \code{data.CatL}. Values are constant across years.}
#' }
#'
#' All three data frames share the same \code{LengthBin} labels:
#' \code{"0-20"}, \code{"21-25"}, \code{"26-30"}, \code{"31-35"},
#' \code{"36-40"}, \code{"41-45"}, \code{"46-50"}, \code{"51-55"},
#' \code{"56-60"}.
#'
#' @usage data(example_data)
#'
#' @examples
#' data(example_data)
#' data.CatL <- example_data$data.CatL
#' data.wgt  <- example_data$data.wgt
#' data.mat  <- example_data$data.mat
#'
"example_data"
