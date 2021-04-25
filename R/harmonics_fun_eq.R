#' Equidistant time series modelling
#'
#' This function enables the user to model different number of cycles per year and calculate a desired number of artificial, equidistant time steps.
#'
#' @param user_vals A vector with numeric values.
#'
#' @param user_dates A Vector with Date objects. See \code{\link[base]{as.Date}}.
#'
#' @param harmonic_deg Numeric. The number of cycles per year (harmonic degree)
#'                     that should be modelled.
#'
#' @param ref_date (optional) A Date object. Default is 1970-01-01.
#'
#' @param new_dates (optional) A vector with Date objects. See \code{\link[base]{as.Date}}. These are the desired artifical time steps.
#'
#' @param return_vals (optional) String. Defines returned variables as one of "fitted","trend", "h1", "h2", (...and so on for chosen number of harmonics) or "all" to get a list of all of them. Default is "fitted".
#'
#' @param print_variance (optional) Boolean. If TRUE, prints explained variance of every harmonic. Default is FALSE.
#'
#' @return A numeric vector with the fitted values or a list of vectors with fitted values, trend and all individual harmonics (in this order).
#'
#' @details To calculate the harmonic fitted curve of a periodic signal,
#'          ordinary least squares regressions are computed using coupled
#'          sine and cosine curves on time-series data. The underlying algorithm
#'          is based on Shumway & Stoffer (2017) equations 4.1 – 4.2. Based on
#'          the derived function,the input of new dates allows to calculate
#'          artificial time steps other than the input dates.
#'
#' @references Shumway, R. H., & Stoffer, D. S. (2017). Time series analysis and its
#'             applications: with R examples. Springer.
#'
#' @examples
#'
#' @export
#'
harmonics_fun_eq <- function(user_vals,
                                user_dates,
                                harmonic_deg = 3,
                                ref_date = as.Date("1970-01-01", format = "%Y-%m-%d"),
                                new_dates =NULL,
                                return_vals="fitted",
                                print_variance=FALSE,
                                #maxgap = NULL
                                ){

  ### For every missing output parameter set the default ###
  if (missing(user_vals)){
    stop("Values must be provided.")
  }
  if (missing(user_dates)){
    stop("Dates must be provided.")
  }
  # Check if dates are in "Date"-format
  if (class(user_dates) != "Date"){
    stop("Dates must be provided as 'Date' objects.")
  }
  if (missing(harmonic_deg)){
    stop("Harmonic degree must be provided.")
  }
  if (missing(ref_date)){
    ref_date <- as.Date("1970-01-01")
  } else if (class(ref_date) != "Date"){
    stop("Reference date must be provided as a 'Date' object.")
  }

  if (all(is.na(user_vals))) {
    print("User values consist of NA values only. Output is same as input")
    return(user_vals)

    # Otherwise apply harmonic analysis
  } else {

    # if (! is.null(maxgap)){
    #   df <- interpolate_ts(dates=user_dates, values=user_vals, maxgap= maxgap)
    #   user_dates = df[1,]
    #   user_vals = df[,2]
    # }

    df_for_reg <- get_df(user_vals, user_dates, harmonic_deg, ref_date)
    coefficients <- apply_regression(df_for_reg, harmonic_deg, print_variance)

    # output shall be modelled time series with original input time steps
    if (is.null(new_dates)){
      fitted <- calculate_fitted(df=df_for_reg, coefs=coefficients, harmonic_deg, return_vals, user_vals, print_variance)}

    # output shall be modelled equidistant time series with new time steps
    if (! is.null(new_dates)) {
      # calculate radians,sines and cosines
      # a dummy is assigned to user_vals to preserve the data frame structure that is used in the case of modelling with original dates
      df_new <- get_df(user_vals=rep(0, length(new_dates)), user_dates=new_dates, harmonic_deg, ref_date)
      # calculate fitted values basedon new dates
      fitted <- calculate_fitted(df=df_new, coefs=coefficients, harmonic_deg, return_vals, user_vals, print_variance)
    }

    return (fitted)
  }
}

