#' Harmonic modelling
#'
#' This function enables the user to model different number of cycles per year
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
#' @return A numeric vector with the fitted values.
#'
#' @details To calculate the harmonic fitted curve of a periodic signal,
#'          ordinary least squares regressions are computed using coupled
#'          sine and cosine curves on time-series data. The underlying algorithm
#'          is based on Shumway & Stoffer (2017) equations 4.1 – 4.2.
#'
#' @references Shumway, R. H., & Stoffer, D. S. (2017). Time series analysis and its
#'             applications: with r examples. Springer.
#'
#' @examples
#'
#' TODO
#'
#'
#'
#' @export

harmonictest <- function(user_vals, user_dates, harmonic_deg = 3, ref_date = as.Date("1970-01-01", format = "%Y-%m-%d"), new_dates =NULL){

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

  # If user vals are only NA, output is same as input

  if (all(is.na(user_vals))) {
    print("User values consist of NA values only. Output is same as input")
    return(user_vals)

    # Otherwise apply harmonic analysis
  } else {

    ### Calculate difference to ref_date in radians ###

    # Start for loop
    for (i in 1:length(user_dates)){
      current_date <- user_dates[i]
      # Calculate the difference in days
      current_diff_days <- as.numeric(difftime(current_date, ref_date), units="days")
      # Convert to years
      current_diff_years <- current_diff_days/365.25
      # Convert to radians
      current_diff_radians <- current_diff_years * 2 * pi
      if (i == 1){
        my_radians <- current_diff_radians
      } else {
        my_radians <- c(my_radians, current_diff_radians)
      }
      rm(i, current_date, current_diff_days, current_diff_years, current_diff_radians)
    }

    ### Caculate sines and cosines ###

    # Define constant
    my_cons <- rep(1, times = length(my_radians))
    # Create sines and cosines data frames with one column for each harmonic
    my_sin <- data.frame(matrix(nrow = length(my_radians), ncol = harmonic_deg))
    my_cos <- data.frame(matrix(nrow = length(my_radians), ncol = harmonic_deg))
    # Create names for df
    sin_names <- rep("sin_", times=harmonic_deg)
    cos_names <- rep("cos_", times=harmonic_deg)
    my_seq <- as.character(seq(1,harmonic_deg, by=1))
    sin_names <- paste(sin_names, my_seq, sep="")
    cos_names <- paste(cos_names, my_seq, sep="")
    # Add names to df
    names(my_sin) <- sin_names
    names(my_cos) <- cos_names
    # Calculate sines and cosines for each harmonic
    for (j in 1:harmonic_deg){
      # Calculate current sines and cosines by multiplying the radians with the current
      # degree of harmonic and then apply the sine/cosine function
      current_sines <- sin(my_radians * j)
      current_cosines <- cos(my_radians * j)
      # Fill data frames with values
      my_sin[,j] <- current_sines
      my_cos[,j] <- current_cosines
      # remove redundant variables
      rm(j, current_cosines, current_sines)
    }
    # Create df from the dependent and all independent variables
    df_for_reg <- cbind(user_vals, my_cons, my_radians, my_sin, my_cos)

    ### Apply Ordinary Least Squares Regression ###

    # Ordinary Least Squares Regression
    my_reg <- stats::lm(formula = user_vals ~ ., data = df_for_reg)
    # Get coefficients (constant is intercept value)
    # For coefficient and radians it's easy ...
    cons_coef <- my_reg$coefficients[1]
    t_coef <- my_reg$coefficients[3]
    # ... for the sines and cosines selecting the right columns is a little more complicated
    # Start with 4 because the first three values are the ndvi, my_cons and my_radians
    # First define the start and stop column for the sines and cosines ...
    sin_start <- 4
    sin_end <- 4 + harmonic_deg - 1
    cos_start <- 4 + harmonic_deg
    cos_end <- 4 + harmonic_deg + harmonic_deg -1
    # ... and then subset the data accordingly
    sin_coef <- my_reg$coefficients[c(sin_start:sin_end)]
    cos_coef <- my_reg$coefficients[c(cos_start:cos_end)]


    #return(c(cons_coef, t_coef, sin_coef, cos_coef))

    if (is.null(new_dates)){
      ### Calculate fitted values ###

      # multiply independent variables with the coefficients
      df_for_reg[,2] <- df_for_reg[,2] * cons_coef
      df_for_reg[,3] <- df_for_reg[,3] * t_coef
      # for loop multiplying the factor for each harmonic degree
      for (k in 1:harmonic_deg){
        # for the sines define i + 3 because the first sine column is at the 4th position
        df_for_reg[,k + 3] <- df_for_reg[,k + 3] * sin_coef[k]
        # for the cosines define i + 3 + harmonic_deg to get to the first cosine column
        df_for_reg[,k + 3 + harmonic_deg] <- df_for_reg[,k + 3 + harmonic_deg] * cos_coef[k]
        # calculate individual harmonics
        #h1 = df_for_reg[,k +3] + df_for_reg[,k+3+harmonic_deg]

        # remove reduntant variables
        rm(k)
      }
      # calculate sum (fitted value) of the multplied independent variables
      fitted <- rowSums(df_for_reg[,c(2:ncol(df_for_reg))], na.rm = TRUE)

      #if (return_harmonics == T){
      #trend = df_for_reg[,2] + df_for_reg[,3]
      # h1 = df_for_reg[,4] + df_for_reg[,7]
      # h2 = df_for_reg[,5] + df_for_reg[,8]
      # h3 = df_for_reg[,6] + df_for_reg[,9]
      #  return (list(fitted, c(trend, h1, h2, h3)))
      # }

    }


    # equidistant time series

    if (! is.null(new_dates)) {
      for (i in 1:length(new_dates)){
        current_date <- new_dates[i]
        # Calculate the difference in days
        current_diff_days <- as.numeric(difftime(current_date, ref_date), units="days")
        # Convert to years
        current_diff_years <- current_diff_days/365.25
        # Convert to radians
        current_diff_radians <- current_diff_years * 2 * pi
        if (i == 1){
          new_radians <- current_diff_radians
        } else {
          new_radians <- c(new_radians, current_diff_radians)
        }
        rm(i, current_date,current_diff_radians)
      }


      ### Caculate sines and cosines ###

      # Define constant
      new_cons <- rep(1, times = length(new_radians))
      # Create sines and cosines data frames with one column for each harmonic
      new_sin <- data.frame(matrix(nrow = length(new_radians), ncol = harmonic_deg))
      new_cos <- data.frame(matrix(nrow = length(new_radians), ncol = harmonic_deg))
      names(new_sin) <- sin_names
      names(new_cos) <- cos_names
      # Calculate sines and cosines for each harmonic
      for (j in 1:harmonic_deg){
        # Calculate current sines and cosines by multiplying the radians with the current
        # degree of harmonic and then apply the sine/cosine function
        current_sines <- sin(new_radians * j)
        current_cosines <- cos(new_radians * j)
        # Fill data frames with values
        new_sin[,j] <- current_sines
        new_cos[,j] <- current_cosines
        # remove redundant variables
        rm(j, current_cosines, current_sines)
      }
      # Create df from the dependent and all independent variables
      df_new <- cbind(new_cons, new_radians, new_sin, new_cos)


      ### Calculate fitted values ###

      # multiply independent variables with the coefficients
      df_new[,1] <- df_new[,1] * cons_coef
      df_new[,2] <- df_new[,2] * t_coef
      # for loop multiplying the factor for each harmonic degree
      for (k in 1:harmonic_deg){
        # for the sines define i + 3 because the first sine column is at the 4th position
        df_new[,k + 2] <- df_new[,k + 2] * sin_coef[k]
        # for the cosines define i + 3 + harmonic_deg to get to the first cosine column
        df_new[,k + 2 + harmonic_deg] <- df_new[,k + 2 + harmonic_deg] * cos_coef[k]
        # remove reduntant variables
        rm(k)
      }
      # calculate sum (fitted value) of the multplied independent variables
      fitted <- rowSums(df_new[,c(1:ncol(df_new))], na.rm = TRUE)
    }
  }

  #if (return_coefs == T){
  #  coefs = c(cons_coef, t_coef, sin_coef, cos_coef)
  #  return(list(fitted, coefs))
  #} else {return(fitted)}
  return(fitted)
}
