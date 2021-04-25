#'Internal functions
#' NULL declaration to suppres R CMD CHECK warning related to tidyverse syntax
#' @keywords internal
#' @noRd

# not yet implemented
# interpolate_ts <- function(dates,values, maxgap=maxgap){
#   df <-  data.frame(dates, values)
#   for(i in 1:(length(dates)-1)){
#     gap <-  as.numeric(difftime(dates[i+1], dates[i], units = "days"))
#     if (gap > maxgap){
#       new_date <-  dates[i] + gap/2
#       new_val <- values[i]+ ((values[i+1]-values[i])/2)
#       df[nrow(df)+1,1] <- as.POSIXct(new_date, format="%Y-%m-%d")
#       df[nrow(df),2] <- new_val
#     }
#   }
#   df <- dplyr::arrange(df, df$dates)
#   return(df)
# }

get_df <- function(user_vals, user_dates, harmonic_deg, ref_date){
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
  #print(colnames(df_for_reg))
  return(df_for_reg)
}


apply_regression <- function(df_for_reg, harmonic_deg, print_variance){
  # Ordinary Least Squares Regression
  my_reg <- stats::lm(formula = user_vals ~ ., data = df_for_reg)
  # Get coefficients (constant is intercept value)
  # For coefficient and radians it's easy ...
  cons_coef <- my_reg$coefficients[1]
  t_coef <- my_reg$coefficients[3]

  if (print_variance == TRUE){
    print(paste0("overall r-squared:", summary(my_reg)$r.squared))
  }

  # ... for the sines and cosines selecting the right columns is a little more complicated
  # Start with 4 because the first three values are the ndvi, my_cons and my_radians
  # First define the start and stop column for the sines and cosines ...
  sin_start <- 4
  sin_end <- 4 + harmonic_deg - 1
  cos_start <- 4 + harmonic_deg
  cos_end <- 4 + harmonic_deg + harmonic_deg -1
  # ... and then subset the data accordingly
  sin_coef <- my_reg$coefficients[c(sin_start:sin_end)]
  #print(sin_coef)
  cos_coef <- my_reg$coefficients[c(cos_start:cos_end)]
  return(list(cons_coef, t_coef, sin_coef, cos_coef))
}


variance <- function(i, user_vals, cos_coef, sin_coef){
  ((length(user_vals)*sqrt(cos_coef[i]^2 +sin_coef[i]^2))/(2*(length(user_vals)-1) * sd(user_vals, na.rm=T)))#/100
}



calculate_fitted <- function(df, coefs, harmonic_deg, return_vals, user_vals, print_variance){
  ### Calculate fitted values ###

  # multiply independent variables with the coefficients
  df[,2] <- df[,2] * coefs[[1]]
  df[,3] <- df[,3] * coefs[[2]]
  # for loop multiplying the factor for each harmonic degree
  for (k in 1:harmonic_deg){
    # for the sines define i + 3 because the first sine column is at the 4th position
    df[,k + 3] <- df[,k + 3] * coefs[[3]][k]
    # for the cosines define i + 3 + harmonic_deg to get to the first cosine column
    df[,k + 3 + harmonic_deg] <- df[,k + 3 + harmonic_deg] * coefs[[4]][k]
    # remove reduntant variables
    rm(k)
  }
  # calculate sum (fitted value) of the multplied independent variables
  fitted <- rowSums(df[,c(2:ncol(df))], na.rm = TRUE)
  trend = df[,2] + df[,3]

  for (i in 1:harmonic_deg){
    assign(paste0("h",i), df[,i+3] + df[,i+3+harmonic_deg])
    if (print_variance == T){
      print(paste0("explained variance h", i,": ", assign(paste0("var", i), variance(i=i,user_vals,cos_coef=coefs[[4]], sin_coef=coefs[[3]]))))
    }

  }

  if (return_vals == "all"){
    return(append(list(fitted, trend), lapply(seq(1,harmonic_deg), function(x) eval(parse(text = paste0("h", x))))))
  } else {
    return(eval(parse(text = return_vals)))
  }
}
