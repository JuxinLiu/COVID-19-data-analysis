# February 11, 2023
# lag selection based on MPSE (mean prediction square errors for DAILY counts)
# cumulative sum was reset to be the new daily counts on the start date
# Piecewise linear regression implemented by segmented()


#### Inputs of the prediction function pwl_pred()
## df is the input data with 5 columns: 
## Column 1: Date (in date format)
## Column 2: daily counts of Y (e.g., deaths)
## Column 3: cumulative counts of Y
## Column 4: daily counts of X (e.g., cases)
## Column 5: cumulative counts of X
## alpha: 1-alpha is the confidence level for confidence bands
## bp_index: the indices of breakpoints (initial values)
## min_lead: minimal value of the lag
## max_lead: maximal value of the lag
## CI: If CI == TRUE, both pointwise and joint confidence bands for out-of-sample predictions are calculated. Otherwise not.


pwl_pred <- function(df, alpha = 0.1, CI = FALSE, bp_index = bp_index, min_lead = 5,
                   max_lead = 21) {
  
  library(timetk)
  library(tidyverse)
  library(tvReg)
  library(dplyr)
  library(ggplot2)
  library(segmented)

  names(df) <- c("Date", "y", "cumulative_y", "x", "cumulative_x")
  # reset the cumulative deaths/cases using the start_date
  df$cumulative_y_reset <- df$cumulative_y - df$cumulative_y[1] 
  df$cumulative_x_reset <- df$cumulative_x - df$cumulative_x[1]
  
  # Return the optimal choice of lead for piecewise linear regression
  pwl_lag <- function(df, bp_index){
    
    #Create a variable to store Mean Squared Prediction Error for each possible lead value
    MSPE_List <- c()
    
    #Loop through all possible leads and record the MSPE
    for (n_ahead in min_lead:max_lead) {
      
      df_pwl <- df %>%
        tk_augment_lags(c(cumulative_y_reset, y), .lags = -n_ahead,
                        .names = c("y_best_lead_pwl", "daily_y_lead"))
      
      D1 <- max(df_pwl$Date) -  2 * max_lead
      
      df_pwl_insample <- df_pwl %>% filter (Date <= D1)
      
      lm_insample <- lm(y_best_lead_pwl ~ cumulative_x_reset, data = df_pwl_insample)
      
      
      #pwl_insample <- segmented(lm_insample,seg.Z = ~ cumulative_x_reset, npsi = nbreakpoint)
      bp_value <- df_pwl_insample$cumulative_x_reset[bp_index]
      
      pwl_insample <- segmented(lm_insample, seg.Z = ~ cumulative_x_reset, psi = bp_value)
      
      newdf1 <- df_pwl %>% filter (Date > D1  & Date  <= D1 + n_ahead) %>% 
        dplyr::select (c("cumulative_x_reset", "y_best_lead_pwl", "daily_y_lead"))
      
      #get predictions for cases outside the insample
      df_pwl_pred1 <- predict(pwl_insample,newdata = newdf1[,1])
      
      #residual <- c(tvLM_insample$residuals, tvLM_df_pred1 - newdf1$death_best_lead_tv)
      #residual <- tvLM_df_pred1 - newdf1$death_best_lead_tv
      
      daily_pred <- tail(df_pwl_pred1, -1) - head(df_pwl_pred1, -1)
      daily_residual <- daily_pred - newdf1$daily_y_lead[-1]
      MSPE <- mean((daily_residual)^2)
      #print(c(n_ahead, MSPE))
      MSPE_List <- c(MSPE_List, MSPE)
    }
    
    #Return the lead with the smallest MSPE
    ans <- which.min(MSPE_List) + (min_lead - 1)
    #print (paste("best lead for pwl:", ans))
    return(ans)
    
  } # end of pwl_lag function
  
  
  n_ahead_pwl <- pwl_lag(df, bp_index)

  
  # including death counts for the best lead number
  df_pwl <- df %>%
    tk_augment_lags(cumulative_y_reset, .lags = -n_ahead_pwl, .names = "y_best_lead_pwl")
  names(df_pwl) <- names(df_pwl) %>% str_replace_all("lag-|lag", "lead")
  
  D1_pwl <- max(df_pwl$Date) - 2*n_ahead_pwl 
  D2_pwl <- max(df_pwl$Date) - n_ahead_pwl + 1
  
  
  df_pwl_insample <- df_pwl %>% filter (Date <= D1_pwl )
  
  lm_insample <- lm(y_best_lead_pwl ~ cumulative_x_reset,
                    data = df_pwl_insample)
  
breakpoint_values <- df_pwl_insample$cumulative_x_reset[bp_index]
  
  
  pwl_insample <- segmented(lm_insample,seg.Z =~cumulative_x_reset, psi = breakpoint_values)
  
  #pwl_insample <- segmented(lm_insample,seg.Z =~cumulative_x_reset)
  newdf1 <- df_pwl %>% filter (Date > D1_pwl & Date < D2_pwl) %>% dplyr::select (c("cumulative_x_reset", "y_best_lead_pwl"))
  
  df_pwl_pred1 <- predict(pwl_insample, newdata = newdf1[,1], interval = "prediction", level = 1 - alpha)
  
  conf_band_lo <- c(rep(NA, length(pwl_insample$fitted)), df_pwl_pred1[,2])
  conf_band_hi <- c(rep(NA, length(pwl_insample$fitted)), df_pwl_pred1[,3])
  
  pred_pwl <- c(pwl_insample$fitted.values, df_pwl_pred1[,1])
  
  n <- dim(df_pwl)[1] # same as dim(df_tv)[1]
  opt_pwl <- data.frame(pred_pwl, conf_band_lo, conf_band_hi)
  
  names(opt_pwl) <- c("predicted_pwl", "pwl_lower", "pwl_upper")
  
  npred <- dim(opt_pwl)[1]
  
  opt_pwl$pred_daily_pwl <- c(NA, tail(opt_pwl$predicted_pwl, npred - 1) - head (opt_pwl$predicted_pwl, npred - 1))
  
  opt_pwl$date <- head(df_pwl$Date, n - n_ahead_pwl) + n_ahead_pwl
  
  #opt_pwl$obs_cum <- df_pwl$cumulative_y_reset[df_pwl$Date %in% opt_pwl$date]
  opt_pwl$observed_cum <- head(df_pwl$y_best_lead_pwl, n - n_ahead_pwl)
  
  opt_pwl$observed_daily <- df_pwl$y[df_pwl$Date %in% opt_pwl$date]
  
  mspe <- mean((tail(opt_pwl$observed_daily, n_ahead_pwl) - tail(opt_pwl$pred_daily_pwl, n_ahead_pwl))^2, na.rm = TRUE)
  
#Generate plots of the output
 
    gg_cumulative <- opt_pwl %>%
    pivot_longer(cols = where(is.numeric)) %>%
    filter(name %in% c("observed_cum", "predicted_pwl", "pwl_lower", "pwl_upper", "date")) %>%
    ggplot(aes(date, value, color = name, linetype = name)) +
    geom_line() +
    geom_vline(xintercept = D2_pwl, color = "orange") +
    
    scale_linetype_manual(values=c("solid", "dashed", "dotted", "dotted")) +
    scale_color_manual(values=c('black', 'red', 'blue','blue')) +
    
    labs(
      title = paste("Reported vs. Predicted Cumulative Deaths,", region, sep = " "),
      subtitle = "Inputs = Cumulative Cases",
      x = "Date",
      y = "Count",
      caption = paste("pwl lag =", n_ahead_pwl, "days", sep = " ")
    )
  
  gg_daily <- opt_pwl %>%
    pivot_longer(cols = where(is.numeric)) %>%
    filter(name %in% c("observed_daily", "pred_daily_pwl")) %>%
    ggplot(aes(date, value, color = name, linetype = name, size = name)) +
    geom_point() + geom_line() + 
    geom_vline(xintercept = D2_pwl, color = "orange") +
    scale_linetype_manual(values=c("solid", "dashed")) +
    scale_color_manual(values=c('black', 'red')) +
    scale_size_manual(values = c(1, 1)) + 
    labs(
      title = paste("Reported vs. Predicted Daily Counts,", region, sep = " "),
      subtitle = "Inputs = Cumulative Cases",
      x = "Date",
      y = "Count",
      caption = paste("pwl lag =", n_ahead_pwl, "days", sep = " ")
    )

  return(list(daily = gg_daily, cumulative = gg_cumulative, mspe = mspe, 
              out = opt_pwl, lag = n_ahead_pwl, D2 = D2_pwl))
  
}
 
