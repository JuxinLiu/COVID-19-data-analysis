# February 13, 2023
# lag selection based on MPSE (mean prediction square errors for daily prediction)
# cumulative sum was reset to be the new daily counts on the start date
# local polynomial


#############################################
## local polynomial regression
#############################################
## df is the input data with 5 columns: 
## Column 1: Date (in date format)
## Column 2: daily counts of Y (e.g., deaths)
## Column 3: cumulative counts of Y
## Column 4: daily counts of X (e.g., cases)
## Column 5: cumulative counts of X
## alpha: 1-alpha is the confidence level for confidence bands
## method: "ll" for local linear and "lc" for local constant
## min_lead: minimal value of the lag
## max_lead: maximal value of the lag
## CI: If CI == TRUE, both pointwise and joint confidence bands for out-of-sample predictions are calculated. Otherwise not.

localpoly <- function(df, alpha = 0.1, method = "ll", min_lead = 5, max_lead = 21, CI = FALSE) {
  
  library(tidyverse)
  library(tvReg)
  library(dplyr)
  library(ggplot2)
  
  names(df) <- c("Date", "y", "cumulative_y", "x", "cumulative_x")
  # reset the cumulative y/cases using the start_date
  df$cumulative_y_reset <- df$cumulative_y - df$cumulative_y[1] 
  df$cumulative_x_reset <- df$cumulative_x - df$cumulative_x[1]
  
  # Return the optimal choice of lead for Tvreg
  tv <- function(df){
  
    #Create a variable to store Mean Squared Prediction Error for each possible lead value
    MSPE_List <- c()
    
    #Loop through all possible leads and record the MSPE
    for (n_ahead_tv in min_lead:max_lead) {
      
      df_tv <- df %>%
        tk_augment_lags(c(cumulative_y_reset, y), .lags = -n_ahead_tv, .names = c("y_best_lead_tv", "daily_y_lead"))
      
      names(df_tv) <- names(df_tv) %>% str_replace_all("lag-|lag", "lead")
      
      D1 <- max(df_tv$Date) -  2 * max_lead
      
      df_tv_insample <- df_tv %>% filter (Date <= D1)
      
      tvLM_insample <- tvLM(y_best_lead_tv ~ cumulative_x_reset, 
                            est = method, 
                            data = df_tv_insample,
                            #tkernel = "Gaussian",
                            singular.ok = FALSE)
      
      newdf1 <- df_tv %>% filter (Date > D1 & Date <= D1 + n_ahead_tv) %>% 
        dplyr::select (c("cumulative_x_reset", "y_best_lead_tv", "daily_y_lead"))
      
      #get predictions for cases outside the insample
      tvLM_df_pred1 <- forecast(tvLM_insample, newdata = as.matrix(newdf1[,1]), n.ahead = n_ahead_tv)
      
      #residual <- c(tvLM_insample$residuals, tvLM_df_pred1 - newdf1$y_best_lead_tv)
      #residual <- tvLM_df_pred1 - newdf1$y_best_lead_tv
      
      daily_pred <- tail(tvLM_df_pred1, -1) - head(tvLM_df_pred1, -1)
      daily_residual <- daily_pred - newdf1$daily_y_lead[-1]
      MSPE <- mean((daily_residual)^2)
      
      MSPE_List <- c(MSPE_List, MSPE)
    }
    
    #Return the lead with the smallest MSPE
    ans <- which.min(MSPE_List) + (min_lead - 1)
    print (paste("best lead for tv:", ans))
    return(ans)
  }
  
  n_ahead_tv <- tv(df)
  
  #predictions based on TvReg
   df_tv <- df %>%
      tk_augment_lags(cumulative_y_reset, .lags = -n_ahead_tv, .names = "y_best_lead_tv")
    names(df_tv) <- names(df_tv) %>% str_replace_all("lag-|lag", "lead")
    
    D1 <- max(df_tv$Date) - 2*n_ahead_tv 
    D2 <- max(df_tv$Date) - n_ahead_tv + 1
    
    df_tv_insample <- df_tv %>% filter (Date <= D1 )
    
    tvLM_insample <- tvLM(y_best_lead_tv ~ cumulative_x_reset, 
                          est = method, 
                          data = df_tv_insample,
                          #tkernel = "Epa",
                          singular.ok = FALSE)
    
    newdf1 <- df_tv %>% filter (Date > D1 & Date < D2) %>% dplyr::select (c("cumulative_x_reset", "y_best_lead_tv"))
    
    tvLM_df_pred1 <- forecast (tvLM_insample, newdata = as.matrix(newdf1[,1]), n.ahead = n_ahead_tv)
    
    
    pred_tv <- c(tvLM_insample$fitted, tvLM_df_pred1)
    
    n <- dim(df_tv)[1] # same as dim(df_tv)[1]
    opt_tv <- data.frame(pred_tv, head(df_tv$y_best_lead_tv, n - n_ahead_tv))
    names(opt_tv) <- c("predicted_tv", "observed_tv")
    n_tv <- length(pred_tv)
    opt_tv$pred_daily_tv <- c(NA, pred_tv[2:n_tv] - pred_tv[1:(n_tv - 1)])
    
    opt_tv$date <- head(df_tv$Date, n - n_ahead_tv) + n_ahead_tv ### date being moved forward
    
    opt_tv$observed_daily <- df_tv$y[df_tv$Date %in% opt_tv$date]
    
    mspe <- mean((tail(opt_tv$observed_daily, n_ahead_tv) - tail(opt_tv$pred_daily_tv, n_ahead_tv))^2, na.rm = TRUE)  
   
     if (CI == TRUE){
       # Pointwise confidence bands
      point_conf <- function (data, n_ahead, B = 200, alpha = 0.05) { # df after the best lead being selected
      #B = 1000 ## number of bootstrap draws
      data_insample <- head (data, dim(data)[1] - n_ahead)
      nobs <- dim(data_insample)[1]
      object <- tvLM(y_best_lead_tv ~ cumulative_x_reset, 
                     data = data_insample,
                     est = method) 
      
      residuals_raw <- object$residuals - mean(object$residuals)
      residuals_b <- matrix(sample (residuals_raw * rnorm(B * nobs), size = B * nobs, replace = TRUE),
                            nrow = B, ncol = nobs)
      y_b <- matrix(object$fitted, nrow = B, ncol = nobs, byrow = TRUE) + residuals_b ## synthetic y
      
      ## out-of-sample prediction 1 with death counts observable
      prediction1 <- matrix(NA, nrow = B, ncol = n_ahead)
      #prediction_daily <- matrix(NA, nrow = B, ncol = n_ahead - 1)
      totobs <- nobs + n_ahead
      
      newdf <- data %>% tail(n_ahead) %>% dplyr::select (c("cumulative_x_reset", "y_best_lead_tv"))
      newdata <- cbind(rep(1, n_ahead), as.matrix(newdf[, 1]))
      pred_raw <- forecast (object, newdata = newdata, n.ahead = n_ahead)
      
      for (k in 1: B) {
        tmp <- tvLM(y_b[k, ] ~ cumulative_x_reset, 
                    data = data_insample,# z = (1:nobs)/nobs,
                    est = "ll") 
        #prediction1[k, ] <- predict(tmp, newdata = newdata, newz = seq(nobs + 1, nobs + n_ahead)/(nobs + n_ahead))
        prediction1[k, ] <- forecast(tmp, newdata = newdata, n.ahead = n_ahead)
        #prediction_daily[k, ] <- tail(prediction1[k, ], -1) -  head(prediction1[k, ], -1)
      }# end of k loop  
      
      
      sd.est <- apply(prediction1, 2, sd)
      Q <- (prediction1 - matrix(pred_raw, nrow=B, ncol=length(pred_raw), byrow=TRUE))/
        matrix(sd.est, nrow=B, ncol=length(sd.est), byrow=TRUE)
      calpha <- apply(Q, 2, function(x){quantile(x, 1-alpha/2, na.rm = TRUE)})
      
      # output
      ans <- list (lower = pred_raw - sd.est * calpha, upper = pred_raw + sd.est * calpha)
    }
    
    # Simultaneous confidence bands
    joint_conf <- function (data, n_ahead, B = 200, alpha = 0.05) { # df after the best lead being selected
      #B = 1000 ## number of bootstrap draws
      data_insample <- head (data, dim(data)[1] - n_ahead)
      nobs <- dim(data_insample)[1]
      object <- tvLM(y_best_lead_tv ~ cumulative_x_reset, 
                     data = data_insample,
                     est = method) 
      
      residuals_raw <- object$residuals - mean(object$residuals)
      residuals_b <- matrix(sample (residuals_raw * rnorm(2 * B * nobs), size = 2 * B * nobs, replace = TRUE),
                            nrow = 2 * B, ncol = nobs)
      y_b <- matrix(object$fitted, nrow = 2 * B, ncol = nobs, byrow = TRUE) + residuals_b ## synthetic y
      
      ## out-of-sample prediction 1 with death counts observable
      prediction1 <- matrix(NA, nrow = 2 * B, ncol = n_ahead)
      totobs <- nobs + n_ahead
      
      newdf <- data %>% tail(n_ahead) %>% dplyr::select (c("cumulative_x_reset", "y_best_lead_tv"))
      newdata <- cbind(rep(1, n_ahead), as.matrix(newdf[, 1]))
      pred_raw <- forecast (object, newdata = newdata, n.ahead = n_ahead)
      
      for (k in 1: (2 * B)) {
        tmp <- tvLM(y_b[k, ] ~ cumulative_x_reset, 
                    data = data_insample, #z = (1:nobs)/nobs,
                    est = "ll") 
        #prediction1[k, ] <- predict(tmp, newdata = newdata, newz = seq(nobs + 1, nobs + n_ahead)/(nobs + n_ahead))
        prediction1[k, ] <- forecast(tmp, newdata = newdata, n.ahead = n_ahead)
      }# end of k loop  
      
      sd.est <- apply(prediction1[1:B,], 2, sd)
      
      # estimate common c_alpha 
      prediction2 <- prediction1[-(1:B), ]
      Q <- abs(prediction2 - 
                 matrix(pred_raw, nrow=B, ncol=length(pred_raw), byrow=TRUE))/
        matrix(sd.est, nrow=B, ncol=length(sd.est), byrow=TRUE)
      Qstar <- apply(Q, 2, max)
      calpha <- quantile(Qstar, 1 - alpha/2, na.rm = TRUE)  
      # output
      ans <- list (lower = pred_raw - sd.est * calpha, upper = pred_raw + sd.est * calpha)
    }
    
    #Add cumulative prediction bands for TVReg
    tv_conf_point <- point_conf(df_tv %>% filter (Date < D2), n_ahead_tv)
    opt_tv$plower <- c(rep(NA, length(tvLM_insample$fitted)), tv_conf_point$lower)
    opt_tv$pupper <-  c(rep(NA, length(tvLM_insample$fitted)), tv_conf_point$upper)
    
    tv_conf_joint <- joint_conf(df_tv %>% filter (Date < D2), n_ahead_tv)
    opt_tv$jlower <- c(rep(NA, length(tvLM_insample$fitted)), tv_conf_joint$lower)
    opt_tv$jupper <-  c(rep(NA, length(tvLM_insample$fitted)), tv_conf_joint$upper)
    
    # Plots
    names(opt_tv) <- c(
      "Predicted cumulative counts",
      "Reported cumulative counts",
      "Predicted daily counts",
      "date",
      "Reported daily counts",
      "Lower Confidence Bound (Pointwise)",            
      "Upper Confidence Bound (Pointwise)",
      "Lower Confidence Bound (Joint)",
      "Upper Confidence Bound (Joint)"
      )
    
    gg_cumulative <- opt_tv %>%
      pivot_longer(cols = where(is.numeric)) %>%
      filter(name %in% c("Lower Confidence Bound (Joint)",
                         "Lower Confidence Bound (Pointwise)", 
                         "Predicted cumulative counts",
                         "Reported cumulative counts",
                         "Upper Confidence Bound (Joint)",
                         "Upper Confidence Bound (Pointwise)"
                        )) %>%
      ggplot(aes(date, value, color = name, linetype = name)) +
      #ylim(range(c(pred, obs), na.rm = TRUE)) +
      geom_line(size = 0.7) +
      geom_vline(xintercept = D2, color = "orange") + 
      scale_linetype_manual(values=c("dotted", "dashed", "dotdash", "solid", "dotted", "dashed")) +
      scale_color_manual(values=c('purple', "blue", "red", "black", "purple", "blue")) +
      labs(
        title = paste0("Reported vs. Predicted Cumulative Counts ", region),
        subtitle = "Inputs = Cumulative Cases",
        x = "Date",
        y = "Count",
        caption = paste(method, "method lag=", n_ahead_tv, sep = " ")
      )
    
     } else {
       
       names(opt_tv) <- c(
         "Predicted cumulative counts", 
         "Reported cumulative counts",
         "Predicted daily counts",
         "date",
         "Reported daily counts"
       )
       
       gg_cumulative <- opt_tv %>%
         pivot_longer(cols = where(is.numeric)) %>%
         filter(name %in% c("Predicted cumulative counts",
                            "Reported cumulative counts")) %>%
         ggplot(aes(date, value, color = name, linetype = name)) +
         #ylim(range(c(pred, obs), na.rm = TRUE)) +
         geom_line(size = 0.7) +
         geom_vline(xintercept = D2, color = "orange") + 
         scale_linetype_manual(values=c("dotted", "solid")) +
         scale_color_manual(values=c("red", "black")) +
         labs(
           title = paste0("Reported vs. Predicted Cumulative Counts ", region),
           subtitle = "Inputs = Cumulative Cases",
           x = "Date",
           y = "Count",
           caption = paste(method, "method lag=", n_ahead_tv, sep = " ")
         )
       
     } # else for CI
    
    gg_daily <- opt_tv %>%
      pivot_longer(cols = where(is.numeric)) %>%
      filter(name %in% c("Predicted daily counts",
                         "Reported daily counts")) %>%
      ggplot(aes(date, value, color = name, linetype = name, size = name)) +
      geom_point() + geom_line() + 
      geom_vline(xintercept = D2, color = "orange") +
      scale_linetype_manual(values=c("dashed", "solid")) +
      scale_color_manual(values=c('Red', 'Black')) +
      scale_size_manual(values = c(1, 0.7)) + 
      labs(
        title = paste("Reported vs. Predicted Daily Counts,", region, sep = " "),
        subtitle = "Inputs = Cumulative Cases",
        x = "Date",
        y = "Count",
        caption = paste(method, "method lag =", n_ahead_tv, "days", sep = " ")
      )
    

 return(list(daily = gg_daily, cumulative = gg_cumulative, mspe = mspe, out = opt_tv, lag = n_ahead_tv, D2 = D2))

    }

