# Bayesian optimization

# Setup
objective_function <- function(kappa_1, kappa_2, kappa_3, I_Ds0) {
  # Update parameters locally
  params <- parameters
  params["kappa_1"] <- kappa_1 
  params["kappa_2"] <- kappa_2
  params["kappa_3"] <- kappa_3
  params["I_Ds0"]   <- I_Ds0 
  
  # Run model safely
  output <- tryCatch({
    as.data.frame(lsoda(
      y = initial_state_values,
      func = HPAI_dyn,
      parms = params,
      times = time
    ))
  }, error = function(e) {
    message("Model run error: ", e$message)
    return(NULL)
  })
  
  if (is.null(output)) return(list(Score = -1e6))
  
  ##### --------- Model predictions 
  # Check herd size at day 118 exists and is valid
  output$herd <- (output$S_L) + (output$I_Lc) + (output$I_Ls) + (output$H_L) + (output$R_L) +
    (output$S_D) + (output$I_Dc) + (output$I_Ds) + (output$H_D) + (output$R_D)
  #alternatively
  #output$herd <- (output$S_L) + (output$I_Lc) + (output$I_Ls) + (output$H_L) + (output$R_Lc) + (output$R_Ls) +
  #  (output$S_D) + (output$I_Dc) + (output$I_Ds) + (output$H_D) + (output$R_Dc) + (output$R_Ds)
  herd_118 <-tail(output$herd, 1)
  
  if (length(herd_118) == 0 || is.na(herd_118)) {
    print("DEBUG: herd_118 is missing or NA")
    return(list(Score = -1e6))  # Penalty for missing data
  }
  
  # 
  # Process predictions to calculate pick day and peak number of cases
  output$daily_new_clinical_lact <- c(NA, diff(output$Lic))
  output$daily_new_subclinical_lact <- c(NA, diff(output$Lis))
  output$daily_new_clinical_dry <- c(NA, diff(output$Dic))
  output$daily_new_subclinical_dry <- c(NA, diff(output$Dis))
  
  #PREDICT P4
  #predict total daily clinical cases
  output$daily_new_clinical <- output$daily_new_clinical_lact + output$daily_new_clinical_dry
  #predict total clinical and subclinical cases
  output$daily_new_total <- output$daily_new_clinical_lact + output$daily_new_subclinical_lact +
    output$daily_new_clinical_dry + output$daily_new_subclinical_dry
  #output$daily_new_total[is.na(output$daily_new_total)] <- 0 ## replace NAs with zero to prevent error in ploting
  
  #PREDICT P1
  # Calculate peak day of predicted clinical (lactating and dry) infected cows
  M_peak_d <- output$time[which.max(output$daily_new_clinical)]
  
  #PREDICT P2
  # Calculate peak number of predicted clinical (lactating and dry) infected cows
  M_peak_n <- max(output$daily_new_clinical[-1])
  
  #PREDICT P3a
  # Calculate total number of ever infected cattle at t=118 
  M_total_inf_n <- tail(output$I_Lc, 1) + tail(output$I_Ls, 1) + tail(output$H_L, 1) +  tail(output$R_L, 1) + 
    tail(output$I_Dc, 1) +  tail(output$I_Ds, 1) +  tail(output$H_D, 1) +  tail(output$R_D, 1)
  
  #PREDICT P3b
  # Calculate seroprevalence of ever infected cattle at t=118 
  M_total_inf_p  <- M_total_inf_n / herd_118 
  
  
  #PREDICT P5a
  # Calculate total number of ever clinically infected cattle at t=118 
  M_clin_n  <- tail(output$I_Lc, 1) + tail(output$H_L, 1) + tail(output$R_Lc, 1) + 
    tail(output$I_Dc, 1) + tail(output$H_D, 1) + tail(output$R_Dc, 1)
  
  #PREDICT P5b
  # Calculate seroprevalence of ever clinically infected cattle at t=118 
  M_clin_p  <- M_clin_n / herd_118 
  
  
  #PREDICT P6a
  # Calculate total number of ever subclinically infected cattle at t=118 
  M_subclin_n  <-  tail(output$I_Ls, 1) + tail(output$R_Ls, 1) + tail(output$I_Ds, 1) + tail(output$R_Ds, 1)
  
  #PREDICT P6b
  # Calculate seroprevalence of ever infected cattle at t=118 
  M_subclin_p  <- M_subclin_n / herd_118 
  
  
  #PREDICT P7
  M_total_inf_lac_n<-	
    tail(output$I_Lc, 1) + 
    tail(output$I_Ls, 1) + 
    tail(output$H_L, 1) + 
    tail(output$R_L, 1)  
  M_total_inf_lac_n
  M_total_inf_lac_p<-M_total_inf_lac_n / (M_total_inf_lac_n + tail(output$S_L, 1) )
  M_total_inf_lac_p
  
  #PREDICT P8
  M_total_inf_dry_n<-	
    tail(output$I_Dc, 1) + 
    tail(output$I_Ds, 1) + 
    tail(output$H_D, 1) + 
    tail(output$R_D, 1)  
  M_total_inf_dry_n
  M_total_inf_dry_p<-M_total_inf_dry_n / (M_total_inf_dry_n + tail(output$S_D, 1) )
  M_total_inf_dry_p
  
  
  
  # NA checks on model outputs
  if (any(is.na(c(M_peak_d, M_peak_n, M_total_inf_p, M_clin_p, M_subclin_p, M_total_inf_lac_p, M_total_inf_dry_p)))) {
    print("DEBUG: one of model predictions P2, P2, P3b, P5b, P6b, P7, P8 are NA") 
    return(list(Score = -1e6))
  }
  
  
  ##### --------- Dynamic and fixed model targets 
  # TARGET T1
  T_peak_d <- 23
  
  # TARGET T2
  T_peak_n <- 121
  
  
  
  # TARGET T3a
  T_total_inf_n <- p.sero_pos * herd_118 #0.89*herd(t=118)
  # TARGET T3b
  T_total_inf_p  <- p.sero_pos #0.89
  
  
  
  ## Dynamic targets based on model herd size at day 118
  #TARGET 5a
  T_clin_n <- p.sero_pos_clin * p.sero_pos * herd_118 #0.15*0.89*herd(t=118)
  #TARGET 5b
  T_clin_p <- p.sero_pos_clin * p.sero_pos #0.133=0.15*0.89
  
  
  
  #TARGET 6a
  T_subclin_n  <- p.sero_pos_subclin * p.sero_pos *  herd_118 #0.85*0.89*herd(t=118)
  #TARGET 6b
  T_subclin_p  <- p.sero_pos_subclin * p.sero_pos # 0.761=485/570
  
  
  
  #TARGET 7 
  #Prevalence of total ever infected  lactating at t=118
  #target=0.929
  V_total_inf_lac_p <- 0.929
  
  #TARGET 8 
  #Prevalence of total ever infected  dry at t=118
  #target=0.405
  V_total_inf_dry_p <-0.405
  
  
  #TARGET T4
  output$observed_cases <- outbreak_data$cases[match(output$time, outbreak_data$day)]
  output$observed_cases[is.na(output$observed_cases)] <- 0 ## replace NAs with zero to prevent error in plotting
  
  E1 <- ((M_peak_d - T_peak_d) / T_peak_d)^2 #T1
  E2 <- ((M_peak_n - T_peak_n) / T_peak_n)^2 #T2
  E3b <- ((M_total_inf_p  - T_total_inf_p) / T_total_inf_p )^2 #T3b
  E5b <- ((M_clin_p - T_clin_p) / T_clin_p)^2 #T5b
  E6b <- ((M_subclin_p  - T_subclin_p ) / T_subclin_p)^2 #T6b
  E7 <- ((M_total_inf_lac_p - V_total_inf_lac_p) / V_total_inf_lac_p)^2 #T7
  E8 <- ((M_total_inf_dry_p  - V_total_inf_dry_p ) / V_total_inf_dry_p)^2 #T8
  
  
  # Daily clinical cases error
  # Assumes observed_cases vector is available and aligned with 'times'
  model_daily <- output$daily_new_clinical[-1]      # remove initial NA for alignment
  obs_daily <- output$observed_cases[-1]                    # observed daily cases vector
  
  
  # If lengths mismatch or NA, return penalty
  if (length(model_daily) != length(obs_daily) || any(is.na(obs_daily))) {
    print("DEBUG: Length mismatch or NA in daily cases")  # <---- Inserted here
    print(sprintf("model_daily length: %d, obs_daily length: %d", length(model_daily), length(obs_daily)))
    print(sprintf("NA in obs_daily? %s", any(is.na(obs_daily))))
    return(list(Score = -1e6))
  }
  
  #The Relative Mean Squared Error (RMSE) is a normalized version of the Mean Squared Error (MSE), 
  #used here because it is comparable to other errors. 
  #calculated by dividing the Root Mean Squared Error (RMSE) by mean of the observed values.
  E4_step1 <- mean(abs(model_daily - obs_daily)^2, na.rm = TRUE) #MSE
  E4_step2 <- sqrt(E4_step1) # root of MSE
  E4_step3 <- mean(obs_daily) #mean of observed values
  E4 <- E4_step2/E4_step3 #RMSE
  
  
  #prediction<-as.data.frame(t(c(M_peak_d, M_peak_n, M_total_inf_p, M_clin_p, M_subclin_p, M_total_inf_lac_p, M_total_inf_dry_p)))
  #colnames(prediction) <- c("M_peak_d", "M_peak_n", "M_total_inf_p", "M_clin_p", "M_subclin_p", "M_total_inf_lac_p", "M_total_inf_dry_p")
  prediction<-c(M_peak_d, M_peak_n, M_total_inf_p, M_clin_p, M_subclin_p, M_total_inf_lac_p, M_total_inf_dry_p)
  
  
  
  # Combine errors (weighted sum)
  total_error <- E1 +E2 + E3b + E5b + E6b + E7 + E8 #+ E4
  
  print(sprintf("kappa_1=%.3f, kappa_2=%.3f, kappa_3=%.3f, I_Ds0=%1.0f, M_peak_d=%1.0f, M_peak_n=%.3f,M_total_inf_p=%.3f, M_clin_p=%.3f, M_subclin_p=%.3f,M_total_inf_lac_p=%.3f, M_total_inf_dry_p=%.3f, total_error=%.4f",
                kappa_1, kappa_2, kappa_3, I_Ds0, M_peak_d, M_peak_n, M_total_inf_p, M_clin_p, M_subclin_p, M_total_inf_lac_p, M_total_inf_dry_p, total_error))
  
  
  
  #return(list(Score = -total_error))  # Negative because optimizer maximizes Score
  return(list(Score = -total_error, Pred = prediction))  # Negative because optimizer maximizes Score
}

# Parameters for Bayesian optimization
n.sero_N<- 637 # number animals serologically tetsed at day 118
n.sero_pos <- 570 #  #(637-67=570); number cows seropositive
n.sero_neg <-67 # (637-570=67) number cows seronegative
n.sero_pos_clin <- 85 #number seropositive with clinical history 
n.sero_pos_subclin <- 485 #number seropositive with NO clinical history (i.e., subclinical)

p.sero_pos <- n.sero_pos/n.sero_N #  0.895=570/637 seroprevalence in the herd at t=118 (9)also final size epidemic)
p.sero_neg <- 1-p.sero_pos # prevalence of seronegatives at t=118 (also proportion susceptible in the herd)
p.sero_pos_clin <- n.sero_pos_clin/n.sero_pos #0.149=85/570 proportion with history of clinical infection among seropositive at t=118
p.sero_pos_subclin <- 1-p.sero_pos_clin #0.851=485/570 proportion with history of subclinical infection among seropositive at t=118


# Run
opt_results <- BayesianOptimization(
  FUN = objective_function,
  bounds = list(kappa_1 = c(1, 50), kappa_2 = c(0.00001, 1), kappa_3 = c(0.0000000001, 1), I_Ds0 = c(1, 42)),
  
  init_points = 10,
  n_iter = 50,
  acq = "ucb",
  kappa = 2.576,
  verbose = TRUE
)

# Export results
best_params <- opt_results$Best_Par
print(best_params)

write.csv(best_params, "best_params.csv", row.names = TRUE)

# write.csv(opt_results$History,"optimization_history.csv", row.names =FALSE)

predictions<-as.data.frame(t(opt_results$Pred))
colnames(predictions) <- c("M_peak_d", "M_peak_n", "M_total_inf_p", "M_clin_p", "M_subclin_p", "M_total_inf_lac_p", "M_total_inf_dry_p")
predictions
# Export to CSV predictions from optimization history
write.csv(predictions,"optimization_predictions.csv", row.names =TRUE)