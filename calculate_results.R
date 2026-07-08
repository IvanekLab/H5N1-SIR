# Calculate results

output$herd <- (output$S_L) + (output$I_Lc) + (output$I_Ls) + (output$H_L) + (output$R_L) +
  (output$S_D) + (output$I_Dc) + (output$I_Ds) + (output$H_D) + (output$R_D)
herd_118 <-tail(output$herd, 1)

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
M_peak_d

#PREDICT P2
# Calculate peak number of predicted clinical (lactating and dry) infected cows
M_peak_n <- max(output$daily_new_clinical[-1])
M_peak_n

#PREDICT P3a over time
# Calculate total number ever infected cattle over time 
output$M_total_inf_n <- (output$I_Lc) + (output$I_Ls) + (output$H_L) + (output$R_L) +
  (output$I_Dc) + (output$I_Ds) + (output$H_D) + (output$R_D)
output$M_total_inf_n

#PREDICT P3b
# Calculate seroprevalence at the end of simulation 
M_total_inf_p  <- tail(output$M_total_inf_n, 1) / herd_118 
M_total_inf_p



#PREDICT P5a
# Calculate total ever clinically infected cattle over time 
output$M_clin_n  <- (output$I_Lc) + (output$H_L) + (output$R_Lc) + 
  (output$I_Dc) + (output$H_D) + (output$R_Dc)
output$M_clin_n

#PREDICT P5b
# Calculate seroprevalence of ever clinical at the end of simulation 
M_clin_p  <- tail(output$M_clin_n,1) / herd_118 
M_clin_p


#PREDICT P6a
# Calculate cumulative total subclinically infected cattle by end of simulation, correcting for survival bias 
output$M_subclin_n  <- (output$I_Ls) + (output$R_Ls) + (output$I_Ds) + (output$R_Ds)
output$M_subclin_n

#PREDICT P6b
#seroprev <- total_infected / N
M_subclin_p  <- tail(output$M_subclin_n,1) / herd_118 
M_subclin_p

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