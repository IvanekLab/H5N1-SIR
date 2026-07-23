# Main model - Based on Slide #1 in Slide Deck (v7b)

# Model code
HPAI_dyn <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    lambda_L <- kappa_1 * beta_c * I_Lc +
      kappa_1 * beta_s * I_Ls +  
      kappa_3 * kappa_2 * beta_c * I_Dc +  
      kappa_3 * kappa_2 * beta_s * I_Ds +  
      kappa_1 * kappa_3 * beta_c * H_L +  
      kappa_3 * kappa_2 * beta_c * H_D
    
    lambda_D <- kappa_2 * beta_c * I_Dc + 
      kappa_2 * beta_s * I_Ds + 
      kappa_3 * kappa_1 * beta_c * I_Lc + 
      kappa_3 * kappa_1 * beta_s * I_Ls +
      kappa_1 * kappa_3 * beta_c * H_L +  
      kappa_3 * kappa_2 * beta_c * H_D
    
    lambda_W <- kappa_1 * kappa_4 * beta_c * I_Lc + 
      kappa_1 * kappa_4 * beta_s * I_Ls + 
      kappa_4 * kappa_2 * beta_c * I_Dc + 
      kappa_4 * kappa_2 * beta_s * I_Ds
    kappa_1 * kappa_4 * beta_c * H_L +  
      kappa_4 * kappa_2 * beta_c * H_D
    
    dS_L <- B + gamma * R_L + f * S_D - lambda_L * S_L - (d + mu1 + mu2) * S_L
    dI_Lc <- f * I_Dc + p * (1-ro) * lambda_L * S_L + alpha_H * H_L - (d + mu1 + mu2 + alpha_c + delta_1 + delta_2) * I_Lc
    dI_Ls <- f * I_Ds + (1 - p) * lambda_L * S_L - (d + mu1 + mu2 + alpha_s) * I_Ls
    dH_L <- f * H_D + p * ro * lambda_L * S_L - (d + mu1 + mu2 + alpha_H + delta_1 + delta_2) * H_L
    dR_L <- f * R_D + alpha_s * I_Ls + alpha_c * I_Lc - (d + mu1 + mu2 + gamma) * R_L
    dR_Lc <- f * R_Dc + alpha_c * I_Lc - (d + mu1 + mu2 + gamma) * R_Lc
    dR_Ls <- f * R_Ds + alpha_s * I_Ls - (d + mu1 + mu2 + gamma) * R_Ls
    
    dS_D <- gamma * R_D + d * S_L - lambda_D * S_D - (f + mu1 + mu2) * S_D
    dI_Dc <- d * I_Lc + p * (1-ro) * lambda_D * S_D + alpha_H * H_D - (f + mu1 + mu2 + alpha_c + delta_1 + delta_2) * I_Dc
    dI_Ds <- d * I_Ls + (1 - p) * lambda_D * S_D - (f + mu1 + mu2 + alpha_s) * I_Ds
    dH_D <- d * H_L + ro * p * lambda_D * S_D - (f + mu1 + mu2 + alpha_H + delta_1 + delta_2) * H_D
    dR_D <- d * R_L + alpha_s * I_Ds + alpha_c * I_Dc - (f + mu1 + mu2 + gamma) * R_D
    dR_Dc <- d * R_Lc + alpha_c * I_Dc - (f + mu1 + mu2 + gamma) * R_Dc
    dR_Ds <- d * R_Ls + alpha_s * I_Ds - (f + mu1 + mu2 + gamma) * R_Ds
    
    dS_W <- (1 - v) * B_W - lambda_W * S_W - z * S_W
    dI_Wc <- pi * lambda_W * S_W - alpha_W * I_Wc - z * I_Wc
    dI_Ws <- (1 - pi) * lambda_W * S_W - alpha_W * I_Ws - z * I_Ws
    dR_W <- alpha_W * (I_Wc + I_Ws) + v * B_W - z * R_W
    
    dDd <- delta_1 * (I_Lc + I_Dc + H_L + H_D)
    dQ <- delta_2 * (I_Lc + I_Dc + H_L + H_D)
    
    dLic <- p * lambda_L * S_L
    dLis <- (1 - p) * lambda_L * S_L
    dDic <- p * lambda_D * S_D
    dDis <- (1 - p) * lambda_D * S_D
    
    dWic <- pi * lambda_W * S_W
    dWr <- alpha_W * (I_Wc + I_Ws)
    
    dFi_L <-lambda_L
    dFi_D <-lambda_D  
    dFi_W <-lambda_W
    
    return(list(c(
      dS_L, dI_Lc, dI_Ls, dH_L, dR_L, dR_Lc, dR_Ls, 
      dS_D, dI_Dc, dI_Ds, dH_D, dR_D, dR_Dc, dR_Ds,
      dS_W, dI_Wc, dI_Ws, dR_W,
      dDd, dQ, 
      
      dLic, dLis, dDic, dDis,
      dWic, dWr,
      
      dFi_L, dFi_D, dFi_W
    )))
  })
}

# Generate outputs
output <- tryCatch({
  as.data.frame(lsoda(
    y = initial_state_values,
    func = HPAI_dyn,
    parms = parameters,
    times = time
  ))
}, error = function(e) {
  message("Model run error: ", e$message)
  return(NULL)
})
