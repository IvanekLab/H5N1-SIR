# Baseline Parameters

N <- 4014 # Herd size incl. 127 heifers
W <- 0    # Number of workers
I_Ds0 <- 3  # Initial subclinically infected dry cows at t=0 (plausible 1-42)
lac <- 0.858  # Proportion of cows lactating (=1-(443 dry +127 heifers)/4014 total)
v <- 0  # Proportion immune workers at start

# Initial state vector (named)
initial_state_values <- c(
  
  # Lactating
  S_L = (lac * N), # Susceptible lactating
  I_Lc = 0, # Infected lactating clinical
  I_Ls = 0, # Infected lactating subclinical
  H_L = 0, # Hospital lactating
  R_L = 0, # Recovered lactating
  R_Lc = 0, # Recovered lactating clinical
  R_Ls = 0, # Recovered lactating subclinical
  
  # Dry
  S_D = ((1 - lac) * N) - I_Ds0, # Susceptible dry
  I_Dc = 0, # Infected dry clinical
  I_Ds = I_Ds0, # Infected dry subclinical
  H_D = 0, # Hospital dry
  R_D = 0, # Recovered dry
  R_Dc = 0, # Recovered dry clinical
  R_Ds = 0, # Recovered dry subclinical
  
  # Workers | TODO: Deprecate
  S_W = (1 - v) * W,
  I_Wc = 0,
  I_Ws = 0,
  R_W = v * W,
  
  # Deaths and culls
  Dd = 0, # Cumulative deaths due to infection
  Q = 0, # Cumulative culls due to infection
  
  # Infections
  Lic = 0, # Cumulative lactating clinical
  Lis = 0, # Cumulative lactating subclinical
  Dic = 0, # Cumulative dry clinical
  Dis = I_Ds0, # Cumulative dry subclinical
  
  # Workers | TODO: Deprecate
  Wic = 0, # Cumulative clinical workers
  Wr = 0, # Cumulative recovered workers
  
  # Force of infection
  Fi_L=0, # Lactating
  Fi_D=0, # Dry
  Fi_W=0 # Workers
)

## Fixed parameters
mu1 <- 0.056 / 365 # Natural mortality rate (NAHMS Dairy 2014)
mu2 <- 0.284 / 365 # Replacement rate excl. mortality (NAHMS Dairy 2014)
B <- N * (mu1 + mu2) # Susceptible cattle introduced during outbreak
beta_c <- 0.000563 # Transmission rate from clinical cow (BO)
beta_s <- 0.000057 # Transmission rate from subclinical cow (BO)
p <- 0.18 # Probability of symptomaticity if clinical (BO)
Dc <- 7.9 # Duration of infectiousness from clinical cow (BO)
Ds <-Dc # Duration of infectiousness from subclinical cow (BO)
Dh <-5.1 # Duration of hospital stay
alpha_c <- 1 / Dc # Recovery rate of clinical cow
alpha_s <- 1/ Ds # Recovery rate of subclinical cow
alpha_H <- 1 / 5.1 # Rate of leaving hospital
ro <- 311/777 # Probability of clinically infected cow being isolated
delta_1 <-0.068 /Dc # Mortality rate of clinical cows
delta_2 <- 0.316 /Dc # Culling rate of clinical cows
Dr <- 5*365 # Duration of protective immunity post infection (~0 to prevent recovery during acute outbreak)
gamma <- 1 / Dr # Rate of immunity loss
d <- 1 / 305 # Rate of lactating cows drying off
f <- 1 / 60 # Rate of dry cows freshening
alpha_W <- 1 / 7 # Recovery rate for workers
pi <- 0.5 # Probability of worker being clinical
z <- 0 # Worker turnover rate
B_W <- z * W # New workers introduced

# Time horizon
horizon <- 118 # Date of followup
time <- seq(1, horizon, by = 1) # Daily timestep

# Parameters
parameters <- c(
  mu1 = mu1,
  mu2 = mu2,
  B = B,
  kappa_1 = 1,  # Infectious lactating cow transmission rate modifier (BO) 
  kappa_2 = 1,  # Infectious dry cow transmission rate modifier (BO)
  kappa_3 = 1,   # Infectious separate barn transmission rate modifier (BO)
  kappa_4 = 1,       # Cow-to-human transmission rate modifier (BO)
  beta_c = beta_c,
  beta_s = beta_s,
  p = p,
  alpha_c = alpha_c,
  alpha_s = alpha_s,
  alpha_H = alpha_H,
  ro = ro,
  delta_1 = delta_1,
  delta_2 = delta_2,
  gamma = gamma,
  d = d,
  f = f,
  v = v,
  alpha_W = alpha_W,
  pi = pi,
  z = z,
  B_W = B_W
)