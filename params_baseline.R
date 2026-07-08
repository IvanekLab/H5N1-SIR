# Baseline Parameters

N <- 4014 
W <- 30    # Number of workers
I_Ds0 <- 3  # Initial subclinicaly infected dry cows at t=0 (plausible 1-42)
lac <- 0.858  
v <- 0  # Proportion immune workers at start

# Initial state vector (named)
initial_state_values <- c(
  S_L = (lac * N), 
  I_Lc = 0,
  I_Ls = 0,
  H_L = 0,
  R_L = 0,
  R_Lc = 0,
  R_Ls = 0,
  
  S_D = ((1 - lac) * N) - I_Ds0,
  I_Dc = 0,
  I_Ds = I_Ds0,
  H_D = 0,
  R_D = 0,
  R_Dc = 0,
  R_Ds = 0,
  
  S_W = (1 - v) * W,
  I_Wc = 0,
  I_Ws = 0,
  R_W = v * W,
  
  Dd = 0,
  Q = 0,
  
  Lic = 0,
  Lis = 0, 
  Dic = 0,
  Dis = I_Ds0,
  
  Wic = 0,
  Wr = 0,
  
  Fi_L=0,
  Fi_D=0,
  Fi_W=0
)

## Fixed parameters
mu1 <- 0.05 / 365
mu2 <- 0.32 / 365
B <- N * (mu1 + mu2)
beta_c <- 0.000563
beta_s <- 0.000057
p <- 0.18
Dc <- 7.9
Ds <-Dc
Dh <-5.1
alpha_c <- 1 / Dc
alpha_s <- 1/ Ds
alpha_H <- 1 / 5.1
ro <- 311/777
sigma <- 1/2 #0=animals are not isolated; 1/2=animals are isolated 1 days after clinical presentation
delta_1 <-0.068 /Dc
delta_2 <- 0.316 /Dc
Dr <- 5*365
gamma <- 1 / Dr ### set to zero to prevent loss of immunity during simulation
d <- 1 / 305
f <- 1 / 60
alpha_W <- 1 / 7
pi <- 0.5
z <- 0
B_W <- z * W

# Time horizon
horizon <- 118
time <- seq(1, horizon, by = 1)

# Parameters
parameters <- c(
  mu1 = mu1,
  mu2 = mu2,
  B = B,
  kappa_1 = 1,  # initial guess, will be optimized 
  kappa_2 = 1,  # initial guess, will be optimized
  kappa_3 = 1,   # initial guess, will be optimized
  kappa_4 = 1,       # initial guess, will be optimized
  beta_c = beta_c,
  beta_s = beta_s,
  p = p,
  alpha_c = alpha_c,
  alpha_s = alpha_s,
  alpha_H = alpha_H,
  ro = ro,
  sigma = sigma,
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