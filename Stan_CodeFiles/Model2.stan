// Model 1 Structure:
// 3 hidden behavior states
// HMM DCRW model in Stan (hidden markov model, first difference correlated random walk)- using to identify hidden/latent states of time series

// Variables:
// Swimming speed - Weibull distribution (parameters alpha and sigma)
// Theta or turning angle - Wrapped Cauchy distribution (parameter rho)
// Dive type- count of deep (D1/N1) or shallow (D2/N2) dives and total count of all dives - Poisson distribution (parameters Pr and lambda)

// Model Modifications:
// We create a "transformed parameter", Pr, which corresponds to the probability of a dive being type D1 (deep) (the probability of being shallow/D1 is simply 1-Pr)
// Pr for each behavior state (S1-S3) is multiplied by lambda (total count of dives per day) and then modeled with a Poisson process for N1 (Pr*lambda) and N2 (1-Pr*lambda)
// We also place an ordering on Pr, such that S1 Pr < S2 Pr < S3 Pr, to assist with convergence and due to biological expectation that turtles have 3 vertical behaviors (shallower, intermediate, deeper)
// We constrain the initial Pr (Pinit) in the parameters section to be between 0.8-1 (for S3)- the S1 Pr and S2 Pr are then derived from this in the transformed section (S1 Pr= 1- S3 Pr and S2 Pr= 1-[S1 Pr + S1 Pr])
// We place semi-vague priors for the alpha parameter, such that S1 alpha > S2 alpha > S3 alpha (prior knowledge has an expectation that turtles transiting have shallow dives and faster speeds)
// in the generated quantities block, we also derive the log likelihood of each data point (log_lik)
// which is needed for the computation of loo and model selection

// create function to apply wrapped Cauchy distribution to rho, since Stan doesn't have one
functions {
  real wrappedCauchy_lpdf(real aPhi, real aRho) {
      return(- log(2*pi()) + log(1-aRho^2) - log(1+aRho^2-2*aRho*cos(aPhi)));
  }
}

data{
  int<lower=1> T; // number of observations
  vector<lower=0>[T] S;  // swim speed
  vector<lower=-pi(), upper=pi()>[T] Theta;  // turning angles
  int<lower=0> N1[T] ; // count of deep dives
  int<lower=0> N2[T] ; // count of shallow dives
}

parameters {
  vector<lower=0, upper=1>[3] rho;  // rho_1, rho_2 (theta)
  
  real<lower=0> alpha[3]; // parameter for weibull dist (S)
  vector<lower=0>[3] sigma; // parameter for weibull dist (S)
  
  real<lower = 0.8,upper=1> Pinit; // constrain to be greater than 0.8
  
  real<lower=0> lambda; // parameter for mean no. of dives per day- could have 2 lambdas (1 for each state) if we expect diff no. per state
  
  simplex[3] P[3];  // transition probability matrix
}


transformed parameters{
  matrix[3, 3] P_matrix; // 3 x 3 matrix
  simplex[3] P_star; // stationary distribution
  
  vector<lower = 0, upper=1>[3] Pr; // define transformed parameter, Pr, to be constrained from 0-1 (it's a probability)

  Pr[1] = 1 - Pinit; // smaller value (0-0.2)
  Pr[2] = 1 - (Pr[1]+Pr[1]); // median value (0-0.4)
  Pr[3] = Pinit; // larger value (0.8-1.0)
  
  for(j in 1:3){ // Si(t+1) | Sj(t) - you are in state Si given that you were in state Sj
  for(i in 1:3){
    P_matrix[i,j]= P[i,j];
  }}
  P_star = to_vector((to_row_vector(rep_vector(1.0, 3))/(diag_matrix(rep_vector(1.0,3)) - P_matrix + rep_matrix(1, 3, 3)))) ; // 3 x 3 matrix filled with 1's
}



model{
  vector[3] log_P[3];
  vector[3] lp;
  vector[3] lp_p1;
  // PRIORS (vague, semi-vague)
  lambda ~ exponential(0.16); // vague
  
  alpha ~ exponential(0.25); //vague
  sigma[1] ~ normal(60,30); // high speed for state 1
  sigma[2] ~ normal(1,1); // low speed for state 2
  sigma[3] ~ normal(15,10); // median speed for state 3
    
  // State Vector
  // transpose the P and take natural log of entries
  for (i in 1:3)
    for (j in 1:3)
      log_P[i,j] = log(P[i, j]);
      
  // First Observation
  // Note: there is no theta (turning angle) for the first observation (it is the angle between two locations)
  for(k in 1:3){
    lp[k] = log(P_star[k]) + 
    weibull_lpdf(S[1] | alpha[k], sigma[k]) + // swim speed
    poisson_lpmf(N1[1] | Pr[k]*lambda) + // count of deeper dives
    poisson_lpmf(N2[1] | (1-Pr[k])*lambda); // count of shallower dives
  }
      
  // forward algorithm implementation
  for(t in 2:T){
    for (k in 1:3) // looping over states
      lp_p1[k] = log_sum_exp(log_P[k] + lp) + 
      weibull_lpdf(S[t] | alpha[k], sigma[k]) + // swim speed
      poisson_lpmf(N1[t] | Pr[k]*lambda) + // count of deeper dives
      poisson_lpmf(N2[t] | (1-Pr[k])*lambda) + // count of shallower dives
      wrappedCauchy_lpdf(Theta[t] | rho[k]);
      lp = lp_p1;
  }
  target += log_sum_exp(lp);
}

// In this block, we compute
// the log likelihood of each data point (log_lik)
// which is needed for the computation of loo later
 
generated quantities {
  real log_lik[3];
  
  // First Observation
  for(k in 1:3){
    log_lik[k] = weibull_lpdf(S[1] | alpha[k], sigma[k]) + // swim speed
    poisson_lpmf(N1[1] | Pr[k]*lambda) + // count of deeper dives
    poisson_lpmf(N2[1] | (1-Pr[k])*lambda); // count of shallower dives
  }
      
  // the rest
  for(t in 2:T){
    for (k in 1:3) // looping over states
      log_lik[k] = weibull_lpdf(S[t] | alpha[k], sigma[k]) + // swim speed
      poisson_lpmf(N1[t] | Pr[k]*lambda) + // count of deeper dives
      poisson_lpmf(N2[t] | (1-Pr[k])*lambda) + // count of shallower dives
      wrappedCauchy_lpdf(Theta[t] | rho[k]);
  }
}
//
  