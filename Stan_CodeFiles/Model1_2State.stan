// Model 2 Structure:
// 2 hidden behavior states
// HMM DCRW model in Stan (hidden markov model, first difference correlated random walk)- using to identify hidden/latent states of time series

// Variables:
// Move Persistence Velocity - Cauchy distribution (parameters mu and sigma)
// Dive type- count of deep (D1/N1) or shallow (D2/N2) dives and total count of all dives - Poisson distribution (parameter lambda)

// Model Modifications:
// We create a "transformed parameter", Pr, which corresponds to the probability of a dive being type D1 (deep) (the probability of being shallow/D1 is simply 1-Pr)
// P for each behavior state (S1-S2) is multiplied by lambda (total count of dives per day) and then modeled with a Poisson process for N1 (Pr*lambda) and N2 (1-Pr*lambda)
// We also place an ordering on Pr, such that S1 Pr < S2 Pr, to assist with convergence and due to biological expectation that turtles have at least 2 vertical behaviors (shallower, deeper)
// We constrain the initial Pr in the parameters section to be between 0.8-1 (for S2)- the S1 Pr is then derived from this in the transformed section (S1 Pr= 1- S3 Pr)
// We place semi-vague priors for the mu parameter, such that S1 mu > S2 mu (prior knowledge has an expectation that turtles transiting have shallow dives and high persistence velocity)
// in the generated quantities block, we also derive the log likelihood of each data point (log_lik)
// which is needed for the computation of loo and model selection

data{
  int<lower=1> T; // number of observations
  vector[T] V;  // move persistence velocity
  int<lower=0> N1[T] ; // count of deep dives
  int<lower=0> N2[T] ; // count of shallow dives
}

parameters {
  vector<lower=0>[2] mu; // parameter for cauchy dist (move persistence)
  vector<lower=0>[2] sigma; // parameter for cauchy dist (move persistence)
  
  real<lower = 0.8,upper=1> Pinit; // constrain to be greater than 0.8
  
  real<lower=0> lambda; // parameter for mean no. of dives per day- could have 2 lambdas (1 for each state) if we expect diff no. per state
  
  simplex[2] P[2];  // transition probability matrix
}


transformed parameters{
  matrix[2, 2] P_matrix; // 2 x 2 matrix
  simplex[2] P_star; // stationary distribution
  
  vector<lower = 0, upper=1>[2] Pr; // define transformed parameter, Pr, to be constrained from 0-1 (it's a probability)

  Pr[1] = 1 - Pinit; // smaller value (0-0.2)
  Pr[2] = Pinit; // larger value (0.8-1.0)
  
  for(j in 1:2){ // Si(t+1) | Sj(t) - you are in state Si given that you were in state Sj
  for(i in 1:2){
    P_matrix[i,j]= P[i,j];
  }}
  P_star = to_vector((to_row_vector(rep_vector(1.0, 2))/(diag_matrix(rep_vector(1.0,2)) - P_matrix + rep_matrix(1, 2, 2)))) ; // 2 x 2 matrix filled with 1's
}



model{
  vector[2] log_P[2];
  vector[2] lp;
  vector[2] lp_p1;
  // PRIORS (vague, semi-vague)
  lambda ~ exponential(0.16); // vague
  
  mu[1] ~ normal(60,30); // higher persistence v for state 1
  mu[2] ~ normal(1,1); // low persistence v for state 2
  sigma ~ exponential(0.25); // vague
    
  // State Vector
  // transpose the P and take natural log of entries
  for (i in 1:2)
    for (j in 1:2)
      log_P[i,j] = log(P[i, j]);
      
  // First Observation
  for(k in 1:2){
    lp[k] = log(P_star[k]) + 
    cauchy_lpdf(V[1] | mu[k], sigma[k]) + //  move persistence
    poisson_lpmf(N1[1] | Pr[k]*lambda) + // count of deeper dives
    poisson_lpmf(N2[1] | (1-Pr[k])*lambda); // count of shallower dives
  }
      
  // forward algorithm implementation
  for(t in 2:T){
    for (k in 1:2) // looping over states
      lp_p1[k] = log_sum_exp(log_P[k] + lp) + 
      cauchy_lpdf(V[t] | mu[k], sigma[k]) + //  move persistence
      poisson_lpmf(N1[t] | Pr[k]*lambda) + // count of deeper dives
      poisson_lpmf(N2[t] | (1-Pr[k])*lambda); // count of shallower dives
      lp = lp_p1;
  }
  target += log_sum_exp(lp);
}

// In this block, we compute
// the log likelihood of each data point (log_lik)
// which is needed for the computation of loo later
 
generated quantities {
  real log_lik[2];
  
  // First Observation
  for(k in 1:2){
    log_lik[k] = cauchy_lpdf(V[1] | mu[k], sigma[k]) + //  move persistence
    poisson_lpmf(N1[1] | Pr[k]*lambda) + // count of deeper dives
    poisson_lpmf(N2[1] | (1-Pr[k])*lambda); // count of shallower dives
  }
      
  // the rest
  for(t in 2:T){
    for (k in 1:2) // looping over states
      log_lik[k] = cauchy_lpdf(V[t] | mu[k], sigma[k]) + //  move persistence
      poisson_lpmf(N1[t] | Pr[k]*lambda) + // count of deeper dives
      poisson_lpmf(N2[t] | (1-Pr[k])*lambda); // count of shallower dives
  }
}
//