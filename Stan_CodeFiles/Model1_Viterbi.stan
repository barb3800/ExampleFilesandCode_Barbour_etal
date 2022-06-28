// Model 1 Structure:
// 3 hidden behavior states
// HMM DCRW model in Stan (hidden markov model, first difference correlated random walk)- using to identify hidden/latent states of time series

// Variables:
// Move Persistence Velocity - Cauchy distribution (parameters mu and sigma)
// Dive type- count of deep (D1/N1) or shallow (D2/N2) dives and total count of all dives - Poisson distribution (parameters Pr and lambda)

// Model Modifications:
// We create a "transformed parameter", Pr, which corresponds to the probability of a dive being type D1 (deep) (the probability of being shallow/D1 is simply 1-Pr)
// Pr for each behavior state (S1-S3) is multiplied by lambda (total count of dives per day) and then modeled with a Poisson process for N1 (Pr*lambda) and N2 (1-Pr*lambda)
// We also place an ordering on Pr, such that S1 Pr < S2 Pr < S3 Pr, to assist with convergence and due to biological expectation that turtles have 3 vertical behaviors (shallower, intermediate, deeper)
// We constrain the initial Pr (Pinit) in the parameters section to be between 0.8-1 (for S3)- the S1 Pr and S2 Pr are then derived from this in the transformed section (S1 Pr= 1- S3 Pr and S2 Pr= 1-[S1 Pr + S1 Pr])
// We place semi-vague priors for the mu parameter, such that S1 mu > S2 mu > S3 mu (prior knowledge has an expectation that turtles transiting have shallow dives and high persistence velocity)
// in the generated quantities block, we apply the Viterbi algorithm to predict the most likely behavioral state (S1-S3) for each data point in our input time series

data{
  int<lower=1> T; // number of observations
  vector[T] V;  // move persistence velocity
  int<lower=0> N1[T] ; // count of deep dives
  int<lower=0> N2[T] ; // count of shallow dives
}

parameters {
  vector<lower=0>[3] mu; // parameter for cauchy dist (move persistence)
  vector<lower=0>[3] sigma; // parameter for cauchy dist (move persistence)
  
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
  
  mu[1] ~ normal(60,30); // higher persistence v for state 1
  mu[2] ~ normal(1,1); // low persistence v for state 2
  mu[3] ~ normal(15,10); // median persistence v for state 3
  sigma ~ exponential(0.25); // vague
    
  // State Vector
  // transpose the P and take natural log of entries
  for (i in 1:3)
    for (j in 1:3)
      log_P[i,j] = log(P[i, j]);
      
  // First Observation
  for(k in 1:3){
    lp[k] = log(P_star[k]) + 
    cauchy_lpdf(V[1] | mu[k], sigma[k]) + //  move persistence
    poisson_lpmf(N1[1] | Pr[k]*lambda) + // count of deeper dives
    poisson_lpmf(N2[1] | (1-Pr[k])*lambda); // count of shallower dives
  }
      
  // forward algorithm implementation
  for(t in 2:T){
    for (k in 1:3) // looping over states
      lp_p1[k] = log_sum_exp(log_P[k] + lp) + 
      cauchy_lpdf(V[t] | mu[k], sigma[k]) + //  move persistence
      poisson_lpmf(N1[t] | Pr[k]*lambda) + // count of deeper dives
      poisson_lpmf(N2[t] | (1-Pr[k])*lambda); // count of shallower dives
      lp = lp_p1;
  }
  target += log_sum_exp(lp);
}

// VITERBI STATE SEQUENCE DETERMINATION
generated quantities {
  
  int<lower=1,upper=3> viterbi[T]; 
  vector[3] log_P[3];
  
  // Viterbi algorithm (most likely state sequence)  
  {
    real max_logp;
    int back_ptr[T, 3];
    real best_logp[T, 3];
    for (i in 1:3)
    for (j in 1:3)
      log_P[i,j] = log(P[i, j]);
      
    for (t in 1:T) {
      if(t==1) {
        for(k in 1:3) 
          best_logp[t,k]= cauchy_lpdf(V[t] | mu[k], sigma[k]) + 
                  poisson_lpmf(N1[t] | Pr[k]*lambda) + 
                  poisson_lpmf(N2[t] | (1-Pr[k])*lambda);
        } else {
        for (k in 1:3) {   // from j to k at step t
          best_logp[t, k] = negative_infinity();
          for (j in 1:3) {
            real logp;
            logp = best_logp[t-1, j] + log_P[j,k] + 
                cauchy_lpdf(V[t] | mu[k], sigma[k]) + 
                poisson_lpmf(N1[t] | Pr[k]*lambda) + 
                poisson_lpmf(N2[t] | (1-Pr[k])*lambda);
            
          if (logp > best_logp[t, k]) {
                back_ptr[t, k] = j;
                best_logp[t, k] = logp;
              }}}}
    }
        
     for(t0 in 1:T) {
        int t = T - t0 + 1;
        if(t==T) {
          max_logp = max(best_logp[t]);
          for (k in 1:3)
            if (best_logp[t, k] == max_logp)
              viterbi[t] = k;
        } else {
            viterbi[t] = back_ptr[t+1, viterbi[t+1]];}
       }
  }
}
//
  