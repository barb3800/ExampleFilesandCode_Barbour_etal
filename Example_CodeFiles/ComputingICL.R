# Computing Integrated Completed Likelihood (ICL) from fitted model

require(magrittr)

viterbi <- extract(model1_example, pars = "viterbi")[[1]]

state.hat <- ceiling(colMeans(viterbi))
mu.hats <- rstan::extract(model1_example,pars=c("mu"))$mu %>% apply(2, median) 
sigma.hats <- rstan::extract(model1_example,pars=c("sigma"))$sigma %>% apply(2, median)
lambda.hat <- median(rstan::extract(model1_example,pars=c("lambda"))$lambda)
Pr.hats <- rstan::extract(model1_example,pars=c("Pr"))$Pr %>% apply(2, median)

logLikeByState <- with(data, 
            dcauchy(V, mu.hats[state.hat], sigma.hats[state.hat], log = TRUE) + 
              dpois(N1, lambda.hat*Pr.hats[state.hat], log = TRUE) + 
              dpois(N2, lambda.hat*(1 - Pr.hats[state.hat]), log = TRUE)) 
              
p <- 3 + #mus
  3 + #sigmas
  1 + #lambda
  3 + #Prs
  6 # transition probabilities

ICL <- -2*sum(logLikeByState) + p * log(data$T)

