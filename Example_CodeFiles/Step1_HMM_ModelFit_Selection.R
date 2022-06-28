## EXAMPLE: FITTING HMM MODELS TO DATA & SELECTING "BEST" MODELS

## Author and Contact Info: Nicole Barbour, nbarbour@umd.edu 
## Code as example for Barbour et al., "Incorporating multidimensional behavior into a dynamic management tool for a critically endangered and migratory species"

##############################################
# SET UP 
##############################################

# load libraries (note: all these packages need to be installed before use)
# you will also need up-to-date versions of R and Rtools installed
## for help with rstan installation, see: https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
library(rstan);library(ggplot2)

# set working directory
setwd("yourfilepathhere")

# read in data (.rda format)
load("./Example_Data/Jitter_Turtle_Data.rda")
## data format:
## NOTE 1: this data is for 3 individual turtles that have had their datetime and lat/long info jittered (e.g., they are not real data...)

## NOTE 2: this data has been pre-processed (see Barbour et al.) to:
## remove erroneous positions 
## have daily/regular positions ("foieGras" package, Jonsen et al. 2019)
## correct observed ("ground") daily speeds for current-influence to determine "true" daily swimming speed (Gaspar et al. 2006)
## have daily turning angles, theta ("bcpa" R package, Gurarie 2014)
## have daily move persistence values, derived from swim speed * cos (turn angle) (Gurarie et al. 2009)
## cluster dives into groups ("D1"- deep or "D2"-shallow) using Dynamic Time Warping (Barbour et al.), with the count of each (D1 or D2) being found for each day/individual combination

#############################################
# VARIABLE EXPLORATION
#############################################

# in this example, we will be exploring both vertical & horizontal movement variables:
## variables include:
## Ref (id name)
## N1 (daily count of D1/deep dives)
## N2 (daily count of D2/shallow dives)
## Total_N (total daily count of dives)
## Speed_km_da (daily displacement, in km)
## Persistence_v (persistence velocity in a given direction, or speed*cos(turn angle))
## Theta (turning angle, or difference in angle between 2 locations)

# each individual is stored as an individual dataframe within a larger list, "turtle_jitter"
str(turtle_jitter)
str(turtle_jitter[[1]]) # you can select an individual with its index number

# it is good practice to visualize your variables of interest!
# especially for HMM's, where you want to specify suitable distributions for each model parameter/variable

# first take a look at the tracks:
lapply(turtle_jitter,function(z){
  ggplot(data=z,aes(lon_j,lat_j))+
    geom_path()+
    geom_point(size=1,color="red")+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position = "none")+
    xlab("Longitude")+ylab("Latitude")+ggtitle(paste("Individual",unique(z$Ref)))
})

# N1 (daily count of D1/deep dives)
lapply(turtle_jitter,function(z){
  ggplot() + 
    geom_bar(data=as.data.frame(table(z$N1)),aes(x=Var1,y=Freq),color="black",stat="identity") +
    theme_classic()+
    theme(legend.position="none")+
    xlab(paste("Count of Deep Dives, ID: ",z$Ref[1]))
})

# N2 (daily count of D1/deep dives)
lapply(turtle_jitter,function(z){
  ggplot() + 
    geom_bar(data=as.data.frame(table(z$N2)),aes(x=Var1,y=Freq),color="black",stat="identity") +
    theme_classic()+
    theme(legend.position="none")+
    xlab(paste("Count of Shallow Dives, ID: ",z$Ref[1]))
})


## Speed_km_da (daily displacement, in km)
lapply(turtle_jitter,function(z){
  ggplot() + 
    geom_histogram(data=z,aes(x=Speed_km_da),bins=30,color="black",fill="pink") +
    theme_classic()+
    theme(legend.position="none")+
    xlab(paste("Swim Speed (km/day), ID: ",z$Ref[1]))
})

## Persistence_v (persistence in a given direction, or speed*cos(turn angle))
lapply(turtle_jitter,function(z){
  ggplot() + 
    geom_histogram(data=z,aes(x=Persistence_v),bins=30,color="black",fill="lightblue") +
    theme_classic()+
    theme(legend.position="none")+
    xlab(paste("Persistence Velocity (km/day), ID: ",z$Ref[1]))
})

## Theta (turning angle, or difference in angle between 2 locations)
lapply(turtle_jitter,function(z){
  ggplot() + 
    geom_histogram(data=z,aes(x=Theta,y=..density..),position="identity",bins=36,color="black",fill=NA) +
    coord_polar(start=0)+
    scale_x_continuous("", limits = c(-3, 3), breaks = seq(-3, 3), labels = seq(-3,3))+
    theme_minimal()+
    theme(legend.position="none")+
    xlab(paste("Turning Angle (radians), ID: ",z$Ref[1]))+ylab("Density")
})

# Observations:
## N1/N2 are count data and skewed- suitable for a Poisson distribution (parameter- lambda)
## Swim speed is continuous/positive and skewed- suitable for a Weibull distribution (parameter- alpha, sigma)
## Persistence velocity is continuous/negative-positive - suitable for a Cauchy distribution (parameter- mu, sigma)
## Turning angle is continuous/negative-positive/circular- suitable for a Wrapped Cauchy distribution (parameter- rho)

## for more possible distributions, see the Stan ref manual: https://mc-stan.org/docs/2_29/functions-reference/
## (discrete and continuous distributions and associated rstan code)

##########################################################################
# FIT MODELS
##########################################################################

# Note: 
## Stan uses C++ and interfaces with R with "rstan"- however, it uses it's own language for the code files itself!
## for help and more info on creating your own Stan models and the Stan language, we recommend:
## The Appendix for Auger-Methe et al. (2020, "An introduction to state-space modeling of ecological time series")
## Leos-Barajas and Michelot (2018, "An Introduction to Animal Movement Modeling with Hidden Markov Models using Stan for Bayesian Inference")
## Stan User's Guide and Language Reference Manual (https://mc-stan.org/users/documentation/)

# check to make sure you can stan models using the code below:
example(stan_model, package = "rstan", run.dontrun = TRUE)
## if this doesn't work, you may need to install rstan from source!

# Step 1: Format your data for Stan 

# store data in "list of lists" for your ids
# select the variables you will use in your model!
# make sure they are named the same as in your stan model...
# Here we choose:
## T - number of observations
## V - move persistence velocity
## Theta - turning angle
## S - swim speed
## N1 - count of D1 (deep) dives
## N2 - count of D2 (shallow) dives

turtle_stan<-c()
for ( i in c(1:length(turtle_jitter))){
  turtle_stan[[i]]<-with(turtle_jitter[[i]], 
                         list(T = nrow(turtle_jitter[[i]]), V = Persistence_v, Theta = Theta,
                              S = Speed_km_da,
                              N1 = N1, N2 = N2))
}

# always check your data to see if it looks ok
str(turtle_stan[[1]])

# Step 2: Run your Models!

# run the below code first to:
## run model in parallel on multiple cores
options(mc.cores=parallel::detectCores())
## prevent recompiling stan program
rstan_options(auto_write=TRUE)

# model options:
## model 1: 3-state HMM with variables for persistence velocity and dive type
## model 2: 3-state HMM with variables for theta, speed, and dive type

## there is a "Viterbi" version for each model- this allows you to predict behavior states for your observed data
## for model selection/comparison, we will use the non-Viterbi versions:

# let's try fitting model1 to ID1 as example:
model1_example<-stan(file="./Stan_CodeFiles/Model1.stan",
                     data=turtle_stan[[1]],
                     chains=4,iter=2000,control=list(adapt_delta=0.95,max_treedepth=15),seed="123")
# NOTE: ignore this warning message if you get it, "In readLines(file, warn = TRUE): incomplete final line found"

# you can get a quick view of the model convergence with "plot":
## use the "pars" argument to see specific parameters...
plot(model1_example,pars=c("mu","sigma","Pr"))

# now let's fit the models to all IDs at once and store in object of lists (model1, model2):
## Note: if you have a lot of individuals (> 20), this can take 30+ minutes!
model1<-c()
for (i in c(1:length(turtle_stan))){
  file.remove("./Stan_CodeFiles/Model1.rds") # this code removes intermediate rds files that can mess up model fitting...make sure it reads to the file path of your Stan code files...
  model1[[i]]<-stan(file="./Stan_CodeFiles/Model1.stan",data=turtle_stan[[i]],
                    chains=4,iter=2000,control=list(adapt_delta=0.95,max_treedepth=15),seed="123")
}

# model 2
model2<-c()
for (i in c(1:length(turtle_stan))){
  file.remove("./Stan_CodeFiles/Model2.rds") # this code removes intermediate rds files that can mess up model fitting...make sure it reads to the file path of your Stan code files...
  model2[[i]]<-stan(file="./Stan_CodeFiles/Model2.stan",data=turtle_stan[[i]],chains=4,iter=2000,control=list(adapt_delta=0.95,max_treedepth=15),seed="123")
}

###################################################################################
# MODEL COMPARISON AND SELECTION
##################################################################################

# NOTE:
# it's important to remember that there is no such thing as a perfect model!
# just "better" ones, and ultimately dependent on your data and research questions
# but if all of your models look reasonable, quantitative comparisons can help you choose!

# We will be using Leave-One-Out Cross Validation (LOO) 
## this is recommended instead of AIC/DIC/BIC for Bayesian models
## see Vehtari et al. 2017, 2019 and the "loo" R package documentation and vignettes

## vignette: https://mc-stan.org/loo/articles/loo2-with-rstan.html

# make sure it's installed!
library(loo)

# In order to use the loo package, we need to extract the log-lik of each observation at each MCMC sample 
## this is added as code in the "generated quantities" section of our non-Viterbi Stan model versions
## also see: https://discourse.mc-stan.org/t/waic-for-any-stanfit-object/10724/5
## example of how to add log lik to code: https://www.r-bloggers.com/2019/05/bayesian-modeling-using-stan-a-case-study/

# Let's compare models for ID1 as example
# make a list of models to compare:
model_list <- list(model1[[1]],model2[[1]])

# then we extract the log likelihood at each MCMC sample and determine the relative effective sample size (MCMC effective sample size/total sample size)
log_lik_list <- lapply(model_list, extract_log_lik)

r_eff_list <- lapply(model_list, function(x) {
  ll_array <- extract_log_lik(x, merge_chains = FALSE)
  relative_eff(exp(ll_array))
})

# we then compute the PSIS-LOO (Pareto-Smoothed Importance Sampling Leave-One-Out) using the array of log-lik values
loo_list <- lapply(1:length(log_lik_list), function(j) {
  loo(log_lik_list[[j]], r_eff = r_eff_list[[j]], cores = 2)
})

# we can also compute the waic scores for comparison:
waic_list <- lapply(1:length(log_lik_list), function(j) {
  waic(log_lik_list[[j]])
})

# now we compare the model using:
## pairwise comparisons between each model and the model with the largest ELPD (model in first row)
## first row: difference between the preferred model and itself
## other rows: difference between the preferred model and other models
comp_loo <- loo_compare(loo_list)
print(comp_loo) # can set simplify=FALSE for more detailed print output
# model 1 is preferred

comp_waic <- loo_compare(waic_list)
print(comp_waic)

# we can also look at the elpd scores for each model/id:
loo_list[[1]]
loo_list[[2]]
# NOTE: model 2 has 33% of Pareto k estimates for ID1 that are "very bad"
## can use k-fold cross validation to improve the PSIS estimate


# we can also compare models using:
## stacking weights, to stack predictive distributions using Bayesian bootstrap
## combines all models by maximizing the leave-one-out predictive density of the combination distribution. 
## aka, it finds the optimal linear combining weights for maximizing the leave-one-out log score
loo_model_weights(log_lik_list, method = 'stacking', r_eff_list = r_eff_list)

# stacked pseudo bayes factor: ratio of 2 weights to how much more weight is given compared to other
## 

# results across the board are in favor of model 1 for ID1- model 2 is problematic

## in order to get behavior state predictions for each location in our tracks, we need to apply the Viterbi algorithm
## this allows us to determine the most likely sequence of behavior states, given our chosen model and data
model_final<-c()
for (i in c(1:length(turtle_stan))){
  file.remove("./Stan_CodeFiles/Model1_Viterbi.rds") # this code removes intermediate rds files that can mess up model fitting...make sure it reads to the file path of your Stan code files...
  model_final[[i]]<-stan(file="./Stan_CodeFiles/Model1_Viterbi.stan",data=turtle_stan[[i]],
                    chains=4,iter=2000,control=list(adapt_delta=0.95,max_treedepth=15),seed="123")
}

# save your fitted models for the next steps: visualizations and diagnostics!
save(model_final,file="./Example_Data/final_HMM_model.rds")

#########################################################################
# HINTS
########################################################################

# to assist with convergence and model fit, you can "play" with various options:
## place an ordering on some parameters for certain states
## use more/less informative priors
## try different distributions or versions for your variables (e.g. using persistence velocity instead of swim speed/theta)
## fit models to your "best" track or to simulated data to get your model working!
## when in doubt, choose parsimony for your models!
## although a pain/time-consuming to fit and compare, fitting models to one individual at a time allows you to cater models to each individual and diagnose convergence issues (natural populations often have variety in behavior not always captured by heirarchical models)
## but hierarchical models can also be created to "borrow strength" across individuals and get population-level predictions





