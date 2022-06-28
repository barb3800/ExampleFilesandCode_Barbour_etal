## EXAMPLE: HMM DIAGNOSTICS, PREDICTIONS, AND VISUALS

## Author and Contact Info: Nicole Barbour, nbarbour@umd.edu 
## Code as example for Barbour et al., "Incorporating multidimensional behavior into a dynamic management tool for a critically endangered and migratory species"

##############################################
# SET UP 
##############################################

# load libraries (note: all these packages need to be installed before use)
library(rstan);library(ggplot2);library(gridExtra);library(dplyr)

# set working directory
setwd("yourfilepathhere")

# read in data (.rda format)
## chosen fitted HMM model
load(file="./Example_Data/final_HMM_model.rds")
## turtle tracking data
load("./Example_Data/Jitter_Turtle_Data.rda")

#############################################
# VISUAL DIAGNOSTICS: CONVERGENCE PLOTS
############################################

# Parameter means/sd's 
# Trace plots of Markov chains (mu, sigma, Pr)
# Histograms for state-level parameters (mu, sigma, Pr)

# add plots to lists
parameter_plot<-c()
trace_plot<-c()
hist_plot<-c()
# individual level trace and parameter plots
for (i in c(1:length(model_final))){
  # plot of parameter means and sd
  parameter_plot[[i]]<-plot(model_final[[i]],pars=c("mu", "sigma","Pr"))+ggtitle(paste("Track",i))
  # traceplots for chains for each parameter
  trace_plot[[i]]<-traceplot(model_final[[i]],pars=c("mu", "sigma","Pr"),inc_warmup=FALSE)+ggtitle(paste("Track",i))
  # histograms of parameters for each state
  hist_plot[[i]]<-stan_hist(model_final[[i]],pars=c("mu", "sigma","Pr"),alpha=0.5)+ggtitle(paste("Track",i))
}

# arrange into a grid to compare plots across ids
# specify no. of columns/rows
marrangeGrob(parameter_plot,nrow=3,ncol=1)
marrangeGrob(trace_plot,nrow=3,ncol=1)
marrangeGrob(hist_plot,nrow=3,ncol=1)

# can also print plots individually:
trace_plot[[1]]

## visualizations so far show good convergence of chains and separation of means/sds for parameters between states

#############################################
# VISUAL DIAGNOSTICS: DENSITY PLOTS
############################################

# density histograms of parameters for each behavior state

# first we create a dataframe of the parameters for each id:
# extract data for each parameter for each state and put into a dataframe
# Note: we also extract the parameter "Rhat"- Rhat > 1.05 indicates model convergence issues
turtle_parameter<-c()
for ( i in c(1:length(model_final))){
  # extract chain estimates for each parameter and state
  mu1<-rstan::extract(model_final[[i]],pars=c("mu"))[[1]][,1]
  mu2<-rstan::extract(model_final[[i]],pars=c("mu"))[[1]][,2]
  mu3<-rstan::extract(model_final[[i]],pars=c("mu"))[[1]][,3]
  sigma1<-rstan::extract(model_final[[i]],pars=c("sigma"))[[1]][,1]
  sigma2<-rstan::extract(model_final[[i]],pars=c("sigma"))[[1]][,2]
  sigma3<-rstan::extract(model_final[[i]],pars=c("sigma"))[[1]][,3]
  lambda<-rstan::extract(model_final[[i]],pars=c("lambda"))[[1]]
  Pr1<-rstan::extract(model_final[[i]],pars=c("Pr"))[[1]][,1]
  Pr2<-rstan::extract(model_final[[i]],pars=c("Pr"))[[1]][,2]
  Pr3<-rstan::extract(model_final[[i]],pars=c("Pr"))[[1]][,3]
  
  # add column for id and for index number
  id<-turtle_jitter[[i]]$Ref[1]
  index_no<-i
  # extract max Rhat for all parameters (remove infinite values)
  Rhat<-summary(model_final[[i]])$summary[,"Rhat"]
  Rhat<-max(Rhat[is.finite(Rhat)],na.rm=TRUE)
  # put parameters into a dataframe and add to list
  turtle_parameter[[i]]<-data.frame(mu1,mu2,mu3,sigma1,sigma2,sigma3,lambda,Pr1,Pr2,Pr3,id,index_no,Rhat)
}

# bind all turtle parameter and jittered data frames together into one
turtle_parameter_all<-do.call("rbind",turtle_parameter)
turtle_jitter_all<-do.call("rbind",turtle_jitter)

## extract mean parameters for each state
lambda<-mean(turtle_parameter_all$lambda)
Pr1<-mean(turtle_parameter_all$Pr1) # prob of D1
Pr2<-mean(turtle_parameter_all$Pr2) # prob of D1
Pr3<-mean(turtle_parameter_all$Pr3) # prob of D1

# find probability of dives either being d1 (deep) or d2 (shallow)
lambda1_d1<-as.integer(round(lambda*Pr1)) # prob of D1 in state 1
lambda2_d1<-as.integer(round(lambda*Pr2)) # prob of D1 in state 2
lambda3_d1<-as.integer(round(lambda*Pr3)) # prob of D1 in state 3

lambda1_d2<-as.integer(round(lambda*(1-Pr1))) # prob of D2 in state 1
lambda2_d2<-as.integer(round(lambda*(1-Pr2))) # prob of D2 in state 2
lambda3_d2<-as.integer(round(lambda*(1-Pr3))) # prob of D2 in state 3

# make data frame of counts (N1- count of D1/deep dives, N2-count of D2/shallow dives)
df_N1<-as.data.frame(table(turtle_jitter_all$N1))
df_N2<-as.data.frame(table(turtle_jitter_all$N2))


# make stepwise grids of the counts (N1/N2) of the different dive types 
# this is the data that will be used to create the plots
stepgrid_N1 <- as.integer(round(seq(0,
                                    max(as.numeric(df_N1$Var1), na.rm = TRUE),
                                    length = max(as.numeric(df_N1$Var1))))) 

stepgrid_N2 <- as.integer(round(seq(0,
                                    max(as.numeric(df_N2$Var1), na.rm = TRUE),
                                    length = max(as.numeric(df_N2$Var1))))) 

# create poisson density curves for each state:

## for N1: poisson curves for states 1-3
df_N1$pois_1<-dpois(x=stepgrid_N1, lambda = lambda1_d1) # prob of D1 in state 1
df_N1$pois_2<-dpois(x=stepgrid_N1, lambda = lambda2_d1) # prob of D1 in state 2
df_N1$pois_3<-dpois(x=stepgrid_N1, lambda = lambda3_d1) # prob of D1 in state 2

df_N2$pois_1<-dpois(x=stepgrid_N2, lambda = lambda1_d2) # prob of D2 in state 1
df_N2$pois_2<-dpois(x=stepgrid_N2, lambda = lambda2_d2) # prob of D2 in state 2
df_N2$pois_3<-dpois(x=stepgrid_N2, lambda = lambda3_d2) # prob of D2 in state 3

# density plot for D1 dives
## state 1- pink, state 2- orange, state 3- blue
ggplot() +
  geom_histogram(data=turtle_jitter_all,aes(x=N1,y=..density..),color="black",fill="grey",bins=20)+
  geom_histogram(color="black",stat="identity") +
  xlab("Daily Count of D1 (Deep) Dives")+ylab("Density")+
  geom_line(data=df_N1, aes(y = pois_1,x=as.numeric(Var1)), col = "#CC79A7",size=1) +
  geom_line(data=df_N1, aes(y = pois_2,x=as.numeric(Var1)), col = "#E69F00",size=1) +
  geom_line(data=df_N1, aes(y = pois_3,x=as.numeric(Var1)), col = "#0072B2",size=1) +
  theme_classic()

# density plot for D2 dives
## state 1- pink, state 2- orange, state 3- blue
ggplot() +
  geom_histogram(data=turtle_jitter_all,aes(x=N2,y=..density..),color="black",fill="grey",bins=20)+
  geom_histogram(color="black",stat="identity") +
  xlab("Daily Count of D2 (Shallow) Dives")+ylab("Density")+
  geom_line(data=df_N2, aes(y = pois_1,x=as.numeric(Var1)), col = "#CC79A7",size=1) +
  geom_line(data=df_N2, aes(y = pois_2,x=as.numeric(Var1)), col = "#E69F00",size=1) +
  geom_line(data=df_N2, aes(y = pois_3,x=as.numeric(Var1)), col = "#0072B2",size=1) +
  theme_classic()

# do something similar to create density plots for persistence velocity:
## extract mean parameters 
location1<-mean(turtle_parameter_all$mu1)
location2<-mean(turtle_parameter_all$mu2)
location3<-mean(turtle_parameter_all$mu3)

scale1<-mean(turtle_parameter_all$sigma1)
scale2<-mean(turtle_parameter_all$sigma2)
scale3<-mean(turtle_parameter_all$sigma3)

# stepwise data of persistence v values
stepgrid <- seq(min(turtle_jitter_all$Persistence_v, na.rm = TRUE),
                max(turtle_jitter_all$Persistence_v, na.rm = TRUE),
                length = 1000)

# density plot: 
ggplot()+
  geom_histogram(data=turtle_jitter_all,aes(x=Persistence_v,y=..density..),color="black",fill="grey",bins=20)+
  xlab("Move Persistence Velocity (km/day)")+ylab("Density")+ stat_function(data=data.frame(x=stepgrid),fun=dcauchy,args=list(location=location1,scale=scale1),colour="#CC79A7",size=1)+
  stat_function(data=data.frame(x=stepgrid),fun=dcauchy,args=list(location=location2,scale=scale2),colour="#E69F00",size=1)+
  stat_function(data=data.frame(x=stepgrid),fun=dcauchy,args=list(location=location3,scale=scale3),colour="#0072B2",size=1)+
  theme(text=element_text(size=10),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

#############################################
# VISUAL DIAGNOSTICS: SPATIO-TEMPORAL PLOTS
############################################

# To further visualize spatial and temporal trends (which is important if you are working with spatio-temporal data!)
## we first need to use the results of our Viterbi algorithm to predict behavior states for each data observation

# we extract the predicted behavior states from our Viterbi algorithm and bind the data back to our original data so we have a dataframe of positions, ids, and predicted behavior states
for (i in c(1:length(model_final))){
  # extract Viterbi state estimates for each observation in time series:
  Viterbi<-rstan::extract(model_final[[i]],pars=c("viterbi"))
  
  # extract state estimates: take mean (round up) of state estimates for each observation from iterations (n=4000)
  state_est<-ceiling(colMeans(Viterbi$viterbi))
  
  # bind state estimates back to data
  turtle_jitter[[i]]$States<-state_est
}

# bind data across ids again to make one dataframe:
turtle_jitter_all<-do.call("rbind",turtle_jitter)

# latitude density plots for each state:
ggplot(turtle_jitter_all)+
  xlab("Latitude (deg)")+
  geom_density(aes(x=lat_j,fill=as.factor(States),color=as.factor(States)),alpha=.4)+
  scale_fill_manual(values=c("#CC79A7","#E69F00","#0072B2"))+
  scale_color_manual(values=c("#CC79A7","#E69F00","#0072B2"))+
  theme_classic()+
  facet_wrap(~States)

# proportion of individuals in each state over time (months):
turtle_jitter_all$Month<-as.factor(format(turtle_jitter_all$date,"%m")) # create  month variable from posixct date column

ggplot(turtle_jitter_all, aes(Month, fill = as.factor(States))) +
  geom_bar(position = "fill")+
  xlab("Month")+ylab("Proportion of Individuals in Each State")+
  scale_fill_manual(values=c("#CC79A7","#E69F00","#0072B2"))+
  theme_classic()

# proportion of different behavior states for each id:
agg<-count(turtle_jitter_all,States,Ref)

ggplot()+
  geom_col(data=agg,aes(x=Ref,y=n,fill=as.factor(States)),position="fill",width=0.5)+
  theme_classic()+
  xlab("ID")+ylab("Proportion")+labs(fill="State")+
  scale_fill_manual(values=c(c("#CC79A7","#E69F00","#0072B2")))

# time series of predicted states for each id: 
## first create a new column with the observation number for each day for each id
turtle_jitter_all<- turtle_jitter_all %>% group_by(Ref) %>% mutate(obs_no=cumsum(rep(1,length(date)))) %>% ungroup()
ggplot()+
  geom_point(data=turtle_jitter_all,aes(x=obs_no,y=as.numeric(States),color=States))+
  geom_path(data=turtle_jitter_all,aes(x=obs_no,y=as.numeric(States)))+
  scale_colour_gradient2(low="#CC79A7",high="#0072B2",mid="#E69F00",midpoint=2)+
  xlab("Time (days)")+ylab("States")+
  theme_classic()+
  facet_wrap(~Ref,scales="free",ncol=3)

#############################################
# QUANTITATIVE DIAGNOSTICS
############################################

# we can also do some quantitative diagnostics to assess our model fit 
# for reference, also see Auger Methe et al. (2020), Appendix S1
## and Leos-Barajas and Michelot 2018
## and a Brief Guide to Stan Warnings: https://mc-stan.org/misc/warnings.html 
## and ShinyStan: https://cran.r-project.org/web/packages/shinystan/shinystan.pdf, http://mc-stan.org/shinystan/articles/shinystan-package.html

# Look at the print out for the model fits with the "print()" function:
## allows you to print model results, obtain quantiles from the marginal posterior distributions of the parameters,
## and get estimated values of R-hat (want value < 1.1 to determine convergence)
print(model_final[[1]]) # example, with id 1

# check HMC specific diagnostics:
# Divergences:  results from step size of sampler being too large to capture features of target dist-- can fix by increasing value of adapt_delta parameter or reparameterize the model
# Tree depth: configuring HMC in STan involves placing max cap on tree depth- fix by increasing max tree depth
# Energy Bayesian Fraction of Missing Information (BFMI): indicates adaptation phase of Markov chains did not explore posterior dist efficiently- fix by reparameterizing model or setting higher iter/warmup values
check_hmc_diagnostics(model_final[[1]])

# can also use "shinystan" package to create interactive plots for posterior samples
library(shinystan)

launch_shinystan_demo() # demo of evaluating shiny stan object

# create vector of observations for posterior predictive checks
Y<-length(turtle_jitter[[1]]$date) # is this right? get error than this should match length of model parameter (e.g. alpha)

# apply to stan object
shinystan_turtle<-launch_shinystan(model_final[[1]])

############################################################
# SAVE DATA
############################################################

# the last step is to save your data with with predicted behavior states:
turtle_jitter_states<-turtle_jitter_all

save(turtle_jitter_states,file="./Example_Data/turtle_jitter_states.rda")


