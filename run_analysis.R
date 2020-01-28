library(dplyr)
library(bbmle)
library(EpiEstim)
library(R0)
library(deSolve)
library(ggplot2)
library(truncnorm)

# specify repository path 
repo_path <- 'C:/Users/Megan/Documents/Github/R0-methods-comparison'
setwd(repo_path)

# source functions for estimating R0 & data simulation
source('methods.R')
source('fit_R0_seq.R')
source('data_simulation.R')


#===== 1. Fit to simulated data with varying levels of random noise =====#

# simulate 250 datasets, assuming average infectious and incubation periods of 6 and 14 days
simulations <- simulate_data(N=250, gamma=1/6, sigma=1/14)

# define the generation time distribution
GTd <- R0::generation.time("gamma", c(20/7, 7.4/7))

# list to store results from simulations
Sim_results <- list()

# fit each method at sequential time points of epidemic growth phase
# note: running this locally in a loop will take a long time.
# where possible, we recommend using parallel computing for this.

# for each of the 3 noise levels
for(i in 1:3){ 
  
  # for each of the 250 simulations
  for(j in 1:ncol(simulations$sims)){ 
    
    # list for storage
    results_i <- list()
    
    # data for fitting - simulation j with noise level i
    sim_ij <- data.frame(week=1:nrow(simulations[[i]]), cases=simulations[[i]][ ,j], 
                         country=paste(names(simulations[i]), j, sep='_'))
    # fit & store results
    results_i[[j]] <- fit_R0_seq(data=sim_ij, mean_GT=20/7, sd_GT=7.4/7, GTd=GTd)
    
  }
  
  # store results from each noise level
  Sim_results[[i]] <- results_i
}




#===== 2. Fit to empirical Zika outbreak data from LAC national surveillance data =====#

# read in LAC national Zika surveillance data
zika_LAC <- readRDS('ZikaCases_LatinAmerica&Caribbean.RDS')

# exclude the countries where cases peaked earlier than 6 weeks (2 generation times)
# (Antigua and Barbuda, Nicaragua, Venezuela)
zika_LAC <- zika_LAC[c(2:22, 24:34, 36:37)]

# list to store results 
LAC_results <- list()

# fit each method at sequential time points of epidemic growth phase
# note: running this locally in a loop will take a long time.
# where possible, we recommend using parallel computing for this.

# for each country 
for(i in 1:length(zika_LAC)){
  
  # extract data for fitting
  data_i <- zika_LAC[[i]]
  
  # define the generation time distribution
  GTd <- R0::generation.time("gamma", c(20/7, 7.4/7))
  
  # fit & store results
  LAC_results[[i]] <- fit_R0_seq(data=data_i, mean_GT=20/7, sd_GT=7.4/7, GTd=GTd)
  
  print(paste("Finished", data_i$country[1], sep=" "))
}
