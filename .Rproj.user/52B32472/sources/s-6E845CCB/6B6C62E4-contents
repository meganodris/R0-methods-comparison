library(deSolve)
library(ggplot2)
library(grid)
library(truncnorm)


#=== Define the SEIR model

# gamma = 1 / average infectious period
# sigma = 1 / average incubation period

model <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    
    dS <- -(R0*gamma) * S * I/(S + E + I + R)
    dE <- (R0*gamma) * S * I/(S + E + I + R) - sigma * E 
    dI <- sigma * E - gamma * I 
    dR <- gamma * I
    
    return(list(c(dS, dE, dI, dR)))
    
  })
}



#=== Functions to generate Poisson and Negative Binomial random errors

# Poisson 
add_noiseP <- function(inc){ 
  res <- round(rpois(length(inc), lambda=inc))
  return(res)
}

# Negative binomial (over-dispersed)
add_noiseNB <- function(inc){ 
  res <- round(rnbinom(length(inc), size=2, mu=inc))
  return(res)
}



#=== Generate 250 simulated epidemic curves with randomly selected parameter values

# number of datasets to be simulated
n <- 250 

# randomly select values of R0, population size (N), and inital number infected (I)
R0v <- round(rtruncnorm(n, a=1.2, b=8, mean=2, sd=2),2) 
Nv <- 10^(round(runif(n, min=3, max=6),0))  
Iv <- round(runif(n, min=1, max=5),0)

# daily time step for 2 years
timestep <- seq(0, 2*365) 

# assuming that 20% of infections result in reported cases 
Reporting_fraction <- 0.20

# dataframes to store simulations
sim_data <- data.frame(matrix(ncol=n, nrow=105))
colnames(sim_data) <- paste("sim", seq(1:n), sep="")
sim_data_P <- sim_data
sim_data_NB <- sim_data

# store parameter values
pars_data <- data.frame(sim=1, R0=R0v[1], N=Nv[1], I0=Iv[1])



#=== Simulate SEIR model & aggregate to weekly incidence data
for(i in 1:n){
  
  # define parameters & initial states
  pars <- c(R0=R0v[i], gamma=1/6, sigma=1/14)
  state_init <- c(S=Nv[i], E=0, I=Iv[i], R=0) 
  
  # run the SEIR model
  out <- ode(y=state_init, times=timestep, func=model, parms=pars, method="rk4") 
  
  # extract the incidence of infections
  incidence <- round(pars[which(names(pars)=="sigma")] * out[,which(colnames(out)=="E")]) 
  
  # aggregate to weekly data
  Inc_weekly <- round(unname(tapply(incidence, (seq_along(incidence)-1) %/% 7, sum))* Reporting_fraction) 
  
  # cut the time series to 1st & last week of reported cases
  indx <- which(Inc_weekly>0)
  Inc_weekly <- Inc_weekly[indx[1]:indx[length(indx)]]
  
  # add reporting noise to the data (after week 1)
  inc_P <- c(Inc_weekly[1], add_noiseP(Inc_weekly[2:length(Inc_weekly)]))
  inc_NB <- c(Inc_weekly[1], add_noiseNB(Inc_weekly[2:length(Inc_weekly)]))
  
  # store weekly incidence data and true parameter values
  sim_data[1:length(Inc_weekly) ,i] <- Inc_weekly
  sim_data_P[1:length(Inc_weekly),i] <- inc_P 
  sim_data_NB[1:length(Inc_weekly),i] <- inc_NB
  pars_data[i, ] <- c(i, R0v[i], Nv[i], Iv[i])
  
}
