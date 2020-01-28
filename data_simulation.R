#=== Function to simulate N epidemic curves with 3 levels of noise:
# - no noise
# - Poisson random errors
# - Negative Binomial random errors

simulate_data <- function(N, gamma, sigma){
  
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
  
  
  
  #=== Generate N simulated epidemic curves with randomly selected parameter values
  
  # randomly select values of R0, population size (Np), and inital number infected (I)
  R0v <- round(rtruncnorm(N, a=1.2, b=8, mean=2, sd=2),2) 
  Npv <- 10^(round(runif(N, min=3, max=6),0))  
  Iv <- round(runif(N, min=1, max=5),0)
  
  # daily time step for 2 years
  timestep <- seq(0, 2*365) 
  
  # assuming that 20% of infections result in reported cases 
  Reporting_fraction <- 0.20
  
  # dataframes to store simulations
  sim_data <- data.frame(matrix(ncol=N, nrow=105))
  colnames(sim_data) <- paste("sim", seq(1:N), sep="")
  sim_data_P <- sim_data
  sim_data_NB <- sim_data
  
  # store parameter values
  pars_data <- data.frame(sim=1, R0=R0v[1], Np=Npv[1], I0=Iv[1])
  
  
  
  #=== Simulate SEIR model & aggregate to weekly incidence data
  for(i in 1:N){
    
    # define parameters & initial states
    pars <- c(R0=R0v[i], gamma=gamma, sigma=sigma)
    state_init <- c(S=Npv[i], E=0, I=Iv[i], R=0) 
    
    # run the SEIR model
    out <- ode(y=state_init, times=timestep, func=model, parms=pars, method="rk4") 
    
    # extract the incidence of infections
    incidence <- round(pars[which(names(pars)=="sigma")] * out[,which(colnames(out)=="E")]) 
    
    # aggregate to weekly data
    Inc_weekly <- round(unname(tapply(incidence, (seq_along(incidence)-1) %/% 7, sum))* Reporting_fraction) 
    
    # cut the time series to begin at time of 1st case
    indx <- which(Inc_weekly>0)
    Inc_weekly <- Inc_weekly[indx[1]:indx[length(indx)]]
    
    # add reporting noise to the data (after week 1)
    inc_P <- c(Inc_weekly[1], add_noiseP(Inc_weekly[2:length(Inc_weekly)]))
    inc_NB <- c(Inc_weekly[1], add_noiseNB(Inc_weekly[2:length(Inc_weekly)]))
    
    # store weekly incidence data and true parameter values
    sim_data[1:length(Inc_weekly) ,i] <- Inc_weekly
    sim_data_P[1:length(Inc_weekly),i] <- inc_P 
    sim_data_NB[1:length(Inc_weekly),i] <- inc_NB
    pars_data[i, ] <- c(i, R0v[i], Npv[i], Iv[i])
    
  }
  
  # return a list of simulated data, with varying noise levels, & parameter values used
  return(list(sims=sim_data, simsP=sim_data_P, simsNB=sim_data_NB, pars=pars_data))
}
  



