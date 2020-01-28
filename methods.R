
#=== 1. Linear exponential growth method
EG_Lin <- function(data){
  
  # log transform cases, add small constant to time-series if zeros are present
  if(any(data$cases==0)){
    data$casesC <- data$cases+0.5
    data$logcases <- log(data$casesC)
  }else{
    data$logcases <- log(data$cases)
  }
  
  # define growth rate, r
  lm_mod <- lm(logcases~ week, data=data)
  r <- coef(lm_mod)[2]
  
  if(r>0){ 
    
    # calculate R by integrating over the GT distribution, g(a)
    int_ga <- integrate(Mz_ga, r=r, lower=0, upper=Inf, subdivisions=1e8)
    R <- 1/int_ga$value
    
    #=== calculate 95% CI for R
    
    # stanard error of r
    se_r <- coef(summary(lm_mod))[, "Std. Error"][2]
    
    # bootstrap sample from t-distribution
    boot <- vector()
    for(k in 1:10000){
      rdf <- rt(1, df=(lm_mod$df.residual))
      r_new <- abs(r + se_r*rdf) 
      int_ga_new <- integrate(Mz_ga, r=r_new, lower=0, upper=Inf, subdivisions=1e8)
      boot[k] <- 1/int_ga_new$value
    }
    results <- c(R, quantile(boot, c(0.025,0.975), na.rm=TRUE))
  }else{
    results <- c(NA, NA, NA)
  }
  return(results)
}



#=== 2. Poisson exponential growth method
EG_P <- function(data){
  
  # define growth rate, r
  P_mod <- glm(cases~ week, data=data, family="poisson")
  r_p <- coef(P_mod)[2]
  
  if(r_p>0){
    
    # calculate R by integrating over the GT distribution, g(a)
    int_ga <- integrate(Mz_ga, r=r_p, lower=0, upper=Inf, subdivisions=1e8)
    R_p <- 1/int_ga$value
    
    #=== calculate 95% CI for R
    
    # standard error of r
    se_r_p <- coef(summary(P_mod))[, "Std. Error"][2]
    
    # bootstrap sample from t-distribution
    R0exp_pois_CI <- vector()
    for(k in 1:10000){
      rdf <- rt(1, df=(P_mod$df.residual))
      r_new <- abs(r_p + se_r_p*rdf)
      int_ga_new <- integrate(Mz_ga, r=r_new, lower=0, upper=Inf, subdivisions=1e8)
      R0exp_pois_CI[k] <- 1/int_ga_new$value
      
    }
    results <- c(R_p, quantile(R0exp_pois_CI, c(0.025, 0.975), na.rm=TRUE))
  }else{
    results <- c(NA, NA, NA)
  }
  return(results)
}


#=== 3. Maximum likelihood exponential growth method

EG_MLE <- function(data){
  
  # define log-likelihood (negative required by mle2)
  ll <- function(I0, r){
    -sum(data$cases*(r*data$week+log(I0))-I0*exp(r*data$week)-lfactorial(data$cases))
  }
  
  # vectors & matrix for storing r, likelihood & CI estimates
  r_est <- vector()
  likel <- vector()
  CI_mle <- data.frame(matrix(ncol=2, nrow=100))
  
  # run MLE from different starting points
  for(m in 1:250){
    
    tryCatch({ # catch errors arising from negative r values and filter these before storing results
      est <- mle2(ll, start=list(I0=runif(1,1,10), r=runif(1,1,10)))
      
      # store values of r, likelihood & CI
      r_est[m] <- coef(est)[2]
      likel[m] <- logLik(est)
      CI_mle[m, ] <- c(confint(est, quietly=TRUE)[2,1], confint(est, quietly=TRUE)[2,2])
      rm(est)
    }, error=function(e){})
  } 
  
  # take maximum likelihood estimates
  r_mle <- r_est[which.max(likel)]
  ci_mle <- CI_mle[which.max(likel), ]
  
  # check the validity of results (i.e. not NA & are above min_r)
  if(is.na(ci_mle[1])+is.na(ci_mle[2])==0 & ci_mle[1]>-0.38){
    
    # calculate R by integrating over the GT distribution, g(a)
    R <- 1/(integrate(Mz_ga, r=r_mle, lower=0, upper=Inf, subdivisions=1e8))$value
    R_CIL <- 1/(integrate(Mz_ga, r=CI_mle[1,1], lower=0, upper=Inf, subdivisions=1e8))$value
    R_CIU <- 1/(integrate(Mz_ga, r=CI_mle[1,2], lower=0, upper=Inf, subdivisions=1e8))$value
    
    results <- c(R, R_CIL, R_CIU)
    
  }else{
    results <- c(NA, NA, NA)
  }
  return(results)
}


#=== 4. EpiEstim method

epiestim <- function(data, mean_si, std_si){
  
  R_EpiEstim <- EpiEstim::estimate_R(data$cases, method="parametric_si", 
                                     config=make_config(list(mean_si=20/7, std_si=7.4/7,
                                                             t_start=2, t_end=as.numeric(nrow(data)))))
  if(exists("R_EpiEstim")==TRUE){
    results <- c(R_EpiEstim$R$`Mean(R)`, R_EpiEstim$R$`Quantile.0.025(R)`, R_EpiEstim$R$`Quantile.0.975(R)`)
  } else {
    results <- c(NA, NA, NA)
  }
  return(results)
}


#=== 5. White and Pagano method

WP <- function(data, GTd){
  
  tryCatch({
    R0_WP <- R0::est.R0.ML(data$cases, GT=GTd, begin=1, end=as.numeric(nrow(data)))
  }, error=function(e){})
  
  if(exists("R0_WP")==TRUE){
    results <- c(R0_WP$R, R0_WP$conf.int)
  }else{
    results <- c(NA, NA, NA)
  }
}