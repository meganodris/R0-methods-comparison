#====== Function to estimate time-constant R0 at sequential time points of the epidemic growth phase ======#

R_estimates <- function(data, path){
  
  data <- as.data.frame(data)
  
  # define a fixed mean and SD of the generation time distribution, GT
  mean_GT <- 20/7
  sd_GT <- 7.4/7
  variance_GT <- sd_GT^2
  GTd <- R0::generation.time("gamma", c(mean_GT, sd_GT))
  
  # define the moment generating function of the GT distribution, g(a)
  Mz_ga <- function(t, r){
    dgamma(t, shape=(mean_GT^2/variance_GT), scale=(variance_GT/mean_GT))*exp(-t*r)
  }
  
  # define minimum value of r that can be integrated in the moment generating function, Mz_ga
  min_r <- -0.38
  
  # dataframe to store results
  results <- data.frame(Country=NA, Nweeks=NA, method=NA, R0=NA, CI_L=NA, CI_U=NA)
  
  # define peak - max time point of highest reported cases
  peak <- max.col(matrix(data$cases,nrow=1),"last") 

  
  # if peak is >=6 weeks after 1st case report, run analysis
  if(peak>=6){
    
    # define sections: from 6 weeks, by 3, up to peak (i.e. 6,9,12,15...)
    sections <- seq(from=6, to=peak, by=3)
    
    count <- 0
    
    # Loop to fit each method to each section with increasing number of time points
    for(j in 1:length(sections)){
  
      # extract section of the epidemic curve for fitting
      s <- data[1:sections[j], ]
      
      #======== 1. Linear Exponential Growth ========#
      
      # log transform cases, add small constant to time-series if zeros are present
      if(any(s$cases==0)){
        s$casesC <- s$cases+0.5
        s$logcases <- log(s$casesC)
      }else{
        s$logcases <- log(s$cases)
      }
      
      # define growth rate, r
      lm_mod <- lm(logcases~ week, data=s)
      r <- coef(lm_mod)[2]
      
      if(r>min_r){ 
        
        # calculate R by integrating over the GT distribution, g(a)
        int_ga <- integrate(Mz_ga, r=r, lower=0, upper=Inf, subdivisions=1e8)
        R <- 1/int_ga$value
        
        if(R>1.00000000000206){ # equivalent of r=0
          
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
          
          # store results
          count <- count+1
          results[count, ] <- c(s$country[1], sections[j], "ExpLin", R, quantile(boot, c(0.025,0.975), na.rm=TRUE))
        }else{
          count <- count+1
          results[count, ] <- c(s$country[1], sections[j], "ExpLin", NA, NA, NA)
        }
      }else{
        count <- count+1
        results[count, ] <- c(s$country[1], sections[j], "ExpLin", NA, NA, NA)
      }
      
      
      #======== 2. Poisson Exponential Growth ========#
      
      # define growth rate, r
      P_mod <- glm(cases~ week, data=s, family="poisson")
      r_p <- coef(P_mod)[2]
      
      if(r_p>min_r){
        
        # calculate R by integrating over the GT distribution, g(a)
        int_ga <- integrate(Mz_ga, r=r_p, lower=0, upper=Inf, subdivisions=1e8)
        R_p <- 1/int_ga$value
        
        if(R_p>1.00000000000206){
          
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
          
          # store results
          count <- count+1
          results[count, ] <- c(s$country[1], Nweeks=sections[j], "ExpPois", R_p, quantile(R0exp_pois_CI, c(0.025, 0.975), na.rm=TRUE))
        }else{
          count <- count+1
          results[count, ] <- c(s$country[1], sections[j], "ExpPois", NA, NA, NA)
        }
      }else{
        count <- count+1
        results[count, ] <- c(s$country[1], sections[j], "ExpPois", NA, NA, NA)
      }
      
      
      #======== 2.1. MLE Exponential Growth ========#
      
      # define log-likelihood (negative required by mle2)
      ll <- function(I0, r){
        -sum(s$cases*(r*s$week+log(I0))-I0*exp(r*s$week)-lfactorial(s$cases))
      }
      
      # vectors & matrix for storing r, likelihood & CI estimates
      r_est <- vector()
      likel <- vector()
      CI_mle <- data.frame(matrix(ncol=2, nrow=100))
      
      # run MLE from different starting points
      for(m in 1:500){
        
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
      if(is.na(ci_mle[1])+is.na(ci_mle[2])==0 & ci_mle[1]>min_r){
    
        # calculate R by integrating over the GT distribution, g(a)
        R <- 1/(integrate(Mz_ga, r=r_mle, lower=0, upper=Inf, subdivisions=1e8))$value
        R_CIL <- 1/(integrate(Mz_ga, r=CI_mle[1,1], lower=0, upper=Inf, subdivisions=1e8))$value
        R_CIU <- 1/(integrate(Mz_ga, r=CI_mle[1,2], lower=0, upper=Inf, subdivisions=1e8))$value
        
        # store results
        count <- count+1
        results[count, ] <- c(s$country[1], sections[j], "MLE_ExpLin", R, R_CIL, R_CIU)
        
      }else{
        count <- count+1
        results[count, ] <- c(s$country[1], sections[j], "MLE_ExpLin", NA, NA, NA)
      }
      
      
      #======== 3. EpiEstim Method (Renewal Equation) ========#
      
      R_EpiEstim <- EpiEstim::estimate_R(data$cases, method="parametric_si", 
                                         config=make_config(list(mean_si=20/7, std_si=7.4/7,
                                                                 t_start=2, t_end=sections[j])))
      
      # store results
      count <- count+1
      if(exists("R_EpiEstim")==TRUE){
        results[count, ] <- c(data$country[1], sections[j], "EpiEstim", R_EpiEstim$R$`Mean(R)`,
                              R_EpiEstim$R$`Quantile.0.025(R)`, R_EpiEstim$R$`Quantile.0.975(R)`)
        rm(R_EpiEstim)
      }else{
        results[count, ] <- c(data$country[1], sections[j], "EpiEstim", NA, NA, NA)
      }
      

      #======== 4. White & Pagano Method ========#
      
      tryCatch({
        R0_WP <- R0::est.R0.ML(s$cases, GT=GTd, begin=1, end=as.numeric(nrow(s)))
      }, error=function(e){})
      
      # store results
      count <- count+1
      if(exists("R0_WP")==TRUE){
        results[count, ] <- c(s$country[1], sections[j], "WP", R0_WP$R, R0_WP$conf.int)
        rm(R0_WP)
      }else{
        results[count, ] <- c(s$country[1], sections[j], "WP", NA, NA, NA)
      }
      
    }
    
    #=== Output the results
    
    setwd(path)
    results$Country <- paste(data$country[1])
    write.csv(results, file=paste(s$country[1], "Results.csv", sep="_"))
  }
}

