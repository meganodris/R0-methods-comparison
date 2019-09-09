#===== Function to estimate R at different stages of the epidemic growth period ====#

R_estimates <- function(data, path){
  
  data <- as.data.frame(data)
  
  # define the moment generating function of the GT distribution, g(a)
  Mz_ga <- function(t, r){
    dgamma(t, shape=(mean_GT^2/variance_GT), scale=(variance_GT/mean_GT))*exp(-t*r)
  }
  
  # define a fixed mean and variance of the generation time, GT
  mean_GT <- 20/7
  sd_GT <- 7.4/7
  variance_GT <- sd_GT^2
  GTd <- R0::generation.time("gamma", c(mean_GT, sd_GT))
  
  # dataframe to store results
  results <- data.frame(Country=NA, Nweeks=NA, method=NA, R0=NA, CI_L=NA, CI_U=NA)
  
  # define peak - furthest point of highest reported cases
  peak <- max.col(matrix(data$cases,nrow=1),"last") 
  
  # if the peak is >6 weeks after 1st case report, run analysis
  if(peak>=6){
    
    # define sections: from 2+ generation intervals
    sections <- seq(from=6, to=peak, by=3)
    
    count <- 0
    
    # Loop to fit each method to each section 
    for(j in 1:length(sections)){
    
      # extract section for fitting
      s <- data[1:sections[j], ]
      
      #======== 1. Linear Exponential Growth ========#
      
      # log transform cases
      s$cases[which(s$cases==0)] <- 1
      s$logcases <- log(s$cases)
      
      # define growth rate, r
      m <- lm(logcases~ week, data=s)
      r <- coef(m)[2]
      
      # calculate R by integrating over the GT distribution, g(a)
      int_ga <- integrate(Mz_ga, r=r, lower=0, upper=Inf, subdivisions=100000000)
      R <- 1/int_ga$value
      
      if(R>1.00000000000206){ # equivalent of r=0
        
        #=== calculate 95% CI for R
        
        # stanard error of r
        se_r <- coef(summary(m))[, "Std. Error"][2]
        
        # bootstrap sample from t-distribution
        boot <- vector()
        for (k in 1:10000){
          # this can take a while!!
          rdf <- rt(1, df=(m$df.residual))
          r_new <- abs(r + se_r*rdf) 
          int_ga_new <- integrate(Mz_ga, r=r_new, lower=0, upper=Inf, subdivisions=100000000)
          boot[k] <- 1/int_ga_new$value
          
        }
        
        
        # store results
        count <- count+1
        results[count, ] <- c(s$country[1], sections[j], "ExpLin", R, quantile(boot, c(0.025,0.975), na.rm=TRUE))
      }else{
        count <- count+1
        results[count, ] <- c(s$country[1], sections[j], "ExpLin", NA, NA, NA)
      }
      
      #======== 2. Poisson Exponential Growth ========#
      
      # define growth rate, r
      mP <- glm(cases~ week, data=s, family="poisson")
      r_p <- coef(mP)[2]
      
      # calculate R by integrating over the GT distribution, g(a)
      int_ga <- integrate(Mz_ga, r=r_p, lower=0, upper=Inf, subdivisions=100000000)
      R_p <- 1/int_ga$value
      
      if(R_p>1.00000000000206){
        
        #=== calculate 95% CI for R
        
        # standard error of r
        se_r_p <- coef(summary(mP))[, "Std. Error"][2]
        
        # bootstrap sample from t-distribution
        R0exp_pois_CI <- vector()
        for (k in 1:10000){
          # simulate exponential growth rate from t-distrbution
          rdf <- rt(1, df=(mP$df.residual))
          r_new <- abs(r_p + se_r_p*rdf)
          int_ga_new <- integrate(Mz_ga, r=r_new, lower=0, upper=Inf, subdivisions=100000000)
          R0exp_pois_CI[k] <- 1/int_ga_new$value
          
        }
        
        # store results
        count <- count+1
        results[count, ] <- c(s$country[1], Nweeks=sections[j], "ExpPois", R_p, quantile(R0exp_pois_CI, c(0.025, 0.975), na.rm=TRUE)) 
      }else{
        count <- count+1
        results[count, ] <- c(s$country[1], sections[j], "ExpPois", NA, NA, NA)
      }
      
      #======== 3. EpiEstim Method (Renewal Equation) ========#
      
      # define start and end of sections
      startG <- 2
      endG <- sections[j]
      
      R_EpiEstim <- EpiEstim::EstimateR(data$cases, T.Start=startG, T.End=endG, method="ParametricSI",
                                        Mean.SI=mean_GT, Std.SI=sd_GT)
      
      # store results
      count <- count+1
      if(exists("R_EpiEstim")==TRUE){
        results[count, ] <- c(data$country[1], sections[j], "EpiEstim", R_EpiEstim$R$`Mean(R)`,
                              R_EpiEstim$R$`Quantile.0.025(R)`, R_EpiEstim$R$`Quantile.0.975(R)`)
      }else{
        results[count, ] <- c(data$country[1], sections[j], "EpiEstim", NA, NA, NA)
      }
      rm(R_EpiEstim)
      
      #======== 4. White & Pagano Method ========#
      
      tryCatch({
        R0_WP <- R0::est.R0.ML(s$cases, GT=GTd)
      }, error=function(e){})
      
      # store results
      count <- count+1
      if(exists("R0_WP")==TRUE){
        results[count, ] <- c(s$country[1], sections[j], "WP", R0_WP$R, R0_WP$conf.int)
      }else{
        results[count, ] <- c(s$country[1], sections[j], "WP", NA, NA, NA)
      }
      rm(R0_WP)
    }
    
    #=== Output the results
    
    setwd(path)
    results$Country <- paste(data$country[1])
    write.csv(results, file=paste(s$country[1], "Results.csv", sep="_"))
  }
}
