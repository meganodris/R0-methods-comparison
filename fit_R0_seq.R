#=== Function to estimate time-constant R0 with each method at sequential ===#
#=== time points in the epidemic growth phase and output results as csv's ===#

fit_R0_seq <- function(data, mean_GT, sd_GT, GTd){
  
  # list to store results
  store <- list()
  
  # method names
  methods <- c("ExpLin", "ExpPois", "MLE_ExpLin", "EpiEstim", "WP", "WT")
  
  # define peak - max time point of highest reported cases
  peak <- which.max(data$cases)[length(which.max(data$cases))] 
  
  if(peak>=6){
    
    # define sections: from 6 weeks, by 3, up to peak (i.e. 6,9,12,15...)
    sections <- seq(from=6, to=peak, by=3)
    
    # Loop to fit each method to each section with increasing number of time points
    for(t in 1:length(sections)){
      
      # extract section of the epidemic curve for fitting
      s_t <- data[1:sections[t], ]
      
      # fit each method
      exp_lin <- EG_Lin(data=s_t, mean_GT=mean_GT, sd_GT=sd_GT)
      exp_pois <- EG_P(data=s_t, mean_GT=mean_GT, sd_GT=sd_GT)
      exp_mle <- EG_MLE(data=s_t, mean_GT=mean_GT, sd_GT=sd_GT)
      epiest <- epiestim(data=s_t, mean_GT=mean_GT, sd_GT=sd_GT)
      wp <- WP(data=s_t, GTd=GTd)
      wt <- WT(data=s_t, GTd=GTd)
      
      # store results
      store[[t]] <- cbind(paste(data$country[1]), sections[t], methods, rbind(exp_lin, exp_pois, exp_mle, epiest, wp, wt))
      colnames(store[[t]]) <- c("Country", "Nweeks", "method", "R0", "CI_L", "CI_U")
    }
    
    # bind results by country and output as csv
    results_all <- do.call("rbind", store)
    return(results_all)
  }
}
