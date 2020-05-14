#=== Function to estimate time-constant R0 with each method at sequential ===#
#=== time points in the epidemic growth phase and output results as csv's ===#

fit_R0_seq <- function(data, mean_GT, sd_GT, GTd, GT_week){
  
  # list to store results
  store <- list()
  
  # method names
  methods <- c("EG_Lin", "EG_P", "EG_MLE", "EpiEstim", "WP", "WT")
  
  # define peak - max time point of highest reported cases
  peak <- which.max(data$cases)[length(which.max(data$cases))] 
  
  if(peak>=6){
    
    # define sections: from 2 generation times (approximated in terms of weeks) on, up to peak
    sections <- seq(from=GT_week*2, to=peak, by=GT_week)
    
    # Loop to fit each method to each section with increasing number of time points
    for(t in 1:length(sections)){
      
      # extract section of the epidemic curve for fitting
      s_t <- data[1:sections[t], ]
      
      # fit each method
      eglin <- EG_Lin(data=s_t, mean_GT=mean_GT, sd_GT=sd_GT)
      egp <- EG_P(data=s_t, mean_GT=mean_GT, sd_GT=sd_GT)
      egmle <- EG_MLE(data=s_t, mean_GT=mean_GT, sd_GT=sd_GT)
      epiest <- epiestim(data=s_t, mean_GT=mean_GT, sd_GT=sd_GT)
      wp <- WP(data=s_t, GTd=GTd)
      wt <- WT(data=s_t, GTd=GTd)
      
      # store results
      store[[t]] <- cbind(paste(data$country[1]), sections[t], methods, rbind(eglin, egp, egmle, epiest, wp, wt), peak)
      colnames(store[[t]]) <- c("Country", "Nweeks", "method", "R0", "CI_L", "CI_U", "peak")
    }
    
    # bind results by country and output as csv
    results_all <- do.call("rbind", store)
    return(results_all)
  }
}
