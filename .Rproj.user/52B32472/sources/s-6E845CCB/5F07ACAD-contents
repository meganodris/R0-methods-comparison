#=== Function to fit each method to increasing number of time points ===#
#=== in the epidemic growth phase and output results as csv's ===#

fit_R0_seq <- function(data){
  
  # list to store results
  store <- list()
  
  # method names
  methods <- c("ExpLin", "ExpPois", "MLE_ExpLin", "EpiEstim", "WP")
  
  # define peak - max time point of highest reported cases
  peak <- max.col(matrix(data$cases,nrow=1),"last") 
  
  if(peak>=6){
    
    # define sections: from 6 weeks, by 3, up to peak (i.e. 6,9,12,15...)
    sections <- seq(from=6, to=peak, by=3)
    
    # Loop to fit each method to each section with increasing number of time points
    for(j in 1:length(sections)){
      
      # extract section of the epidemic curve for fitting
      s_j <- data[1:sections[j], ]
      
      # fit each method
      exp_lin <- EG_Lin(data=s_j)
      exp_pois <- EG_P(data=s_j)
      exp_mle <- EG_MLE(data=s_j)
      epiest <- epiestim(data=s_j, mean_si=mean_GT, std_si=sd_GT)
      wp <- WP(data=s_j, GTd=GTd)
      
      # store results
      store[[j]] <- cbind(paste(data$country[1]), sections[j], methods, rbind(exp_lin, exp_pois, exp_mle, epiest, wp))
      colnames(store[[j]]) <- c("Country", "Nweeks", "method", "R0", "CI_L", "CI_U")
    }
    
    # bind results by country and output as csv
    results_all <- do.call("rbind", store)
    write.csv(results_all, paste(data$country[1], "Results.csv", sep="_"))
  }
}
