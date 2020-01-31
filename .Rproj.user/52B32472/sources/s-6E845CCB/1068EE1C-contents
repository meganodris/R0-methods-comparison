
# root mean squared error
rmse <- function(error){
  round(sqrt(mean(error^2)), digits=2)
}

# mean absolute error
mae <- function(error){
  round(mean(abs(error)), digits=2)
}

# summarise performance
performance_metrics <- function(trueR0, estimates, max_weeks){
  
  # unlist & bind R0 estimates
  estimates <- do.call('rbind', estimates)
  estimates$sim <- substr(estimates$Country, 5, 7)
  
  # merge with actual/true R0 values
  trueR0$sim <- as.character(trueR0$sim)
  estimates <- left_join(estimates, trueR0, by='sim')
  colnames(estimates)[c(4,8)] <- c('R0', 'R0_true')
  
  # subset analyses to the early epidemic growth phase, e.g. 15 weeks
  estimates <- subset(estimates, estimates$Nweeks<=max_weeks)
  
  # remove sections of simulations with any missing R0 estimates (for systematic comparison across methods)
  missing <- as.data.frame(cbind(sim=estimates$sim[which(is.na(estimates$R0))],
                                 Nweeks=estimates$Nweeks[which(is.na(estimates$R0))]))
  missing <- unique(missing[c("sim", "Nweeks")])
  for(m in 1:nrow(missing)){
    estimates <- estimates[!(estimates$sim==missing$sim[m] & estimates$Nweeks==missing$Nweeks[m]), ]
  }
  
  #=== Calculate performance metrics
  
  estimates$bias <- estimates$R0 - estimates$R0_true # bias
  
  estimates$ciW <- estimates$CI_U - estimates$CI_L # width of CI / uncertainty
  
  # coverage by 95% CI's (1=yes, 0=no)
  estimates$covYN <- NA
  for(c in 1:nrow(estimates)){
    if(estimates$CI_U[c]>=estimates$R0_true[c] & estimates$CI_L[c]<=estimates$R0_true[c]){
      estimates$covYN[c] <- 1
    }else{
      estimates$covYN[c] <- 0
    }
  }

  # summarize method performance at each stage
  metrics <- data.frame(Nweeks=NA, method=NA, r2=NA, slope=NA, pcc=NA, scc=NA, rmse=NA, mae=NA,
                        CIwidth=NA, cov=NA, abs_bias=NA, freq_O=NA, freq_U=NA, mag_O=NA, mag_U=NA)
  
  indx <- 0 # index for storing results
  for(t in unique(estimates$Nweeks)){ # subset data for time t
    
    ests_t <- estimates[estimates$Nweeks==t, ] 
    
    for(method in unique(ests_t$method)){ # subset for each method
      
      ests_tm <- ests_t[ests_t$method==method, ]
      
      # linear regression for slope and r-squared
      lin_mod <- lm(R0 ~ R0_true, data=ests_tm)
      r2 <- format(round(summary(lin_mod)$r.squared, digits=2))
      m <- format(round(coef(lin_mod)[2], digits=2))
      
      # pearsons & spearmans correlation coefficients
      pears <- format(round(cor(ests_tm$R0_true, ests_tm$R0, method="pearson"), digits=2))
      spear <- format(round(cor(ests_tm$R0_true, ests_tm$R0, method="spearman"), digits=2))
      
      # errors 
      errors <- ests_tm$R0_true - ests_tm$R0
      
      # store results
      indx <- indx+1
      metrics[indx, ] <- c(paste(t), paste(method), r2, m, pears, spear, rmse(errors), mae(errors),
                           mean(ests_tm$ciW), (sum(ests_tm$covYN)/nrow(ests_tm)), mean(abs(ests_tm$bias)),
                           (nrow(ests_tm[ests_tm$bias>0, ])/nrow(ests_tm)), (nrow(ests_tm[ests_tm$bias<0, ])/nrow(ests_tm)),
                           mean(ests_tm$bias[ests_tm$bias>0]), mean(ests_tm$bias[ests_tm$bias<0]))
    }
  }
  return(list(metrics_summ=metrics, metrics_indiv=estimates))
}


# function to return a vector of metric names
get_metrics <- function(){
  
  m_names <- c('Rsquared', 'Slope', 'PCC', 'SCC', 'RMSE', 'MAE', 'Uncertainty', 'Coverage',
               'Bias', 'Freq_overest', 'Freq_underest', 'Bias_overest', 'Bias_underest')
  return(m_names)
}

# metrics summary plot
summary_plot <- function(metrics_summ, include){
  
  # reformat & keep only metrics of interest 
  colnames(metrics_summ)[3:15] <- get_metrics()
  met <- tidyr::gather(metrics_summ, key='metrics', value='value', Rsquared:Bias_underest)
  met <- met[met$metrics %in% include, ]
  
  # plot
  met$Nweeks <- factor(met$Nweeks, levels=c(6,9,12,15))
  met$method <- factor(met$method, levels=c('ExpLin', 'ExpPois', 'MLE_ExpLin', 'EpiEstim', 'WP', 'WT'))
  plotM <- ggplot(met, aes(Nweeks, as.numeric(value), colour=method, shape=method))+ geom_point()+ 
    facet_wrap(~metrics, scales="free_y")+ ylab('Performance')+ xlab('Number of weeks')+
    scale_colour_manual(values=c("royalblue1", "violetred1", "lawngreen", "orange1", "turquoise1", "purple"))
  
  return(plotM)
}


# bias distribution plot: estimated R0 - actual R0
bias_plot <- function(metrics_indiv){
  
  metrics_indiv$method <- factor(metrics_indiv$method, levels=c('ExpLin', 'ExpPois', 'MLE_ExpLin', 'EpiEstim', 'WP', 'WT'))
  plotB <- ggplot(metrics_indiv)+ ylab('')+ scale_x_continuous(expand=c(0, 0))+ theme_minimal()+
    geom_density_ridges(aes(x=bias, y=method, fill=method, col=method), scale=0.95, alpha=0.5, bandwidth=0.25)+
    facet_grid(~Nweeks)+ geom_vline(xintercept=0)+ xlab('bias')+ 
    scale_fill_manual(values=c("royalblue1", "violetred1", "lawngreen", "orange1", "turquoise1", "purple"))+
    scale_colour_manual(values=c("royalblue1", "violetred1", "lawngreen", "orange1", "turquoise1", "purple"))+
    theme(legend.position='right', axis.text.y=element_blank(), axis.ticks.y=element_blank())+ 
    scale_y_discrete(expand=expand_scale(mult=c(0.01, 0.01)), limits=rev(levels(metrics_indiv$method)))
  
  return(plotB)  
}
