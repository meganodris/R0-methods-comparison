#===== Plotting Zika R0 Estimates from Case Surveillance =====#
#========== Data in Latin America & the Caribbean ==========#

library(ggplot2)
library(cowplot)


#=== Read in case surveillance data
setwd("C:/Users/modrisco/Desktop/Zika R0/Final-for sharing")
c_data <- readRDS("ZikaCases_LatinAmerica&Caribbean.RDS") 


#=== Read in R0 estimates
path <- "Q:/R0 analysis/Results_Observed"
setwd(path)
file_names <- list.files(path=path, pattern="\\.csv$") 
results <- lapply(file_names, read.csv)


# exclude the countries where cases peaked earlier than 6 weeks (2 generation times)
# (Antigua and Barbuda, Nicaragua, Venezuela)
c_data <- c_data[c(2:22, 24:34, 36:37)]



#=== Plot R0 results against case data per country
plots <- list()

for(i in 1:length(c_data)){
  
  # extract case data & R0 estimates per country
  df <- c_data[[i]]
  r0 <- results[[i]]
  r0$method <- factor(r0$method, levels=c("ExpLin", "ExpPois", "EpiEstim", "WP"))
  
  # define peak
  peak <- max.col(matrix(df$cases,nrow=1),"last") 
  
  # define sections from 2+ generation times
  sections <- seq(from=6, to=peak, by=3)
  
  # mid points for plotting
  r0$mpt[which(r0$Nweeks==6)] <- 3
  r0$mpt[which(r0$Nweeks!=6)] <- r0$Nweeks[which(r0$Nweeks!=6)]-1
  
  # case time series
  p1 <- ggplot(df, aes(week, cases))+ ggtitle(df$country[1])+ theme_grey()+
    geom_rect(xmin=-Inf, xmax=peak+0.5, ymin=-Inf, ymax=Inf, fill="lemonchiffon", alpha=0.5)+
    ylab("cases")+ theme(axis.line=element_line(), axis.title.x=element_blank(), axis.text.x=element_blank(),
                         axis.ticks.x=element_blank())+ 
    geom_vline(xintercept=sections+0.5, linetype="dashed", colour="grey")+ 
    geom_vline(xintercept=peak+0.5, linetype="dashed")+ geom_point(colour="blue")
  
  # R0 estimates
  p2 <- ggplot(r0, aes(mpt, R0, ymin=CI_L, ymax=CI_U, group=method))+ theme_grey()+
    geom_rect(xmin=-Inf, xmax=sections[length(sections)]+0.5, ymin=-Inf, ymax=Inf, fill="white")+
    geom_point(aes(colour=method), position=position_dodge(width=2), size=2)+ ylab("R0")+  
    geom_errorbar(aes(colour=method), position=position_dodge(width=2))+ xlab("Weeks")+
    theme(axis.line=element_line(), legend.position="none", plot.margin=unit(c(0,5.5,5.5,5.5), "pt"),
          axis.title.x=element_blank())+ xlim(1,max(df$week))+
    scale_color_manual(values=c("royalblue1", "violetred1", "lawngreen", "orange1"))+
    geom_vline(xintercept=sections+0.5, linetype="dashed", colour="grey")
  
  # bind the 2 plots & store
  p <- plot_grid(p1,p2, nrow=2, align="v", rel_widths=c(1.3,1))
  plots[[i]] <- p
  
  # expect quite a few warnings about missing values where it was not possible to estimate R0
}

# check the plots
plots[[6]]
