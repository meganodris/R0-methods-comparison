#=== Function to plot R0 estimates against case data
library(cowplot)
library(ggplot2)

plotting_R0 <- function(case_data, R0_ests){
  
  # define peak
  peak <- which.max(case_data$cases)[length(which.max(case_data$cases))] 
  
  # define sections from 2+ generation times
  sections <- seq(from=6, to=peak, by=3)
  
  # define mid-points for plotting
  R0_ests$mpt[which(R0_ests$Nweeks==6)] <- 3
  R0_ests$mpt[which(R0_ests$Nweeks!=6)] <- R0_ests$Nweeks[which(R0_ests$Nweeks!=6)]-1
  
  # order methods
  R0_ests$method <- factor(R0_ests$method, levels=c("EG_Lin", "EG_P", "EG_MLE", "EpiEstim", "WP", 'WT', "BR"))
  
  # plot of case numbers over time
  p1 <- ggplot(case_data, aes(week, cases))+ ggtitle(case_data$country[1])+ theme_grey()+ 
    geom_rect(xmin=-Inf, xmax=peak+0.5, ymin=-Inf, ymax=Inf, fill="lemonchiffon", alpha=0.5)+
    theme(axis.line=element_line(), axis.title.x=element_blank(),axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+ 
    geom_bar(stat="identity", fill="blue", alpha=0.6, width=1)+ xlab("cases")+ 
    geom_vline(xintercept=sections+0.5, linetype="dashed", col="grey50")+ 
    geom_vline(xintercept=peak+0.5)+ scale_y_continuous(expand=c(0,1))
  
  # plot of R0 estimates over time
  p2 <- ggplot(R0_ests, aes(mpt, R0, ymin=CI_L, ymax=CI_U, group=method))+ theme_grey()+
    geom_rect(xmin=-Inf, xmax=sections[length(sections)]+0.5, ymin=-Inf, ymax=Inf, fill="white")+
    geom_point(aes(colour=method), position=position_dodge(width=2), size=2)+ ylab("R0")+  
    geom_errorbar(aes(colour=method), position=position_dodge(width=2))+ xlab("Weeks")+
    theme(axis.line=element_line(), legend.position="none", plot.margin=unit(c(0,5.5,5.5,5.5), "pt"), 
          axis.title.x=element_blank())+ xlim(1,max(case_data$week))+
    scale_color_manual(values=c("royalblue1", "violetred1", "lawngreen", "orange1", "turquoise1", 'purple', 'grey60'),
                       labels=c("EG_Lin", "EG_P", "EG_MLE", "EpiEstim", "WP", 'WT', 'BR'))+
    geom_vline(xintercept=sections+0.5, linetype="dashed", colour="grey")
  
  # bind the 2 plots & return
  final_plot <- plot_grid(p1,p2, nrow=2, align="v", rel_widths=c(1.3,1))
  return(final_plot)
}
