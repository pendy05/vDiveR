#plot correlation between entropy and total variant incidence
plot_correlation_of_entropy <- function(df, host = 1 , alpha = 1/3, size = 3, ylabel = "k-mer entropy (bits)\n", xlabel = "\nTotal variants (%)", ymax = ceiling(max(df$entropy)) ,ybreak=0.5){
    plot2<-ggplot(df)+geom_point(mapping = aes(x=totalVariants.incidence,y=entropy),alpha=alpha,size=size)+
        labs(y = ylabel,x= xlabel)+
        scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20))+
        scale_y_continuous(limits = c(0, ymax), breaks = seq(0, ymax, ybreak))+
        theme_classic()+
        theme(
            panel.border = element_rect(colour = "#000000", fill=NA, size=1)
        )

    if (host == 1){ #single host
        #plot the scatter plot with density
        plot2
    }else{ #multiple host
        plot2+facet_grid_sc(rows = vars(df$host),space = "free",switch = "x")
    }

}
