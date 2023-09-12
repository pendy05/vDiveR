#' Time Distribution of Sequences Plot
#'
#' This function plots the time distribution of provided sequences in the form of bar plot with 'Month' as x-axis and
#' 'Number of Sequences' as y-axis. Aside from the plot, this function also returns a dataframe with 2 columns: 'Date' and 'Number of sequences'.
#' The input dataframe of this function is obtainable from metadata_extraction(), with NCBI Protein / GISAID EpiCoV FASTA file as input.
#'
#' @param metadata a dataframe with 3 columns, 'ID', 'country', and 'date'
#' @param base_size word size in plot
#' @param date_format date format of the input dataframe
#' @param date_break date break for the scale_x_date
#' @param scale plot counts or log scale the data
#' @param only_plot logical, return only plot or dataframe info as well, default FALSE
#'
#' @return A single plot or a list with 2 elements (a plot followed by a dataframe, default)
#' @examples time_plot <- plot_time(metadata)$plot
#' @examples time_df <- plot_time(metadata)$df
#' @importFrom ggplot2 scale_x_date
#' @importFrom ggplot2 geom_bar
#' @importFrom scales date_format
#' @export
plot_time <- function(metadata, date_format = "%Y-%m-%d", base_size=8,
                      date_break = "2 month", scale = "count", only_plot = F){
    Month <- NULL
    metadata$Month <- as.Date(cut(as.Date(metadata$Date, format = date_format),
                                  breaks = "month"))

    p <- ggplot(data = metadata, aes(x = Month)) + geom_bar() +
        ylab('Number of Sequences') +
        scale_x_date(date_breaks = date_break, labels = date_format("%Y-%b"))  +
        theme_classic(base_size = base_size) +
        theme(axis.text.x = element_text(angle = 90, vjust = .5))

    if (scale == 'log') {
        p <- p + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                 labels = trans_format("log10", math_format(10^.x)))
    }

    if (only_plot) {
        return(p)
    } else {
        temporal <- metadata[,c('country', 'date')]
        temporal$count <- 1
        temporal$date <- as.Date(temporal$date, format = date_format)
        temporal <- aggregate(temporal$count, by=list(temporal$date), sum)
        colnames(temporal) <- c('Date', 'Number of Sequences')

        return(list(plot=p,df=temporal))
    }
}

