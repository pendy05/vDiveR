#' plot the time of collection via metadata
#'
#' This function plot the time of collection via metadata
#'
#'
#' @param meta sequence metadata from 'metadata_extraction' function or user provided, the format of metadata can refer to built-in data, "metadata".
#' @return  A plot
#' @examples filepath <- system.file('inst/extdata','prot_NCBIProtein.fasta', package = 'vDiveR')
#' @examples meta <- metadata_extraction(filepath, 'ncbi')
#' @examples plot_time(meta)
#' @examples plot_time(metadata)
#' @importFrom ggplot2 scale_x_date
#' @importFrom ggplot2 geom_bar
#' @importFrom scales date_format
#' @export
plot_time <- function(meta){
    meta$Month <- as.Date(cut(as.Date(meta$date, format = "%Y-%m-%d"), breaks = "month"))
    p <- ggplot(data = meta, aes(x = Month)) + geom_bar() + ylab('Number of protein sequence records') +
        scale_x_date(date_breaks = "2 month", labels = date_format("%Y-%b"))  + theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = .5))
    return(p)
}
