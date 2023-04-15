#' plot the geographical location via metadata
#'
#' This function plot the cumulative sequence number in geographical via metadata
#'
#'
#' @param meta sequence metadata from 'metadataExtraction' function or user provided, the format of metadata can refer to built-in data, "metadata".
#' @return  A plot
#' @examples filepath <- system.file('inst/extdata','prot_NCBIProtein.fasta', package = 'vDiveR')
#' @examples meta <- metadataExtraction(filepath, 'ncbi')
#' @examples plot_worldmap(meta)
#' @examples plot_worldmap(metadata)
#' @importFrom ggplot2 geom_polygon scale_fill_gradient map_data
#' @importFrom dplyr left_join
#' @export
plot_worldmap <- function(meta){
    meta$country[meta$country == "DRC"] = "Democratic Republic of the Congo"
    meta$country[meta$country == "NewCaledonia"] = "New Caledonia"
    meta$country[meta$country == "Northern Ireland"] = "New Caledonia"
    meta$country[meta$country %in% c("England","Scotland","Wales")] = "UK"
    countrylist <- data.frame(table(meta$country))
    colnames(countrylist) <- c('region','count')

    world_map <- map_data("world")
    p <- ggplot(world_map, aes(x = long, y = lat, group = group)) + geom_polygon(fill="lightgray", colour = "#888888")
    pathogens.map <- left_join(countrylist, world_map, by = "region")
    p <- p + geom_polygon(data = pathogens.map, aes(fill = count), color = "#888888") +
        scale_fill_gradient(low = "#FFFFFF", high = "#E63F00", name = 'Number of sequences') +
        theme(plot.background = element_rect(fill = "transparent", colour = NA),
              panel.border = element_blank(), panel.grid = element_blank(),
              axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
              legend.position = "right") + theme_classic()
    return(p)
}
