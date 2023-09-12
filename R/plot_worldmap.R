#' Geographical Distribution of Sequences Plot
#'
#' This function plots a worldmap and color the affected geographical region(s)
#' from light (lower) to dark (higher), depends on the cumulative number of sequences.
#' Aside from the plot, this function also returns a dataframe with 2 columns: 'Country' and 'Number of Sequences'.
#' The input dataframe of this function is obtainable from metadata_extraction(), with NCBI
#' Protein / GISAID EpiCoV FASTA file as input.

#' @param meta a dataframe with 3 columns, 'ID', 'country', and 'date'
#' @param base_size word size in plot
#' @param only_plot logical, return only plot or dataframe info as well, default FALSE
#'
#' @return A single plot or a list with 2 elements (a plot followed by a dataframe)
#' @examples geographical_plot <- plot_worldmap(metadata)$plot
#' @examples geographical_df <- plot_worldmap(metadata)$df
#' @importFrom ggplot2 geom_polygon scale_fill_gradient map_data
#' @importFrom dplyr left_join
#' @export
plot_worldmap <- function(meta, base_size=8, only_plot = F){
    long <- lat <- group <- count <- NULL

    meta <- refineCountry(meta)
    # meta$country[meta$country == "DRC"] = "Democratic Republic of the Congo"
    # meta$country[meta$country == "NewCaledonia"] = "New Caledonia"
    # meta$country[meta$country == "Northern Ireland"] = "New Caledonia"
    # meta$country[meta$country %in% c("England","Scotland","Wales")] = "UK"

    countrylist <- data.frame(table(meta$Country))
    colnames(countrylist) <- c('region','count')

    world_map <- map_data("world")
    p <- ggplot(world_map, aes(x = long, y = lat, group = group)) +
        geom_polygon(fill="lightgray", colour = "#888888")
    pathogens.map <- left_join(countrylist, world_map, by = "region")

    p <- p + geom_polygon(data = pathogens.map, aes(fill = count), color = "#888888") +
        scale_fill_gradient(low = "#FFFFFF", high = "#E63F00", name = 'Number of Sequences') +
        theme(plot.background = element_rect(fill = "transparent", colour = NA),
              panel.border = element_blank(), panel.grid = element_blank(),
              axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
              legend.position = "right") +
        theme_classic(base_size=base_size)

    colnames(countrylist) <- c('Country','Number of Sequences')

    if (only_plot) {
        return(p)
    } else {
        return(list(plot=p, df=countrylist))
    }

}


refineCountry <- function(metatable){
    metatable$Country[metatable$Country == "DRC"] = "Democratic Republic of the Congo"
    metatable$Country[metatable$Country == "NewCaledonia"] = "New Caledonia"
    metatable$Country[metatable$Country == "Northern Ireland"] = "New Caledonia"
    metatable$Country[metatable$Country %in% c("England","Scotland","Wales")] = "UK"
    metatable$Country[metatable$Country %in% c("Shangahi", "Xinjiang","Sichuan", "Guangdong","Shannxi", "Chongqing", "Inner_Mongolia","Shenzhen", "Wuhan",
                                               "Fujian", "Inner Mongolia", "Tianjing", "Hebei","Jiangsu", "Shandong", "Zhejiang",
                                               "Liaoning","Shanxi", "Henan", "Chongqin", "Yunnan", "Beijing","Heilongjiang",
                                               "Hunan", "Guangxi","Ningxia","Jilin","Tibet","Hainan", "Macao",
                                               "Jiangxi","Qinghai", "Hubei", "Gansu", "Anhui" ,"Guizhou" )] = "China"
    return(metatable)
}
