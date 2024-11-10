#' Geographical Distribution of Sequences Plot
#'
#' This function plots a world map and color the affected geographical region(s)
#' from light (lower) to dark (higher), depends on the cumulative number of sequences.
#' Aside from the plot, this function also returns a dataframe with 2 columns: 'Region' and 'Number of Sequences'.
#' The input dataframe of this function is obtainable from metadata_extraction(), with NCBI
#' Protein / GISAID (EpiFlu/EpiCoV/EpiPox/EpiArbo) FASTA file as input.

#' @param metadata a dataframe with 3 columns, 'ID', 'region', and 'date'
#' @param base_size word size in plot
#'
#' @return A list with 2 elements (a plot followed by a dataframe)
#' @examples geographical_plot <- plot_world_map(metadata)$plot
#' @examples geographical_df <- plot_world_map(metadata)$df
#' @importFrom ggplot2 geom_polygon scale_fill_gradient map_data
#' @importFrom dplyr left_join %>% group_by ungroup slice 
#' @importFrom stringr str_to_title
#' @export
plot_world_map <- function(metadata, base_size=8){
    long <- lat <- group <- count <- city_ascii <- region.y <- region <- ID <- date <- NULL
    if (nrow(metadata) < 1){
        error_msg <- paste("No records found in the metadata dataframe.")
        return(list(plot = NULL, df = error_msg))
    }
    if (!all(c("ID", "region", "date") %in% colnames(metadata))) {
        stop("The metadata dataframe must contain 'ID', 'region', and 'date' columns.")
    }

  
    #============= data preparation section - map city to region ========================#
    build_in_path <- system.file("extdata", "city_mapper.csv", package = "vDiveR")
    city2region <- utils::read.csv(build_in_path, stringsAsFactors = FALSE) 
    city2region_unique <- city2region %>%
        dplyr::group_by(city_ascii) %>%
        dplyr::slice(1) %>%  # if more than 1 match, keep the first occurrence
        dplyr::ungroup()

    metadata$region <- tolower(metadata$region)
    city2region_unique$city_ascii <- tolower(city2region_unique$city_ascii)

    # Replace the 'region' column in metadata with the corresponding 'region' value of the matched 'city_ascii' from city_mapper
    metadata <- metadata %>%
        rowwise() %>%
        mutate(
            region = if_else(
                region %in% city2region_unique$city_ascii,  # Check if region matches any city in city_mapper
                city2region_unique$region[match(region, city2region_unique$city_ascii)],  # Replace with corresponding region if city matches
                region  # Otherwise, leave it as is
            )
        ) %>%
        ungroup() %>%
        as.data.frame()

    #============ data preparation section - map data region to ggplot2 region ===========#
    world_map <- ggplot2::map_data("world")
    world_region_list <- unique(world_map$region) # Get unique regions from the world map data

    # Apply the matching function to the 'region' column
    matched_regions <- match_region_to_target(metadata$region, world_region_list)
    # Update 'region' column with matched regions where available
    metadata$region <- ifelse(is.na(matched_regions), metadata$region, matched_regions)
    # Check unique regions in the updated 'region' column
    unique_regions <- unique(metadata$region)
    # Identify any regions not in the target ggplot world map region list
    regions_not_in_target <- setdiff(unique_regions,  world_region_list)
    # Convert 'region' to factor with levels from target list
    metadata$region <- factor(metadata$region, levels =  world_region_list)

    if (length(regions_not_in_target) > 0) {
         warning_info <- paste('The following regions could not be matched and need manual review:\n',as.character(regions_not_in_target), 
         '.\nPlease ensure the provided regions are part of the ggplot2 world map regions. You may check the valid regions via:\n unique(ggplot2::map_data("world")$region)', collapse = ", ")
        warning(warning_info)
    }


    #================= plotting and tabulating section =====================#
    region_list <- data.frame(table(metadata$region))
    colnames(region_list) <- c('region','count')

    p <- ggplot(world_map, aes(x = long, y = lat, group = group)) +
        geom_polygon(fill="lightgray", colour = "#888888")
    pathogens.map <- dplyr::left_join(region_list, world_map, by = "region")

    p <- p + geom_polygon(data = pathogens.map, aes(fill = count), color = "#888888") +
        scale_fill_gradient(low = "#FFFFFF", high = "#E63F00", name = 'Number of Sequences') +
        theme(plot.background = element_rect(fill = "transparent", colour = NA),
              panel.border = element_blank(), panel.grid = element_blank(),
              axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
              legend.position = "right") +
        theme_classic(base_size=base_size)

    
    colnames(region_list) <- c('Region','Number of Sequences')

    return(list(plot=p, df=region_list))

}

#' Map data region to ggplot2 region
#'
#' This function maps the provided data region to ggplot2 region.
#' @param region_column a character vector of data region
#' @param region_list a character vector of ggplot2 region
#' @param max_dist maximum distance for string matching
#' @importFrom stringdist stringdist
match_region_to_target <- function(region_column, region_list, max_dist = 3) {
    matched_regions <- sapply(region_column, function(x) {
        # Compute string distances
        distances <- stringdist::stringdist(tolower(x), tolower(region_list), method = "lv")
        min_index <- which.min(distances)
        min_dist <- distances[min_index]

        # Check if the minimum distance is within the acceptable threshold
        if (min_dist <= max_dist) {
            # Return the target list element with matching capitalization
            region_list[min_index]
        } else {
            # Return NA if no close match is found
            NA
        }
    })
    return(matched_regions)
}