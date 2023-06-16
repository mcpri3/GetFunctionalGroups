#' Function that read and clean GBIF occurrence file 
#'
#' @param ffile name of the GBIF occurrence file to read
#' @param uncertainty resolution in meters to filter data from, for ex., 1000 means that coordinates with uncertainty higher than 1000m are removed
#' @param cols.tokeep names of columns to keep from the GBIF file 
#' @param year.min minimum year to filter data from, for ex., 2000 means that occurences before year 2000 are removed 
#'
#' @return
#' @export
#'
#' @examples
clean_occurGBIF = function(ffile, uncertainty = 1000, year.min = 2000,
                           cols.tokeep = c('species', 'eventDate','day','month','year','X','Y', 'coordinateUncertaintyInMeters')) {
  
  occur <- readr::read_delim(here::here(paste0("data/raw-data/GBIF/Chordata/",ffile)), 
                             delim = ";", escape_double = FALSE, trim_ws = TRUE)
  occur <- occur[occur$taxonRank == 'SPECIES',]
  occur <- occur[occur$coordinateUncertaintyInMeters <= uncertainty,]
  occur <- occur[occur$year >= year.min,]
  
  if (nrow(occur) > 0) {
  if (sum(cols.tokeep %in% colnames(occur)) == length(cols.tokeep)) {
    occur = occur[occur$distance_to_iucn == 0,]
    occur <- occur[, cols.tokeep]
    
  } else {
    
    xy <- reshape2::colsplit(occur$geometry, "\\|", c('X','Y'))
    xy <- sf::st_as_sf(xy, coords = c('X','Y'), crs = 4258)
    xy <- sf::st_transform(xy, crs = 4326)
    occur <- cbind(occur, sf::st_coordinates(xy))
    occur <- occur[, cols.tokeep]
  }
  }
  return(occur)
}
