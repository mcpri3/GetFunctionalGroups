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
clean_occurGBIF = function(ffile, code = NULL, uncertainty = 1000, year.min = 2000,
                           cols.tokeep = c('species', 'eventDate','day','month','year','X','Y', 'coordinateUncertaintyInMeters')) {

  occur <- readr::read_delim(here::here(paste0("data/raw-data/GBIF/Chordata/",ffile)), 
                             delim = ";", escape_double = FALSE, trim_ws = TRUE)
  if (ncol(occur) == 1 ) {
    occur <- readr::read_csv(here::here(paste0("data/raw-data/GBIF/Chordata/",ffile)))
    }
  occur <- occur[occur$taxonRank == 'SPECIES',]
  occur <- occur[occur$coordinateUncertaintyInMeters <= uncertainty,]
  occur <- occur[occur$year >= year.min,]
  
  if (nrow(occur) > 0) {
    
  # if (sum(cols.tokeep %in% colnames(occur)) == length(cols.tokeep)) {
  #   occur = occur[occur$distance_to_iucn == 0,]
  #   occur <- occur[, cols.tokeep]
  # } else {
  #   
  #   xy <- reshape2::colsplit(occur$geometry, "\\|", c('X','Y'))
  #   xy <- sf::st_as_sf(xy, coords = c('X','Y'), crs = 3035)
  #   xy <- sf::st_transform(xy, crs = 4326)
  #   occur <- cbind(occur, sf::st_coordinates(xy))
  #   occur <- occur[, cols.tokeep]
  # }
    
    if (length(unique(occur$distance_to_iucn)) == 1) {
      
      iucn.shape <- try(sf::st_read(here::here(paste0('data/raw-data/IUCN/vectors_ETRS89/', code, '/', code, ".shp"))))
      if (sum(class(iucn.shape) %in% 'try-error') == 0) {
      sp.pt <- sf::st_as_sf(occur, coords = c('X','Y'), crs = 4326)
      sp.pt <- sf::st_transform(sp.pt, sf::st_crs(iucn.shape))
      sp.pt$distance_to_iucn <- as.numeric(sf::st_distance(sp.pt, iucn.shape))
      occur$distance_to_iucn <- sp.pt$distance_to_iucn
      occur = occur[occur$distance_to_iucn == 0,]
      }
      occur <- occur[, cols.tokeep]
    } else {
      occur = occur[occur$distance_to_iucn == 0,]
      occur <- occur[, cols.tokeep]
    }

}
  return(occur)
}
