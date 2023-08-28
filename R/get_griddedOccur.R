#' Function to projet GBIF occurence on a grid 
#'
#' @param sp.occur dataframe of geolocated occurrences 
#' @param year.todo year of occurrences to keep 
#' @param grid grid in polygon format
#' @param grid.rast grid in raster format 
#' @param nlim minimum number of occurrences per cell to proceed the analyses 
#' @param crs.sp.occur CRS of geolocated occurrences in sp.occur
#' @param colcoords.sp.occur vector of names of the XY (coordinates) columns in sp.occur
#' @param yearly whether to produce yearly occurrence grids, default = TRUE
#' @param save.raster whether to save the grids as raster (in addition to the geopackages), default = FALSE
#'
#' @return Geopackage of the grid per species with the number of occurrences per cell. 
#' Text file is also generated to provide info on whether at least one cell contains at least "nlim" occurrences. 
#' SIG layers are saved in data/derived-data/SIG/Occurrence/AllYearsCombined or data/derived-data/SIG/Occurrence/Yearly 
#' Textfile is saved in data/derived-data/SIG/Occurrence/
#' @export
#'
#' @examples
get_griddedOccur <- function(sp.occur, year.todo = c(2010:2022), 
                             grid = rr.poly, grid.rast = rr, nlim = 25, crs.sp.occur = 4326, 
                             colcoords.sp.occur = c('X', 'Y'), yearly = F, save.raster = T) {
  
  sum.up = data.frame()
  
  if (nrow(sp.occur) > 0) {

    if (yearly) {

    lst.year <- year.todo[year.todo %in% unique(sp.occur$year)]

    for (y in lst.year) {

      sp.occur.y <- sp.occur[sp.occur$year == y,]

      if(nrow(sp.occur.y) >= nlim) {

        sp.occur.y <- sf::st_as_sf(sp.occur.y, coords = colcoords.sp.occur, crs = crs.sp.occur)
        sp.occur.y <- sf::st_transform(sp.occur.y, sf::st_crs(grid))
        inter <- sf::st_intersects(grid, sp.occur.y)
        inter <- lapply(inter, length)
        grid$Nocc <- unlist(inter)
        
        if (sum(grid$Nocc!=0) >= nlim) {
        sf::st_write(grid, here::here(paste0('data/derived-data/SIG/Occurrence/Yearly/GBIFOccurrenceData_France_',
                                             gsub(' ','_',unique(sp.occur$species)),'_Res1000m_',y , '.gpkg')), driver = 'GPKG', delete_layer = T)

        if (save.raster) {
          grid.rast <- terra::rasterize(grid, grid.rast, field = 'Nocc')
          grid.rast$species <- unique(sp.occur.y$species)
          grid.rast$year <- y
          terra::writeRaster(grid.rast, here::here(paste0('data/derived-data/SIG/Occurrence/Yearly/GBIFOccurrenceData_France_',
                                                          gsub(' ','_',unique(sp.occur$species)),'_Res1000m_', y,'.tif')))
        }

        sum.up <- rbind(sum.up, data.frame(species = unique(sp.occur$species), year = y , AllYearsCombined = NA, Doable = 1))
        }
      }
    }
  }

    # All years combined
    if (nrow(sp.occur) >= nlim) {

      sp.occur <- sp.occur[sp.occur$year %in% year.todo,]
      sp.occur <- sf::st_as_sf(sp.occur, coords = colcoords.sp.occur, crs = crs.sp.occur)
      sp.occur <- sf::st_transform(sp.occur, sf::st_crs(grid))
      inter <- sf::st_intersects(grid, sp.occur)
      inter <- lapply(inter, length)
      grid$Nocc <- unlist(inter)
      
      if (sum(grid$Nocc!=0) >= nlim) {
      sf::st_write(grid, here::here(paste0('data/derived-data/SIG/Occurrence/AllYearsCombined/GBIFOccurrenceData_France_',
                                           gsub(' ','_',unique(sp.occur$species)),'_Res1000m_',
                                           paste0(min(year.todo), '-', max(year.todo)), '.gpkg')), driver = 'GPKG', delete_layer = T)

      if (save.raster) {
        grid.rast <- terra::rasterize(grid, grid.rast, field = 'Nocc')
        grid.rast$species <- unique(sp.occur$species)
        grid.rast$year <- paste0(min(year.todo), '-', max(year.todo))
        terra::writeRaster(grid.rast, here::here(paste0('data/derived-data/SIG/Occurrence/AllYearsCombined/GBIFOccurrenceData_France_',
                                                        gsub(' ','_',unique(sp.occur$species)),'_Res1000m_',
                                                        paste0(min(year.todo), '-', max(year.todo)) ,'.tif')), overwrite = T)
      }
      sum.up <- rbind(sum.up, data.frame(species = unique(sp.occur$species), year = NA , AllYearsCombined = paste0(min(year.todo), '-', max(year.todo)), Doable = 1))
      }
    }
  }
  
  openxlsx::write.xlsx(sum.up, here::here(paste0('data/derived-data/SIG/Occurrence/Evaluation/Evaluation_', gsub(' ','_',unique(sp.occur$species)), '.xlsx')))

  }
