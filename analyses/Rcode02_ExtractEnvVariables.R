
#########################################################
######### Extraction all over the study area ############
#########################################################

##########################
#### Read env layers #####
##########################

# Climate var. 
lst.files <- list.files(here::here('data/raw-data/SIG/Chelsa/'))
r.chelsa <- terra::rast(here::here(paste0('data/raw-data/SIG/Chelsa/', lst.files)))
names(r.chelsa) <- unlist(lapply(strsplit(lst.files, '_'), function(x) {return(x[4])})) 

# Linear structures 
lst.files <- list.files(here::here('data/raw-data/SIG/LinearStructures/'))
r.lin.struc <- terra::rast(here::here(paste0('data/raw-data/SIG/LinearStructures/', lst.files)))
r.lin.struc <- terra::subset(r.lin.struc, subset = grep('prop_tot', names(r.lin.struc)))
names(r.lin.struc) <- unlist(lapply(strsplit(lst.files, '_'), function(x) {return(x[1])})) 

# Topography
lst.files <- list.files(here::here('data/raw-data/SIG/Topography/'))
r.topo <- terra::rast(here::here(paste0('data/raw-data/SIG/Topography/', lst.files)))
names(r.topo)[names(r.topo) == 'prop_tot'] <- 'LakeProportion' 
r.topo <- terra::subset(r.topo, subset = c(1:3))

# Land composition
r.land.syst <- terra::rast(here::here('data/raw-data/SIG/LandComposition/LandCover_Yue_WithConiferousForest_WithNA.tif'))
names(r.land.syst) <- 'LandSystem'
r.waw <- terra::rast(here::here('data/raw-data/SIG/LandComposition/WAWProportion_2012-2018_LandCopernicus.tif'))
names(r.waw) <- 'WAW_prop'

# Project all layers into the same CRS 
r.chelsa <- terra::project(r.chelsa, 'epsg:3035')
r.lin.struc <- terra::project(r.lin.struc, 'epsg:3035')
r.topo <- terra::project(r.topo, 'epsg:3035')
r.land.syst <- terra::project(r.land.syst, 'epsg:3035', method = 'mode')
r.waw <- terra::project(r.waw, 'epsg:3035')

##########################
###### Extraction ########
##########################
# Read grid 
grid.fr <- sf::st_read(here::here('data/raw-data/SIG/Grids/ReferenceGrid_France_bin_1000m.gpkg'))
grid.fr.c <- sf::st_centroid(grid.fr)
grid.fr.c <- sf::st_transform(grid.fr.c, sf::st_crs(r.chelsa))

# Extract climatic var.
var <- terra::extract(x = r.chelsa, y = grid.fr.c)
colnames(var)[colnames(var)!='ID'] <- paste0('climatic.', colnames(var)[colnames(var)!='ID'])

# Extract linear structure var.
var2 <- terra::extract(x = r.lin.struc, y = grid.fr.c)
colnames(var2)[colnames(var2)!='ID'] <- paste0('lin.struct.', colnames(var2)[colnames(var2)!='ID'])

# Extract topo var.
var3 <- terra::extract(x = r.topo, y = grid.fr.c)
colnames(var3)[colnames(var3)!='ID'] <- paste0('topo.', colnames(var3)[colnames(var3)!='ID'])

# Extract land system var.
var4 <- terra::extract(x = r.land.syst, y = grid.fr.c)
interm <- as.list(var4$LandSystem)
lb <- unique(var4$LandSystem)
interm <- lapply(interm, function(x){return( lb %in% x)})
toupdate <- matrix(unlist(interm), ncol = length(lb), byrow = T)
colnames(toupdate) = lb
var4 <- cbind(var4, toupdate)
var4 <- var4[, !is.na(colnames(var4))]
colnames(var4)[colnames(var4)!='ID'] <- paste0('land.compo.', colnames(var4)[colnames(var4)!='ID'])
var5 <- terra::extract(x = r.waw, y = grid.fr.c, layer = 1)
colnames(var5)[colnames(var5)!='ID'] <- paste0('land.compo.', colnames(var5)[colnames(var5)!='ID'])
var5$land.compo.WAW_prop <- as.numeric(paste0(var5$land.compo.WAW_prop))

# Put all together 
var <- dplyr::left_join(var, var2, by = 'ID')
var <- dplyr::left_join(var, var3, by = 'ID')
var <- dplyr::left_join(var, var4, by = 'ID')
var <- dplyr::left_join(var, var5, by = 'ID')
var <- var[, colnames(var)!='ID']
rownames(var) <- paste0('pixel_', c(1:nrow(var)))

saveRDS(var, here::here('data/derived-data/EnvironmentalVariables_France_Res1000m'))

#################################################
########## Prep. occurrence file ################
#################################################
# Read list of SNAP species 
lst.sp <- openxlsx::read.xlsx(here::here('data/derived-data/FinalList-of-SNAP-Vertebrate-Species_ActT-Diet-ForagS-NestH-Morpho-HabPref-DispDTraits.xlsx'))

final <- data.frame()
for (s in unique(lst.sp$LB_NOM_VALIDE_SPE_LEVEL_SYN)) {
  
  s <- gsub(' ','_', s)
  
  occur <- sf::st_read(here::here(paste0('data/derived-data/SIG/Occurrence/AllYearsCombined/GBIFOccurrenceData_France_',s,'_Res1000m_2010-2020.gpkg')))
  occur$Nocc[occur$Nocc > 0] <- 1
  occur <- sf::st_drop_geometry(occur)
  
  if (class(try(cbind(final, data.frame(Nocc = occur$Nocc)))) == 'try-error') {
    final <- data.frame(Nocc = occur$Nocc)
  } else {
    final <- cbind(final, data.frame(Nocc = occur$Nocc))
  }
  colnames(final)[ncol(final)] <- s
}

rownames(final) <- paste0('pixel_', c(1:nrow(final)))
saveRDS(final, here::here('data/derived-data/SNAP-Vertebrate-Species-GBIFOccurrenceData_France_Res1000m_2010-2020'))
