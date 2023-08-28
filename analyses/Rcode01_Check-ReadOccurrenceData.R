library(dplyr)
library(ggplot2)

##################################################
# Prep. GBIF occurrence files 
##################################################
# List GBIF files 
lst.files <- list.files(here::here('data/raw-data/GBIF/Chordata/'), recursive = T) 
lst.files <- lst.files[grep('occurencesGBIF', lst.files)]
lst.files.df <- data.frame(file.nme = lst.files)
lst.files <- unlist(lapply(as.list(lst.files), get_speciesname))
lst.files.df$SpeciesName <- lst.files

# Read list of SNAP species 
vert.snap <- readxl::read_excel(here::here("data/raw-data/SNAP-Vertebrate-Species_ActT-Diet-ForagS-NestH-Morpho-HabPref-DispD-MoveMod-LifeHist-PressureTraits.xlsx"))
# Manually checked which one are missing and manually correct one 
vert.snap$LB_NOM_VALIDE_SPE_LEVEL_SYN[vert.snap$LB_NOM_VALIDE_SPE_LEVEL_SYN == 'Phylloscopus sibilatrix'] <- 'Phylloscopus sibillatrix'

# Read, clean and combine occurrence files 
lst.files.df <- lst.files.df[lst.files.df$SpeciesName %in% vert.snap$LB_NOM_VALIDE_SPE_LEVEL_SYN,]
tobind <- vert.snap[, c('LB_NOM_VALIDE_SPE_LEVEL_SYN', "Code")]
lst.files.df <- left_join(lst.files.df, tobind, by = c("SpeciesName" = "LB_NOM_VALIDE_SPE_LEVEL_SYN"))
lst.files.df = distinct(lst.files.df)
# occur <- group_by(lst.files.df, file.nme) %>% do(clean_occurGBIF(ffile = .$file.nme, year.min = 2000, uncertainty = 1000)) %>% data.frame # filter : resolution < 1000m, year after 2000, Species rank 

###################################################################################
# Create occurrence raster files from GBIF occurrence data (yearly & years mixed)
###################################################################################
# Prep. grid 
rr <- terra::rast(here::here('data/raw-data/SIG/Grids/ReferenceGrid_Europe_bin_1000m.tif'))
france <- sf::st_read(here::here('data/raw-data/SIG/Grids/CasestudyOutlines.gpkg'), layer = 'France')
france <- sf::st_transform(france, sf::st_crs(rr))
rr <- terra::crop(rr, france)
# rr.poly <- terra::as.polygons(rr, dissolve = F)
# rr.poly <- rr.poly[rr.poly$layer == 1, ]
# rr.poly <- sf::st_as_sf(rr.poly)
# ddist <- sf::st_distance(x = rr.poly, y = france)
# rr.poly <- rr.poly[as.numeric(ddist[,1]) == 0,]
# sf::st_write(rr.poly, here::here('data/raw-data/SIG/Grids/ReferenceGrid_France_bin_1000m.gpkg'), driver = 'GPKG', delete_layer = T)
rr.poly <- sf::st_read(here::here('data/raw-data/SIG/Grids/ReferenceGrid_France_bin_1000m.gpkg'))
# group_by(occur, species) %>% do(get_griddedOccur(sp.occur = .))

group_by(lst.files.df, file.nme) %>% do(clean_occurGBIF(ffile = .$file.nme, code = .$Code, year.min = 2010, uncertainty = 1000) %>% get_griddedOccur(sp.occur = ., year.todo = c(2010:2020), nlim = 5))

############################################################
# Create occurrence raster files from INPN occurrence data 
############################################################
# Load data summary
load(here::here("data/raw-data/INPN/ListEsp_2022-10-10.RData")) # Species list
# lst.sp <- final$SpeciesName[final$AllyrsComb == 1]
# missing <- vert.snap[!(vert.snap$LB_NOM_VALIDE_SPE_LEVEL_SYN %in% lst.sp), ]
missing = vert.snap
nlim = 5

ListEsp <- ListEsp[ListEsp$cd_ref %in% missing$CD_REF_SPE_LEVEL, ]
inpn <- data.frame()
for (s in unique(ListEsp$cd_ref)) {
  subdata <- ListEsp[ListEsp$cd_ref == s,]
  inpn <- rbind(inpn, data.frame(cd_ref = s, lb_nom_valide = unique(subdata$lb_nom_valide), NbOcc = sum(subdata$NbOcc)))
}
inpn <- inpn[inpn$NbOcc >= nlim,]

# Load occurrence
load("~/Documents/LECA/NaturaConnect/Rprojects/02_GetFunctionalGroups/data/raw-data/INPN/Data_2022-10-10.RData")
Data <- Data[Data$cd_ref %in% inpn$cd_ref, ]
# Filter time period
Data$date_debut <- as.Date(Data$date_debut, format = "%Y-%m-%d")
Data$date_fin <- as.Date(Data$date_fin, format = "%Y-%m-%d")
Data <- Data[Data$date_debut >= as.Date("2010-01-01", format = "%Y-%m-%d"), ]
Data <- Data[Data$date_fin <= as.Date("2020-12-31", format = "%Y-%m-%d"), ]
# Check num. occurrences
Nocc.inpn <- data.frame(cd_ref = names(table(Data$cd_ref)), Nocc = c(table(Data$cd_ref)))
Nocc.inpn <- Nocc.inpn[Nocc.inpn$Nocc >= nlim, ]
Data <- Data[Data$cd_ref %in% Nocc.inpn$cd_ref, ]
# Load geometries
load("~/Documents/LECA/NaturaConnect/Rprojects/02_GetFunctionalGroups/data/raw-data/INPN/GeomCom_2022-10-10.RData")
GeomCom <- GeomCom[GeomCom$cd_sig_geometrie %in% Data$cd_sig_geometrie, ]
GeomCom.c <- sf::st_centroid(GeomCom)
GeomCom.c <- GeomCom.c[, c("cd_sig_geometrie", 'epsg_local')]
GeomCom.c$X <- sf::st_coordinates(GeomCom.c)[, 1]
GeomCom.c$Y <- sf::st_coordinates(GeomCom.c)[, 2]
GeomCom.c <- sf::st_drop_geometry(GeomCom.c)
Data <- left_join(Data,GeomCom.c, by = 'cd_sig_geometrie')

todo <- unique(Data$cd_ref)
sum.up <- data.frame()
year.todo = c(2010:2020)
grid = rr.poly
grid.rast = rr

for (s in todo) {
  
  lb.sp <- missing$LB_NOM_VALIDE_SPE_LEVEL_SYN[missing$CD_REF_SPE_LEVEL == s]
  sp.occur <- Data[Data$cd_ref == s,]
  sp.occur <- sf::st_as_sf(sp.occur, coords = c('X','Y'), crs = unique(sp.occur$epsg_local))
  
  # # Check distance to IUCN
  # code <- missing$Code[missing$CD_REF_SPE_LEVEL == s]
  # iucn.shape <- try(sf::st_read(here::here(paste0('data/raw-data/IUCN/vectors_ETRS89/', code, '/', code, ".shp"))))
  # if (sum(class(iucn.shape) %in% 'try-error') == 0) {
  #   sp.occur <- sf::st_transform(sp.occur, sf::st_crs(iucn.shape))
  #   sp.occur$distance_to_iucn <- as.numeric(sf::st_distance(sp.occur, iucn.shape))
  #   sp.occur = sp.occur[sp.occur$distance_to_iucn == 0,]
  # }
  #
  # if (nrow(sp.occur)>0) {
  
  sp.occur <- sf::st_transform(sp.occur, sf::st_crs(grid))
  inter <- sf::st_intersects(grid, sp.occur)
  inter <- lapply(inter, length)
  grid$Nocc <- unlist(inter)
  
  if (sum(grid$Nocc!=0) >= nlim) {
    sf::st_write(grid, here::here(paste0('data/derived-data/SIG/Occurrence/AllYearsCombined/INPNOccurrenceData_France_',
                                         gsub(' ','_', lb.sp),'_Res1000m_',
                                         paste0(min(year.todo), '-', max(year.todo)), '.gpkg')), driver = 'GPKG', delete_layer = T)
    
    
    grid.rast <- terra::rasterize(grid, grid.rast, field = 'Nocc')
    grid.rast$species <- unique(sp.occur$species)
    grid.rast$year <- paste0(min(year.todo), '-', max(year.todo))
    terra::writeRaster(grid.rast, here::here(paste0('data/derived-data/SIG/Occurrence/AllYearsCombined/INPNOccurrenceData_France_',
                                                    gsub(' ','_', lb.sp),'_Res1000m_',
                                                    paste0(min(year.todo), '-', max(year.todo)) ,'.tif')), overwrite = T)
    
    sum.up <- rbind(sum.up, data.frame(SpeciesName = lb.sp, AllyrsComb = 1 , Nyrs = 0, yrs.min = NA, yrs.max = NA))
  }
  # }
}
openxlsx::write.xlsx(sum.up, here::here('data/derived-data/SIG/Occurrence/Evaluation/Evaluation_INPN.xlsx'))

#########################################################################
# Combine both raster when species are available for both dataset
#########################################################################

#############################################
# Check species doable from GBIF
#############################################
lst.files <- list.files(here::here('data/derived-data/SIG/Occurrence/Evaluation'))
lst.files <- lst.files[lst.files != c('Evaluation_INPN.xlsx')]
lst.files <- lst.files[lst.files != c('Evaluation_IUCN.xlsx')]

lst.files <- lst.files[grep('Evaluation', lst.files)]

final <- data.frame()

for (f in lst.files) {
  
  eval <- readxl::read_excel(here::here(paste0('data/derived-data/SIG/Occurrence/Evaluation/', f)))
  sp <- strsplit(f, 'Evaluation_')[[1]][2]
  sp <- gsub('.xlsx','', sp)
  sp <- gsub('_',' ', sp)
  
  if (nrow(eval) == 0) {
    final <- rbind(final, data.frame(SpeciesName = sp, AllyrsComb = 0, Nyrs = 0, yrs.min = NA, yrs.max = NA))  
  } else {
    
    subeval = eval[!is.na(eval$year),]
    
    if (nrow(subeval)>0) {
      
      yrs.min = min(na.omit(eval$year))
      yrs.max = max(na.omit(eval$year))
      
    } else {
      yrs.min = NA
      yrs.max = NA
    }
    
    final <- rbind(final, data.frame(SpeciesName = sp, AllyrsComb = 1,
                                     Nyrs = nrow(eval)-1, yrs.min = yrs.min, 
                                     yrs.max = yrs.max))  
  }
}

sum(final$AllyrsComb)

#############################################
# Add species doable from INPN 
#############################################
sum.up <- openxlsx::read.xlsx(here::here('data/derived-data/SIG/Occurrence/Evaluation/Evaluation_INPN.xlsx'))
final <- final[final$AllyrsComb == 1, ]
final <- rbind(final, sum.up)

#############################################
# Combine both 
#############################################
double <- table(final$SpeciesName)
double <- double[double == 2]
year.todo <- c(2010:2020)

for (lb.sp in names(double)) {
  
  grid.gbif <- terra::rast(here::here(paste0('data/derived-data/SIG/Occurrence/AllYearsCombined/GBIFOccurrenceData_France_',
                                             gsub(' ','_', lb.sp),'_Res1000m_',
                                             paste0(min(year.todo), '-', max(year.todo)) ,'.tif')))
  file.copy(here::here(paste0('data/derived-data/SIG/Occurrence/AllYearsCombined/GBIFOccurrenceData_France_',
                              gsub(' ','_', lb.sp),'_Res1000m_',
                              paste0(min(year.todo), '-', max(year.todo)) ,'.tif')), here::here(paste0('data/derived-data/SIG/Occurrence/AllYearsCombined/KeptButNotUsedAnymore/GBIFOccurrenceData_France_',
                                                                                                       gsub(' ','_', lb.sp),'_Res1000m_',
                                                                                                       paste0(min(year.todo), '-', max(year.todo)) ,'.tif')))
  
  shape.gbif <- sf::st_read(here::here(paste0('data/derived-data/SIG/Occurrence/AllYearsCombined/GBIFOccurrenceData_France_',
                                              gsub(' ','_', lb.sp),'_Res1000m_',
                                              paste0(min(year.todo), '-', max(year.todo)), '.gpkg')))
  file.copy(here::here(paste0('data/derived-data/SIG/Occurrence/AllYearsCombined/GBIFOccurrenceData_France_',
                              gsub(' ','_', lb.sp),'_Res1000m_',
                              paste0(min(year.todo), '-', max(year.todo)), '.gpkg')), here::here(paste0('data/derived-data/SIG/Occurrence/AllYearsCombined/KeptButNotUsedAnymore/GBIFOccurrenceData_France_',
                                                                                                        gsub(' ','_', lb.sp),'_Res1000m_',
                                                                                                        paste0(min(year.todo), '-', max(year.todo)), '.gpkg')))
  
  grid.inpn <- terra::rast(here::here(paste0('data/derived-data/SIG/Occurrence/AllYearsCombined/INPNOccurrenceData_France_',
                                             gsub(' ','_', lb.sp),'_Res1000m_',
                                             paste0(min(year.todo), '-', max(year.todo)) ,'.tif')))
  file.copy(here::here(paste0('data/derived-data/SIG/Occurrence/AllYearsCombined/INPNOccurrenceData_France_',
                              gsub(' ','_', lb.sp),'_Res1000m_',
                              paste0(min(year.todo), '-', max(year.todo)) ,'.tif')), here::here(paste0('data/derived-data/SIG/Occurrence/AllYearsCombined/KeptButNotUsedAnymore/INPNOccurrenceData_France_',
                                                                                                       gsub(' ','_', lb.sp),'_Res1000m_',
                                                                                                       paste0(min(year.todo), '-', max(year.todo)) ,'.tif')))
  
  
  shape.inpn <- sf::st_read(here::here(paste0('data/derived-data/SIG/Occurrence/AllYearsCombined/INPNOccurrenceData_France_',
                                              gsub(' ','_', lb.sp),'_Res1000m_',
                                              paste0(min(year.todo), '-', max(year.todo)), '.gpkg')))
  file.copy(here::here(paste0('data/derived-data/SIG/Occurrence/AllYearsCombined/INPNOccurrenceData_France_',
                              gsub(' ','_', lb.sp),'_Res1000m_',
                              paste0(min(year.todo), '-', max(year.todo)), '.gpkg')), here::here(paste0('data/derived-data/SIG/Occurrence/AllYearsCombined/KeptButNotUsedAnymore/INPNOccurrenceData_France_',
                                                                                                        gsub(' ','_', lb.sp),'_Res1000m_',
                                                                                                        paste0(min(year.todo), '-', max(year.todo)), '.gpkg')))
  shape.both <- shape.gbif
  shape.both$Nocc <- shape.gbif$Nocc + shape.inpn$Nocc    
  shape.both$Nocc[shape.both$Nocc > 1] <- 1
  grid.both <- grid.gbif$Nocc + grid.inpn$Nocc 
  grid.both[grid.both>1] <- 1
  grid.both$species <- lb.sp
  grid.both$year <- paste0(min(year.todo), '-', max(year.todo))
  
  terra::writeRaster(grid.both, here::here(paste0('data/derived-data/SIG/Occurrence/AllYearsCombined/GBIF-INPNOccurrenceData_France_',
                                                  gsub(' ','_', lb.sp),'_Res1000m_',
                                                  paste0(min(year.todo), '-', max(year.todo)) ,'.tif')), overwrite = T)
  sf::st_write(shape.both, here::here(paste0('data/derived-data/SIG/Occurrence/AllYearsCombined/GBIF-INPNOccurrenceData_France_',
                                             gsub(' ','_', lb.sp),'_Res1000m_',
                                             paste0(min(year.todo), '-', max(year.todo)), '.gpkg')), driver = 'GPKG', delete_layer = T)
  
  unlink(here::here(paste0('data/derived-data/SIG/Occurrence/AllYearsCombined/INPNOccurrenceData_France_',
                           gsub(' ','_', lb.sp),'_Res1000m_',
                           paste0(min(year.todo), '-', max(year.todo)), '.gpkg')))
  unlink(here::here(paste0('data/derived-data/SIG/Occurrence/AllYearsCombined/GBIFOccurrenceData_France_',
                           gsub(' ','_', lb.sp),'_Res1000m_',
                           paste0(min(year.todo), '-', max(year.todo)) ,'.tif')))
  unlink(here::here(paste0('data/derived-data/SIG/Occurrence/AllYearsCombined/GBIFOccurrenceData_France_',
                           gsub(' ','_', lb.sp),'_Res1000m_',
                           paste0(min(year.todo), '-', max(year.todo)), '.gpkg')))
  unlink(here::here(paste0('data/derived-data/SIG/Occurrence/AllYearsCombined/INPNOccurrenceData_France_',
                           gsub(' ','_', lb.sp),'_Res1000m_',
                           paste0(min(year.todo), '-', max(year.todo)) ,'.tif')))
  
}

########################################################
# Add IUCN ranges 
########################################################
final <- dplyr::distinct(final)

#IUCN range (see code below)
sum.up2 <- openxlsx::read.xlsx(here::here('data/derived-data/SIG/Occurrence/Evaluation/Evaluation_IUCN.xlsx'))
final <- rbind(final, sum.up2)

vert.snap.trait <- vert.snap[!is.na(vert.snap$dispersal_km),] #species for which we have DD 

lst.sp <- final$SpeciesName[final$SpeciesName %in% vert.snap.trait$LB_NOM_VALIDE_SPE_LEVEL_SYN]
missing <- vert.snap[!(vert.snap$LB_NOM_VALIDE_SPE_LEVEL_SYN %in% lst.sp), ]
vert.snap.final <- vert.snap[vert.snap$LB_NOM_VALIDE_SPE_LEVEL_SYN %in% lst.sp, ]
openxlsx::write.xlsx(vert.snap.final, here::here('data/derived-data/FinalList-of-SNAP-Vertebrate-Species_ActT-Diet-ForagS-NestH-Morpho-HabPref-DispD-MoveMod-LifeHist-PressureTraits.xlsx'))
openxlsx::write.xlsx(missing, here::here('data/derived-data/DataDeficientSpecies/DataDeficient-SNAP-Vertebrate-Species.xlsx'))

missing.trait <- missing[missing$LB_NOM_VALIDE_SPE_LEVEL_SYN %in% final$SpeciesName,]
openxlsx::write.xlsx(missing.trait, here::here('data/derived-data/DataDeficientSpecies/DataDeficientTraits-SNAP-Vertebrate-Species.xlsx'))
missing.occur <- missing[missing$LB_NOM_VALIDE_SPE_LEVEL_SYN %in% vert.snap.trait$LB_NOM_VALIDE_SPE_LEVEL_SYN,]
openxlsx::write.xlsx(missing.occur, here::here('data/derived-data/DataDeficientSpecies/DataDeficientOccurrences-SNAP-Vertebrate-Species.xlsx'))
missing.both <-  missing[!(missing$LB_NOM_VALIDE_SPE_LEVEL_SYN %in% c(missing.trait$LB_NOM_VALIDE_SPE_LEVEL_SYN, missing.occur$LB_NOM_VALIDE_SPE_LEVEL_SYN)), ]
openxlsx::write.xlsx(missing.both, here::here('data/derived-data/DataDeficientSpecies/DataDeficientTraits&Occurrences-SNAP-Vertebrate-Species.xlsx'))



df.plot <- data.frame(nms = names(sort(table(missing.trait$ORDRE))), val = sort(table(missing.trait$ORDRE)))
df.plot$nms <- factor(df.plot$nms, levels = df.plot$nms)
p1 <- ggplot(df.plot) + 
  geom_col(aes(x = nms, y = val.Freq)) + 
  ylab('Number of species') +
  xlab('Order') + 
  theme(axis.text.x = element_text(angle = 45,  hjust=1, size = 10), axis.text.y = element_text(size = 10)) +
  scale_y_continuous(breaks = seq(1, 35, by = 2)) +
  ggtitle(paste0('Trait DD - Total number of missing species = ', nrow(missing.trait)))
ggsave(p1, filename = 'DataDeficientforTraits-SNAP-Vertebrate-Species.png',path = here::here('figures/'))

df.plot <- data.frame(nms = names(sort(table(missing.occur$ORDRE))), val = sort(table(missing.occur$ORDRE)))
df.plot$nms <- factor(df.plot$nms, levels = df.plot$nms)
p2 <- ggplot(df.plot) + 
  geom_col(aes(x = nms, y = val.Freq)) + 
  ylab('Number of species') +
  xlab('Order') + 
  theme(axis.text.x = element_text(angle = 45,  hjust=1, size = 10), axis.text.y = element_text(size = 10)) +
  scale_y_continuous(breaks = seq(1, 35, by = 2)) +
  ggtitle(paste0('Occurrence DD - Total number of missing species = ', nrow(missing.occur)))
ggsave(p2, filename = 'DataDeficientforOccurrence-SNAP-Vertebrate-Species.png',path = here::here('figures/'))

df.plot <- data.frame(nms = names(sort(table(missing.both$ORDRE))), val = sort(table(missing.both$ORDRE)))
df.plot$nms <- factor(df.plot$nms, levels = df.plot$nms)
p3 <- ggplot(df.plot) + 
  geom_col(aes(x = nms, y = val)) + 
  ylab('Number of species') +
  xlab('Order') + 
  theme(axis.text.x = element_text(angle = 45,  hjust=1, size = 10), axis.text.y = element_text(size = 10)) +
  scale_y_continuous(breaks = seq(1, 35, by = 2)) +
  ggtitle(paste0('Trait & Occurrence DD - Total number of missing species = ', nrow(missing.both)))
ggsave(p3, filename = 'DataDeficientforOccurrence&Trait-SNAP-Vertebrate-Species.png',path = here::here('figures/'))

######################################################
########## Check IUCN range ##########################
######################################################
missing.occur <- openxlsx::read.xlsx(here::here('data/derived-data/DataDeficientSpecies/DataDeficientOccurrences-SNAP-Vertebrate-Species.xlsx'))
lst.files <- list.files(here::here('data/raw-data/IUCN/vectors_ETRS89/'))
lst.files <- lst.files[lst.files %in% missing.occur$Code]

# Read. grid 
rr <- terra::rast(here::here('data/raw-data/SIG/Grids/ReferenceGrid_Europe_bin_1000m.tif'))
france <- sf::st_read(here::here('data/raw-data/SIG/Grids/CasestudyOutlines.gpkg'), layer = 'France')
france <- sf::st_transform(france, sf::st_crs(rr))
rr <- terra::crop(rr, france)
rr.poly <- sf::st_read(here::here('data/raw-data/SIG/Grids/ReferenceGrid_France_bin_1000m.gpkg'))
nlim <- 5

sum.up <- data.frame()

for (l in lst.files) {
  
  lb.sp <- missing.occur$LB_NOM_VALIDE_SPE_LEVEL_SYN[missing.occur$Code == l]
  
  sp.occur <- sf::st_read(here::here(paste0('data/raw-data/IUCN/vectors_ETRS89/', l, '/', l, '.shp')))
  sp.occur <- sf::st_transform(sp.occur, sf::st_crs(rr.poly))
  inter <- sf::st_intersects(rr.poly, sp.occur)
  inter <- lapply(inter, length)
  rr.poly$Nocc <- unlist(inter)
  
  if (sum(rr.poly$Nocc!=0) >= nlim) {
    
    sf::st_write(rr.poly, here::here(paste0('data/derived-data/SIG/Occurrence/AllYearsCombined/IUCNOccurrenceData_France_',
                                            gsub(' ','_', lb.sp),'_Res1000m_NoTime.gpkg')), driver = 'GPKG', delete_layer = T)
    
    
    rr <- terra::rasterize(rr.poly, rr, field = 'Nocc')
    rr$species <- lb.sp
    rr$year <- NA
    terra::writeRaster(rr, here::here(paste0('data/derived-data/SIG/Occurrence/AllYearsCombined/IUCNOccurrenceData_France_',
                                             gsub(' ','_', lb.sp),'_Res1000m_NoTime', '.tif')), overwrite = T)
    
    sum.up <- rbind(sum.up, data.frame(SpeciesName = lb.sp, AllyrsComb = 1 , Nyrs = 0, yrs.min = NA, yrs.max = NA))
  }
  
}
openxlsx::write.xlsx(sum.up, here::here('data/derived-data/SIG/Occurrence/Evaluation/Evaluation_IUCN.xlsx'))
