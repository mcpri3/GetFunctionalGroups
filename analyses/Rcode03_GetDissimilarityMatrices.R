library(ggplot2)
library(ggrepel)
library(utils)
library(foreach)
library(ecospat)
library(ade4)

#################################
# Run PCAe to get niche overlap
#################################

# Read data
env <- readRDS(here::here('data/derived-data/EnvironmentalVariables_France_Res1000m'))
env$climatic.swe[is.na(env$climatic.swe)] <- 0
env <- env[, colnames(env) != 'climatic.prsd']
env <- env[, colnames(env) != 'climatic.cdd']
env <- env[, colnames(env) != 'climatic.fd']
env <- env[, colnames(env) != 'land.compo.LandSystem']
env <- env[, colnames(env) != 'topo.SDalti']
env <- env[, colnames(env) != 'climatic.gdd5']
env <- env[, colnames(env) != 'climatic.swe']
env <- env[, colnames(env) != 'climatic.scd']

occur <- readRDS(here::here('data/derived-data/SNAP-Vertebrate-Species-GBIF-INPN-IUCNOccurrenceData_France_Res1000m_2010-2020'))
# Get taxonomic groups 
vert.snap <- openxlsx::read.xlsx(here::here('data/derived-data/FinalList-of-SNAP-Vertebrate-Species_ActT-Diet-ForagS-NestH-Morpho-HabPref-DispD-MoveMod-LifeHist-PressureTraits.xlsx'))
vert.snap <- vert.snap[, c('GROUP2_INPN', 'LB_NOM_VALIDE_SPE_LEVEL_SYN', 'mov.mode.swimmer')]
vert.snap <- dplyr::distinct(vert.snap)
vert.snap$METAGROUP <- ifelse(vert.snap$GROUP2_INPN %in% c('Amphibiens', 'Reptiles'), 'Anura-Squamata-Testudines-Urodela', vert.snap$GROUP2_INPN)
vert.snap$METAGROUP <- gsub('Mammifères', 'Mammalia', vert.snap$METAGROUP)
vert.snap$METAGROUP <- gsub('Oiseaux', 'Aves', vert.snap$METAGROUP)

############## Get ePCA ###################
############## PCA 1 ###################
# Selection of variables 
idx <- c(grep('land.compo', colnames(env)), grep('lin.struct', colnames(env)))
subenv <- env[, idx]

## Calculate PCA for all environment
pca.env <- ade4::dudi.hillsmith(subenv, scannf = F, nf = 2)
var.coord <- pca.env$c1
var.lb <- lapply(strsplit(rownames(var.coord), '.', fixed = T), function(x){return(x[length(x)])})

# Plot
p1 <- ggplot() + 
  geom_point(aes(x = var.coord$CS1, y = var.coord$CS2)) +
  xlim(-1 ,1) +
  ylim(-1, 1) + 
  geom_segment(aes(x = 0, y = 0, xend = var.coord$CS1, yend = var.coord$CS2),
                 arrow = arrow(length = unit(0.3, "cm"))) +
  geom_label_repel(aes(x=var.coord$CS1, y = var.coord$CS2, label = var.lb), max.overlaps = 80, col = '#696969') +
  xlab('Axis 1') + 
  ylab('Axis 2') +
  # theme_minimal() +
  ggtitle('Environmental PCA on landscape composition and linear structures')
p1
ggsave(p1, path = here::here('figures/'), filename = 'PCAe_LandCompo_LinStruct.png', width = 17, height = 7.4)


############## PCA 2 ###################
# Selection of variables 
idx <- c(grep('climatic.', colnames(env)), grep('topo.', colnames(env)))
subenv <- env[, idx]

## Calculate PCA for all environment
pca.env <- dudi.hillsmith(subenv, scannf = F, nf = 2)
var.coord <- pca.env$c1
var.lb <- lapply(strsplit(rownames(var.coord), '.', fixed = T), function(x){return(x[length(x)])})

p2 <- ggplot() + 
  geom_point(aes(x = var.coord$CS1, y = var.coord$CS2)) +
  xlim(-1 ,1) +
  ylim(-1, 1) + 
  geom_segment(aes(x = 0, y = 0, xend = var.coord$CS1, yend = var.coord$CS2),
               arrow = arrow(length = unit(0.3, "cm"))) +
  geom_label_repel(aes(x=var.coord$CS1, y = var.coord$CS2, label = var.lb), max.overlaps = 80, col = '#696969') +
  xlab('Axis 1') + 
  ylab('Axis 2') +
  # theme_minimal() +
  ggtitle('Environmental PCA on climatic conditions and topography')
p2
ggsave(p2, path = here::here('figures/'), filename = 'PCAe_Climatic_Topo.png', width = 17, height = 7.4)

# # Run ePCA using FATE 
lst.acp <- list(c('land.compo', 'lin.struct'), c('topo', 'climatic'))
nms.group.acp <- c('Land cover', 'Abiotic conditions')
names(lst.acp) <- nms.group.acp

for (g in unique(vert.snap$METAGROUP)) {
  
  sp.tokeep <- vert.snap$LB_NOM_VALIDE_SPE_LEVEL_SYN[vert.snap$METAGROUP == g]
  sp.tokeep <- gsub(' ', '_', sp.tokeep)
  suboccur <- occur[, sp.tokeep]
  
for (l in nms.group.acp) {
  
  col.id <- lst.acp[names(lst.acp) == l][[1]]
  col.tokeep <- c()
  for (c in col.id) {
    col.tokeep <- c(col.tokeep, grep(c, colnames(env)))
  }
  
  subenv <- env[, col.tokeep]
  mat.dist.l <- PRE_FATE.speciesDistanceOverlap(mat.overlap.object = list(tab.dom.PA = suboccur, tab.env = subenv), nme.pca = gsub(' ','-', l))
  saveRDS(mat.dist.l, here::here(paste0('data/derived-data/DissimilarityMatrices/SNAP-Vertebrate-Species_PairwiseNicheOverlap_', g,'_', gsub(' ', '-', l ))))
}
}

#####################################################
# Run Gower D to get pairwise trait dissimilarity
#####################################################
# Read data
traits <- openxlsx::read.xlsx(here::here('data/derived-data/FinalList-of-SNAP-Vertebrate-Species_ActT-Diet-ForagS-NestH-Morpho-HabPref-DispD-MoveMod-LifeHist-PressureTraits.xlsx'))
traits$METAGROUP <- ifelse(traits$GROUP2_INPN %in% c('Amphibiens', 'Reptiles'), 'Anura-Squamata-Testudines-Urodela', traits$GROUP2_INPN)
traits$METAGROUP <- gsub('Mammifères', 'Mammalia', traits$METAGROUP)
traits$METAGROUP <- gsub('Oiseaux', 'Aves', traits$METAGROUP)

lst.traits <- c('morpho.BM.or.BL', 'diet.breadth','act.time', 'nest.hab.breadth','hab.pref.breadth','forag.strat', 'dispersal_km', 
                'mov.mode.crawler', 'mov.mode.flier', 'mov.mode.walker', 'mov.mode.swimmer', 'life.hist.offspring_per_year_n', 'pressure')

# Scale BM and BL and combine them in one column 
traits$morpho.BL.Herptiles <- ifelse(traits$GROUP2_INPN %in% c('Amphibiens', 'Reptiles'), traits$morpho.Body_length_max_.cm.,
                                 NA)
traits$morpho.BL.Herptiles <- as.numeric(traits$morpho.BL.Herptiles)
traits$morpho.BL.Herptiles <- scale(traits$morpho.BL.Herptiles)

traits$morpho.BM.BirdsMammals <- ifelse(traits$GROUP2_INPN %in% c('Amphibiens', 'Reptiles'), NA,
                                 traits$morpho.BodyMass)
traits$morpho.BM.BirdsMammals <- scale(traits$morpho.BM.BirdsMammals)
traits$morpho.BM.or.BL <- ifelse(traits$GROUP2_INPN %in% c('Amphibiens', 'Reptiles'), traits$morpho.BL.Herptiles, 
                                 traits$morpho.BM.BirdsMammals) 
# Select columns
idx.tokeep <- c()
for (l in lst.traits) {
  idx.tokeep <- c(idx.tokeep, grep(l, colnames(traits), fixed = T))
}

df <- traits[, idx.tokeep]
df$species <- gsub(' ','_', traits$LB_NOM_VALIDE_SPE_LEVEL_SYN)
df$METAGROUP <- traits$METAGROUP
df <- dplyr::distinct(df)

# Run dissimilarity matrix estimation 
lst.traits2 <- list('life.hist', c('forag.strat','diet'), c('dispersal','mov.mode','act.time'), c('nest.hab', 'hab.pref'), 'morpho',  'pressure', 'aquatic')
nms.group.trait <- c('Reproductive traits', 'Foraging behaviour', 'Movement behaviour', 'Habitat requirement', 'Morphology', 'Vulnerability to pressures', 'Aquatic')
names(lst.traits2) <- nms.group.trait

for (g in unique(df$METAGROUP)) {
  df2 <- df
  if (g %in% c('Mammalia', 'Aves')) {colnames(df2)[colnames(df2) == 'mov.mode.swimmer'] <- 'aquatic'}
  
for (l in nms.group.trait) {
  
  col.id <- lst.traits2[names(lst.traits2) == l][[1]]
  col.tokeep <- c()
  for (c in col.id) {
    col.tokeep <- c(col.tokeep, grep(c, colnames(df2)))
  }
  
  subdf <- as.data.frame(df2[df2$METAGROUP == g, col.tokeep])
  mat.trait.dis <- try(as.matrix(FD::gowdis(subdf)))
  if (sum(class(mat.trait.dis) %in% 'try-error') == 0) {
  rownames(mat.trait.dis) <- gsub(' ','_', df$species[df$METAGROUP == g])
  colnames(mat.trait.dis) <- gsub(' ','_', df$species[df$METAGROUP == g])
  saveRDS(mat.trait.dis, here::here(paste0('data/derived-data/DissimilarityMatrices/SNAP-Vertebrate-Species_PairwiseTraitsDissimilarity_', g, '_', gsub(' ', '-' ,l ))))
  }
  }
}

#####################################################
# Combine matrices 
#####################################################
rm(list = ls())

# Mode 1 : traits and environment weight proportional to their number = 0.125 each (=1/8 = 1/num. of mat)
lst.mat.files <- list.files(here::here('data/derived-data/DissimilarityMatrices/'))
lst.mat.files <- lst.mat.files[grep('Pairwise', lst.mat.files)]

if (length(grep('Combined', lst.mat.files)) != 0) {lst.mat.files <- lst.mat.files[-grep('Combined', lst.mat.files)]}
if (length(grep('AllTraits', lst.mat.files)) != 0) {lst.mat.files <- lst.mat.files[-grep('AllTraits', lst.mat.files)]} 

for (g in c('Aves', 'Mammalia', 'Anura-Squamata-Testudines-Urodela')) {
  
  lst.mat.files.gp <- lst.mat.files[grep(g, lst.mat.files)]
  remove(list = paste0('group.', g))
  
  # Weight calculation
  ntrait <- length(grep('Traits', lst.mat.files.gp))
  nenv <- length(grep('Niche', lst.mat.files.gp))
  wtrait <- ntrait/(ntrait + nenv)
  wenv <- nenv/(ntrait + nenv)
  wtrait <- wtrait/ntrait 
  wenv <- wenv/nenv 
  
  w.vec <- c()
  
  for (l in lst.mat.files.gp) {
  mmat <- readRDS(paste0(here::here('data/derived-data/DissimilarityMatrices/', l)))
  mmat <- as.matrix(mmat)
  assign(paste0('group.', g, '.', which(l == lst.mat.files.gp)), mmat)
  w.val = ifelse(length(grep('Trait', l)) == 1, wtrait, wenv)
  w.vec <- c(w.vec, w.val)
  }
  
  lys <- ls()[grep(g, ls())]
  assign(paste0('group.', g), mget(x = lys))
  remove(list = lys)  
  
  
  mat.dist <- PRE_FATE.speciesDistanceCombine(list.mat.dist = get(paste0('group.', g)), opt.weights = w.vec)
  saveRDS(mat.dist, here::here(paste0('data/derived-data/DissimilarityMatrices/SNAP-Vertebrate-Species_PairwiseCombinedDissimilarity_TraitsEnv_', g, '_Weight-EqualTraitsEnv')))
  
}

# Mode 2 : traits and environment weight half and half 
lst.mat.files <- list.files(here::here('data/derived-data/DissimilarityMatrices/'))
lst.mat.files <- lst.mat.files[grep('Pairwise', lst.mat.files)]

if (length(grep('Combined', lst.mat.files)) != 0) {lst.mat.files <- lst.mat.files[-grep('Combined', lst.mat.files)]}
if (length(grep('AllTraits', lst.mat.files)) != 0) {lst.mat.files <- lst.mat.files[-grep('AllTraits', lst.mat.files)]} 

for (g in c('Aves', 'Mammalia', 'Anura-Squamata-Testudines-Urodela')) {
  
  lst.mat.files.gp <- lst.mat.files[grep(g, lst.mat.files)]
  remove(list = paste0('group.', g))
  
  # Weight calculation
  ntrait <- length(grep('Traits', lst.mat.files.gp))
  nenv <- length(grep('Niche', lst.mat.files.gp))
  wtrait <- 1/2
  wenv <- 1/2
  wtrait <- wtrait/ntrait 
  wenv <- wenv/nenv 
  
  w.vec <- c()
  
  for (l in lst.mat.files.gp) {
    mmat <- readRDS(paste0(here::here('data/derived-data/DissimilarityMatrices/', l)))
    mmat <- as.matrix(mmat)
    assign(paste0('group.', g, '.', which(l == lst.mat.files.gp)), mmat)
    w.val = ifelse(length(grep('Trait', l)) == 1, wtrait, wenv)
    w.vec <- c(w.vec, w.val)
  }
  
  lys <- ls()[grep(g, ls())]
  assign(paste0('group.', g), mget(x = lys))
  remove(list = lys)  
  
  
  mat.dist <- PRE_FATE.speciesDistanceCombine(list.mat.dist = get(paste0('group.', g)), opt.weights = w.vec)
  saveRDS(mat.dist, here::here(paste0('data/derived-data/DissimilarityMatrices/SNAP-Vertebrate-Species_PairwiseCombinedDissimilarity_TraitsEnv_', g, '_Weight-HalfTraits-HalfEnv')))
  
}

