#################################
# Run PCAe to get niche overlap
#################################

# Read data
env <- readRDS(here::here('data/derived-data/EnvironmentalVariables_France_Res1000m'))
env <- env[, colnames(env) != 'climatic.swe']
env <- env[, colnames(env) != 'land.compo.LandSystem']
occur <- readRDS(here::here('data/derived-data/SNAP-Vertebrate-Species-GBIFOccurrenceData_France_Res1000m_2010-2020'))

# Run PCA 
lst.acp <- c('climatic','lin.struct', 'topo', 'land.compo')

for (l in lst.acp) {
  sub.env <- env[, grep(l, colnames(env), fixed = T)]
  sub.env <- na.omit(sub.env)
  mat.dist.l <- PRE_FATE.speciesDistanceOverlap(mat.overlap.object = list(tab.dom.PA = occur, tab.env = sub.env))
  assign(paste0('niche.overlap.dissim.', l), mat.dist.l)
  saveRDS(mat.dist.l, here::here(paste0('data/derived-data/SNAP-Vertebrate-Species_PairwiseNicheOverlap_', l )))
}

#####################################################
# Run Gower D to get pairwise trait dissimilarity
#####################################################
# Read data
traits <- openxlsx::read.xlsx(here::here('data/derived-data/FinalList-of-SNAP-Vertebrate-Species_ActT-Diet-ForagS-NestH-Morpho-HabPref-DispDTraits.xlsx'))
lst.traits <- c('morpho.BM.or.BL', 'diet','act.time', 'nest.hab','hab.pref','forag.strat', 'dispersal_km')

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
df$species <- gsub(' ','_', traits$LB_NOM_VALIDE_SPE_LEVEL)

# Run dissimilarity matrix estimation 
mat.dist <- PRE_FATE.speciesDistanceTraits(mat.traits = df)
saveRDS(mat.dist, here::here('data/derived-data/SNAP-Vertebrate-Species_PairwiseTraitsDissimilarity'))

#####################################################
# Combine matrices 
#####################################################

lst.mat.files <- list.files(here::here('data/derived-data/'))
lst.mat.files <- lst.mat.files[grep('Pairwise', lst.mat.files)]

for (l in lst.mat.files) {
  mmat <- readRDS(paste0(here::here('data/derived-data/', l)))
  mmat <- as.matrix(mmat)
  assign(gsub('SNAP-Vertebrate-Species_', '', l), mmat)
}

lst.mat <- list(PairwiseNicheOverlap_climatic, PairwiseNicheOverlap_land.compo, PairwiseNicheOverlap_lin.struct, PairwiseNicheOverlap_topo , PairwiseTraitsDissimilarity)
mat.dist <- PRE_FATE.speciesDistanceCombine(list.mat.dist = lst.mat)
saveRDS(mat.dist, here::here('data/derived-data/SNAP-Vertebrate-Species_PairwiseCombinedDissimilarity_TraitsEnv'))

#####################################################
# Get groups 
#####################################################
mat.dist <- readRDS(here::here('data/derived-data/SNAP-Vertebrate-Species_PairwiseCombinedDissimilarity_TraitsEnv'))
dendro <- PRE_FATE.speciesClustering_step1(mat.species.DIST = list(mat.dist), opt.no_clust_max = 15)
groups <- PRE_FATE.speciesClustering_step2(clust.dendrograms = dendro$clust.dendrograms, mat.species.DIST = mat.dist, no.clusters = 13)
