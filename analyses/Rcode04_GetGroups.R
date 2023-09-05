################################
# Get groups 
################################
rm(list = ls())

################################################################
################### MAMMALIA ###################
################################################################

combi.doable <- data.frame()
# fullcolor = grDevices::colors()[grep('dark', grDevices::colors())]
g = 'Mammalia'

############## Traits, env equal ############## 
m = 'EqualTraitsEnv'
mat.dist <- readRDS(here::here(paste0('data/derived-data/DissimilarityMatrices/SNAP-Vertebrate-Species_PairwiseCombinedDissimilarity_TraitsEnv_',g,'_Weight-', m)))
dendro <- PRE_FATE.speciesClustering_step1(mat.species.DIST = list(mat.dist), opt.no_clust_max = 20)
dt.tree <- ape::as.phylo(dendro$clust.dendrograms[[1]])
plot(dt.tree, type = "unrooted", no.margin = TRUE)

# Get groups
k = 11 # 6 et 11 sont ok 
groups <- cutree(dendro$clust.dendrograms[[1]], k = k)
# Plot
# colors <- sample(fullcolor, length(unique(groups)))
colors <- viridisLite::turbo(length(unique(groups)))
colors <-  colorspace::darken(colors, amount = 0.2)
plot(dt.tree, type = "cladogram", no.margin = TRUE, tip.color = colors[groups], cex = 0.8)
# Groups in a DF 
groups <- data.frame(species = names(groups), cluster.id = c(groups))
groups$species <- gsub('_', ' ', groups$species)
lst.sp <- openxlsx::read.xlsx(here::here('data/derived-data/FinalList-of-SNAP-Vertebrate-Species_ActT-Diet-ForagS-NestH-Morpho-HabPref-DispD-MoveMod-LifeHist-PressureTraits.xlsx'))
lst.sp.g <- lst.sp[lst.sp$GROUP2_INPN == 'Mammifères',]
lst.sp.g <- dplyr::left_join(lst.sp.g, groups, c("LB_NOM_VALIDE_SPE_LEVEL_SYN" = "species"))
openxlsx::write.xlsx(lst.sp.g, here::here(paste0('outputs/FinalList-of-SNAP-Vertebrate-Species_ActT-Diet-ForagS-NestH-Morpho-HabPref-DispD-MoveMod-LifeHist-PressureTraits_',g,'_GroupID_K=',
                                                 k, '_Weight-',m,'.xlsx')))
saveRDS(dendro$clust.dendrograms[[1]], here::here(paste0('outputs/HierarchicalClustering_Dendrogram_', g,'_Weight-', m)))

combi.doable <- rbind(combi.doable, data.frame(group = g, Nclus = k, mode = m))

############## Traits, env half half ############## 
m = 'HalfTraits-HalfEnv'
mat.dist <- readRDS(here::here(paste0('data/derived-data/DissimilarityMatrices/SNAP-Vertebrate-Species_PairwiseCombinedDissimilarity_TraitsEnv_',g,'_Weight-', m)))
dendro <- PRE_FATE.speciesClustering_step1(mat.species.DIST = list(mat.dist), opt.no_clust_max = 20)
# Get tree
dt.tree <- ape::as.phylo(dendro$clust.dendrograms[[1]])
plot(dt.tree, type = "unrooted", no.margin = TRUE)
# Get groups
k = 11 # 4 et 11 ok 
groups <- cutree(dendro$clust.dendrograms[[1]], k = k)
# Plot
colors <- viridisLite::turbo(length(unique(groups)))
colors <-  colorspace::darken(colors, amount = 0.2)
plot(dt.tree, type = "cladogram", no.margin = TRUE, tip.color = colors[groups], cex = 0.8)
# Groups in a DF 
groups <- data.frame(species = names(groups), cluster.id = c(groups))
# groups <- PRE_FATE.speciesClustering_step2(clust.dendrograms = dendro$clust.dendrograms, mat.species.DIST = mat.dist, no.clusters = 13)
groups$species <- gsub('_', ' ', groups$species)
lst.sp <- openxlsx::read.xlsx(here::here('data/derived-data/FinalList-of-SNAP-Vertebrate-Species_ActT-Diet-ForagS-NestH-Morpho-HabPref-DispD-MoveMod-LifeHist-PressureTraits.xlsx'))
lst.sp.g <- lst.sp[lst.sp$GROUP2_INPN == 'Mammifères',]
lst.sp.g <- dplyr::left_join(lst.sp.g, groups, c("LB_NOM_VALIDE_SPE_LEVEL_SYN" = "species"))
openxlsx::write.xlsx(lst.sp.g, here::here(paste0('outputs/FinalList-of-SNAP-Vertebrate-Species_ActT-Diet-ForagS-NestH-Morpho-HabPref-DispD-MoveMod-LifeHist-PressureTraits_',g,'_GroupID_K=',
                                                 k, '_Weight-',m,'.xlsx')))
saveRDS(dendro$clust.dendrograms[[1]], here::here(paste0('outputs/HierarchicalClustering_Dendrogram_', g,'_Weight-', m)))

combi.doable <- rbind(combi.doable, data.frame(group = g, Nclus = k, mode = m))


################################################################
################### AVES ###################
################################################################

g = 'Aves'

############## Traits, env equal ############## 
m = 'EqualTraitsEnv'
mat.dist <- readRDS(here::here(paste0('data/derived-data/DissimilarityMatrices/SNAP-Vertebrate-Species_PairwiseCombinedDissimilarity_TraitsEnv_',g,'_Weight-', m)))
dendro <- PRE_FATE.speciesClustering_step1(mat.species.DIST = list(mat.dist), opt.no_clust_max = 25)
# Get tree
dt.tree <- ape::as.phylo(dendro$clust.dendrograms[[1]])
plot(dt.tree, type = "unrooted", no.margin = TRUE)

# Get groups
k = 19 # 13 et 19 mais pas facile
groups <- cutree(dendro$clust.dendrograms[[1]], k = k)
# Plot
# colors <- sample(fullcolor, length(unique(groups)))
colors <- viridisLite::turbo(length(unique(groups)))
colors <-  colorspace::darken(colors, amount = 0.2)
plot(dt.tree, type = "cladogram", no.margin = TRUE, tip.color = colors[groups], cex = 0.8)
# Groups in a DF 
groups <- data.frame(species = names(groups), cluster.id = c(groups))
groups$species <- gsub('_', ' ', groups$species)
lst.sp <- openxlsx::read.xlsx(here::here('data/derived-data/FinalList-of-SNAP-Vertebrate-Species_ActT-Diet-ForagS-NestH-Morpho-HabPref-DispD-MoveMod-LifeHist-PressureTraits.xlsx'))
lst.sp.g <- lst.sp[lst.sp$GROUP2_INPN == 'Oiseaux',]
lst.sp.g <- dplyr::left_join(lst.sp.g, groups, c("LB_NOM_VALIDE_SPE_LEVEL_SYN" = "species"))
openxlsx::write.xlsx(lst.sp.g, here::here(paste0('outputs/FinalList-of-SNAP-Vertebrate-Species_ActT-Diet-ForagS-NestH-Morpho-HabPref-DispD-MoveMod-LifeHist-PressureTraits_',g,'_GroupID_K=',
                                                 k, '_Weight-',m,'.xlsx')))
saveRDS(dendro$clust.dendrograms[[1]], here::here(paste0('outputs/HierarchicalClustering_Dendrogram_', g,'_Weight-', m)))

combi.doable <- rbind(combi.doable, data.frame(group = g, Nclus = k, mode = m))

############## Traits, env half half ############## 
m = 'HalfTraits-HalfEnv'
mat.dist <- readRDS(here::here(paste0('data/derived-data/DissimilarityMatrices/SNAP-Vertebrate-Species_PairwiseCombinedDissimilarity_TraitsEnv_',g,'_Weight-', m)))
dendro <- PRE_FATE.speciesClustering_step1(mat.species.DIST = list(mat.dist), opt.no_clust_max = 25)
# Get tree
dt.tree <- ape::as.phylo(dendro$clust.dendrograms[[1]])
plot(dt.tree, type = "unrooted", no.margin = TRUE)
# Get groups
k = 16 #11 et 16 
groups <- cutree(dendro$clust.dendrograms[[1]], k = k)
# Plot
colors <- viridisLite::turbo(length(unique(groups)))
colors <-  colorspace::darken(colors, amount = 0.2)
plot(dt.tree, type = "cladogram", no.margin = TRUE, tip.color = colors[groups], cex = 0.8)
# Groups in a DF 
groups <- data.frame(species = names(groups), cluster.id = c(groups))
# groups <- PRE_FATE.speciesClustering_step2(clust.dendrograms = dendro$clust.dendrograms, mat.species.DIST = mat.dist, no.clusters = 13)
groups$species <- gsub('_', ' ', groups$species)
lst.sp <- openxlsx::read.xlsx(here::here('data/derived-data/FinalList-of-SNAP-Vertebrate-Species_ActT-Diet-ForagS-NestH-Morpho-HabPref-DispD-MoveMod-LifeHist-PressureTraits.xlsx'))
lst.sp.g <- lst.sp[lst.sp$GROUP2_INPN == 'Oiseaux',]
lst.sp.g <- dplyr::left_join(lst.sp.g, groups, c("LB_NOM_VALIDE_SPE_LEVEL_SYN" = "species"))
openxlsx::write.xlsx(lst.sp.g, here::here(paste0('outputs/FinalList-of-SNAP-Vertebrate-Species_ActT-Diet-ForagS-NestH-Morpho-HabPref-DispD-MoveMod-LifeHist-PressureTraits_',g,'_GroupID_K=',
                                                 k, '_Weight-',m,'.xlsx')))
saveRDS(dendro$clust.dendrograms[[1]], here::here(paste0('outputs/HierarchicalClustering_Dendrogram_', g,'_Weight-', m)))

combi.doable <- rbind(combi.doable, data.frame(group = g, Nclus = k, mode = m))

################### AMPHIBIANS - REPTILES  ###################
g = 'Anura-Squamata-Testudines-Urodela'

############## Traits, env equal ############## 
m = 'EqualTraitsEnv'
mat.dist <- readRDS(here::here(paste0('data/derived-data/DissimilarityMatrices/SNAP-Vertebrate-Species_PairwiseCombinedDissimilarity_TraitsEnv_',g,'_Weight-', m)))
dendro <- PRE_FATE.speciesClustering_step1(mat.species.DIST = list(mat.dist), opt.no_clust_max = 15)
# Get tree
dt.tree <- ape::as.phylo(dendro$clust.dendrograms[[1]])
plot(dt.tree, type = "unrooted", no.margin = TRUE)

# Get groups
k = 11 # 7 ou 11
groups <- cutree(dendro$clust.dendrograms[[1]], k = k)
# Plot
# colors <- sample(fullcolor, length(unique(groups)))
colors <- viridisLite::turbo(length(unique(groups)))
colors <-  colorspace::darken(colors, amount = 0.2)
plot(dt.tree, type = "cladogram", no.margin = TRUE, tip.color = colors[groups], cex = 0.8)
# Groups in a DF 
groups <- data.frame(species = names(groups), cluster.id = c(groups))
groups$species <- gsub('_', ' ', groups$species)
lst.sp <- openxlsx::read.xlsx(here::here('data/derived-data/FinalList-of-SNAP-Vertebrate-Species_ActT-Diet-ForagS-NestH-Morpho-HabPref-DispD-MoveMod-LifeHist-PressureTraits.xlsx'))
lst.sp.g <- lst.sp[lst.sp$GROUP2_INPN == 'Amphibiens' | lst.sp$GROUP2_INPN == 'Reptiles'  ,]
lst.sp.g <- dplyr::left_join(lst.sp.g, groups, c("LB_NOM_VALIDE_SPE_LEVEL_SYN" = "species"))
openxlsx::write.xlsx(lst.sp.g, here::here(paste0('outputs/FinalList-of-SNAP-Vertebrate-Species_ActT-Diet-ForagS-NestH-Morpho-HabPref-DispD-MoveMod-LifeHist-PressureTraits_',g,'_GroupID_K=',
                                                 k, '_Weight-',m,'.xlsx')))
saveRDS(dendro$clust.dendrograms[[1]], here::here(paste0('outputs/HierarchicalClustering_Dendrogram_', g,'_Weight-', m)))

combi.doable <- rbind(combi.doable, data.frame(group = g, Nclus = k, mode = m))

############## Traits, env half half ############## 
m = 'HalfTraits-HalfEnv'
mat.dist <- readRDS(here::here(paste0('data/derived-data/DissimilarityMatrices/SNAP-Vertebrate-Species_PairwiseCombinedDissimilarity_TraitsEnv_',g,'_Weight-', m)))
dendro <- PRE_FATE.speciesClustering_step1(mat.species.DIST = list(mat.dist), opt.no_clust_max = 15)
# Get tree
dt.tree <- ape::as.phylo(dendro$clust.dendrograms[[1]])
plot(dt.tree, type = "unrooted", no.margin = TRUE)
# Get groups
k = 7 # 3 et 7 ok ; 10 possible mais espèces seules 
groups <- cutree(dendro$clust.dendrograms[[1]], k = k)
# Plot
colors <- viridisLite::turbo(length(unique(groups)))
colors <-  colorspace::darken(colors, amount = 0.2)
plot(dt.tree, type = "cladogram", no.margin = TRUE, tip.color = colors[groups], cex = 0.8)
# Groups in a DF 
groups <- data.frame(species = names(groups), cluster.id = c(groups))
# groups <- PRE_FATE.speciesClustering_step2(clust.dendrograms = dendro$clust.dendrograms, mat.species.DIST = mat.dist, no.clusters = 13)
groups$species <- gsub('_', ' ', groups$species)
lst.sp <- openxlsx::read.xlsx(here::here('data/derived-data/FinalList-of-SNAP-Vertebrate-Species_ActT-Diet-ForagS-NestH-Morpho-HabPref-DispD-MoveMod-LifeHist-PressureTraits.xlsx'))
lst.sp.g <- lst.sp[lst.sp$GROUP2_INPN == 'Amphibiens' | lst.sp$GROUP2_INPN == 'Reptiles'  ,]
lst.sp.g <- dplyr::left_join(lst.sp.g, groups, c("LB_NOM_VALIDE_SPE_LEVEL_SYN" = "species"))
openxlsx::write.xlsx(lst.sp.g, here::here(paste0('outputs/FinalList-of-SNAP-Vertebrate-Species_ActT-Diet-ForagS-NestH-Morpho-HabPref-DispD-MoveMod-LifeHist-PressureTraits_',g,'_GroupID_K=',
                                                 k, '_Weight-',m,'.xlsx')))
saveRDS(dendro$clust.dendrograms[[1]], here::here(paste0('outputs/HierarchicalClustering_Dendrogram_', g,'_Weight-', m)))

combi.doable <- rbind(combi.doable, data.frame(group = g, Nclus = k, mode = m))

openxlsx::write.xlsx(combi.doable, here::here('outputs/List-of-clustering-schemes.xlsx'))

################################
# Get group summary
################################
combi.doable <- openxlsx::read.xlsx(here::here('outputs/List-of-clustering-schemes.xlsx'))

for (i in 1:nrow(combi.doable)) {
  g <- combi.doable$group[i]
  k <- combi.doable$Nclus[i]
  w <- combi.doable$mode[i]
  rmarkdown::render(here::here('outputs/SummaryFG.rmd'), output_format = 'pdf_document', 
                        output_file = sprintf("FunctionalGroups-for_%s_K=%s_Weight-%s.pdf", g, k, w), 
                        params = list(group = g, Nclus = k, mode = w))
}

