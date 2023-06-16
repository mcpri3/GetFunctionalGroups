get.mean.slope <- function(df) {
  return(data.frame(ID = unique(df$ID), mean = mean(na.omit(df$slope))))
}
get.mean.alti <- function(df) {
  return(data.frame(ID = unique(df$ID), mean = mean(na.omit(df$france_mnt_2154))))
}
get.sd.alti <- function(df) {
  return(data.frame(ID = unique(df$ID), sd = sd(na.omit(df$france_mnt_2154))))
}

get.slope.mean <- function(df) {
  library(dplyr)
  val <- terra::extract(x = slope, y = grid.fr[grid.fr$groupID == unique(df$groupID),])
  val <- group_by(val, ID) %>% do(get.mean.slope(.))
  return(data.frame(polyID =  grid.fr$polyID[grid.fr$groupID == unique(df$groupID)], slope.mean = val$mean))
}

get.alti.cv <- function(df) {
  library(dplyr)
  val <- terra::extract(x = alti, y = grid.fr[grid.fr$groupID == unique(df$groupID),])
  
  val.mn <- group_by(val, ID) %>% do(get.mean.alti(.))
  val.sd <- group_by(val, ID) %>% do(get.sd.alti(.))
  
  val <- left_join(val.mn, val.sd, by = 'ID')
  val$cv <- val$sd/val$mean*100
  return(data.frame(polyID =  grid.fr$polyID[grid.fr$groupID == unique(df$groupID)], alti.cv = val$cv))
}