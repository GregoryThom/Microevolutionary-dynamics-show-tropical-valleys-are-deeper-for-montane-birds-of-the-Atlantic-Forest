library(ips)
library(parallel)
library(lme4)
library(MuMIn)
library(nlme)
library(geiger)
setwd("~/Dropbox/1_pos_DOC/desenvolvimento/1_paper_22sp/Seasonality_elevation_model/")
random_points <- read.csv("1_random_points_across_AF.csv")
occurances <- read.csv("1_coordinates_gbif.csv")


###extracting seazonality and elevation for latitudinal bands

min_lat <- as.integer(max(random_points$latitude))
max_lat <- as.integer(min(random_points$latitude))
lats <- c(max_lat:min_lat)
#lat=-20

lat_bands <- function(lat) {
  tab <- random_points[grep(lat, random_points$latitude),]
  dat <- data.frame(matrix(ncol = 3, nrow = 1))
  x <- c("lat", "saz", "elev")
  colnames(dat) <- x
  dat$lat <- lat
  dat$saz <- mean(tab$CHELSA_bio10_04, na.rm = T)/100
  dat$elev <- c(max(tab$a.elevation)-min(tab$a.elevation))
  return(dat)
  }

system.time({
  randon_table <- mclapply(lats, lat_bands, mc.cores = 1)
})
randon_table <- do.call(rbind.data.frame, randon_table)
randon_table <- randon_table[order(randon_table$lat, decreasing = T),]


#####extracting species elevation range per latitude, removing the top and botton 2.5%
species <- unique(occurances$Species)
#sp=species[17]
#i=-29
elev_sp <- function(sp) {
  occs <- occurances[grep(sp, occurances$Species),]
  sp_lat=NULL
  for (i in lats) {
    tab <- occs[grep(i, occs$Latitude),]
    if (nrow(tab) > 30) {  
    dat <- data.frame(matrix(ncol = 4, nrow = 1))
    x <- c("lat", "min_elev", "max_elev", "range")
    colnames(dat) <- x
    dat$lat <- i
    dat$min_elev <- min(tab$Elevation[tab$Elevation > quantile(tab$Elevation, 0.05, na.rm = T)],na.rm = T)
    dat$max_elev <- max(tab$Elevation[tab$Elevation < quantile(tab$Elevation, 0.95, na.rm = T)],na.rm = T)
    dat$range <- dat$max_elev[1]-dat$min_elev[1]
    sp_lat <- rbind(sp_lat, dat)
    }
  }
  sp_lat$sp <- sp
  
  if (nrow(sp_lat)==0) {
    return(NULL)
  } else { if (nrow(sp_lat) > 4) {
    return(sp_lat) 
  } else {
    return(NULL)
  }
  }
    
}



system.time({
  occs_table <- mclapply(species, elev_sp, mc.cores = 1)
})
occs_table <- do.call(rbind.data.frame, occs_table)



####Combining tables
occs_table$saz <- NA
occs_table$elev_avai <- NA

for (i in 1:nrow(occs_table)) {
  for (x in 1:nrow(randon_table)) {
      if (randon_table$lat[x] == occs_table$lat[i]) {
      occs_table$saz[i] <- randon_table$saz[x]
      occs_table$elev_avai[i] <- randon_table$elev[x]
    }
    }
}


occs_table$occup_prop <- c(occs_table$range/occs_table$elev_avai)*100

write.csv(occs_table,"1_main_table.csv")
occs_table2 <- occs_table

occs_table2$saz <- sqrt(occs_table2$saz)
occs_table2$saz <- scale(occs_table2$saz, scale=T)

occs_table2$elev_avai <- sqrt(occs_table2$elev_avai)
occs_table2$elev_avai <- scale(occs_table2$elev_avai, scale=T)


library(nlme)
# #models occupied elevation GLM
# full <- glm(occup_prop ~ 1 + sp + elev_avai + lat, data=occs_table)
# lat <- glm(occup_prop ~ 1 + sp + lat, data=occs_table)
# elev <- glm(occup_prop ~ 1 + sp + elev_avai, data=occs_table)

#models occupied elevation Mixed model LME
full_prop <- lme(occup_prop ~ 1 + elev_avai * saz, random = ~1|sp, data=occs_table2)
saz_prop <- lme(occup_prop ~ 1 + saz, random = ~1|sp, data=occs_table2)
elev_prop <- lme(occup_prop ~ 1 + elev_avai, random = ~1|sp, data=occs_table2)

full_prop <- update(full_prop, method="ML")
saz_prop <- update(saz_prop, method="ML")
elev_prop <- update(elev_prop, method="ML")

ac <- AIC(full_prop, saz_prop, elev_prop)
ac <- ac[order(ac$AIC),]
aicw(as.numeric(ac[[2]]))

anova(full_prop)
summary(full_prop)
r.squaredGLMM(full_prop)

anova(full_prop, saz_prop)


#models elevation range Mixed model LME
full_range <- lme(range ~ 1 + elev_avai * saz, random = ~1|sp, data=occs_table2)
saz_range <- lme(range ~ 1 + saz, random = ~1|sp, data=occs_table2)
elev_range <- lme(range ~ 1 + elev_avai, random = ~1|sp, data=occs_table2)

full_range <- update(full_range, method="ML")
saz_range <- update(saz_range, method="ML")
elev_range <- update(elev_range, method="ML")

ac <- AIC(full_range, saz_range, elev_range)
ac <- ac[order(ac$AIC),]
ac
aicw(as.numeric(ac[[2]]))
anova(full_range)
summary(full_range)
r.squaredGLMM(full_range)

anova(full_range, elev_range)

#models minimum elevation Mixed model LME
full_min <- lme(min_elev ~ 1 + elev_avai * saz, random = ~1|sp, data=occs_table2)
saz_min <- lme(min_elev ~ 1 + saz, random = ~1|sp, data=occs_table2)
elev_min <- lme(min_elev ~ 1 + elev_avai, random = ~1|sp, data=occs_table2)

full_min <- update(full_min, method="ML")
saz_min <- update(saz_min, method="ML")
elev_min <- update(elev_min, method="ML")

ac <- AIC(full_min, saz_min, elev_min)
ac <- ac[order(ac$AIC),]
ac
aicw(as.numeric(ac[[2]]))
anova(saz_min)
summary(saz_min)
r.squaredGLMM(saz_min)

anova(full_min, saz_min)

#models maximum elevation Mixed model LME
full_max <- lme(max_elev ~ 1 + elev_avai * saz, random = ~1|sp, data=occs_table2)
saz_max <- lme(max_elev ~ 1 + saz, random = ~1|sp, data=occs_table2)
elev_max <- lme(max_elev ~ 1 + elev_avai, random = ~1|sp, data=occs_table2)

full_max <- update(full_max, method="ML")
saz_max <- update(saz_max, method="ML")
elev_max <- update(elev_max, method="ML")

ac <- AIC(full_max, saz_max, elev_max)
ac <- ac[order(ac$AIC),]
ac
aicw(as.numeric(ac[[2]]))
anova(elev_max)
summary(elev_max)
r.squaredGLMM(elev_max)

anova(full_max,elev_max)




