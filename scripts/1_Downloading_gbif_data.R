###########building ocurances for elevation analyses ######################################

library(spocc)
library(rgeos)
library(dplyr)
library(wallace)
library(elevatr)
library(maptools)
library(rgeos)
library(dplyr)
library(rworldmap)
library(BBmisc)
library(lattice)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(ggthemes)
source(system.file('shiny/funcs', 'functions.R', package = 'wallace'))

#newmap <- getMap(resolution = "coarse")  # different resolutions available
#predsProj <- stack("../1_niche_models/World_clim_var_res05.gri")

setwd("~/Dropbox (ecoevounifesp)/1_pos_DOC/desenvolvimento/1_paper_22sp/teste_correlacao_latitude_vs_altitude/")
sp_names <- read.csv("1_species_names.txt", header = F)
#sp <- as.character(sp_names[1,1])
#sp <- "Castanozoster thoracicus"
for(sp in sp_names$V1){ 
# query selected database for occurrence records
results <- spocc::occ(query = sp, from = "gbif", limit = 5000, has_coords = TRUE)
# retrieve data table from spocc object
results.data <- results[["gbif"]]$data[[formatSpName(sp)]]
# remove rows with duplicate coordinates
occs.dups <- duplicated(results.data[c('longitude', 'latitude')])
occs <- results.data[!occs.dups,]
# give all records a unique ID
occs$occID <- row.names(occs)
occs2 <- cbind(occs$name, occs$longitude, occs$latitude, occs$stateProvince, occs$locality, occs$year, occs$month, occs$identifier)
colnames(occs2) <- c("name", "longitude", "latitude", "state", "locality", "year", "month", "identifier")
write.csv(occs2, paste(sp, "_occurances.csv", sep = ""))
}

#####adding elevation
All_sp <- NULL

for(sp in sp_names$V1){ 
  
  spe <- read.csv(paste(sp, "_occurances.csv", sep = ""), header = T)
  
  All_sp <- rbind(All_sp, spe)
  
  
}

  d <- cbind(All_sp$longitude, All_sp$latitude)
  ll_prj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  a <- get_elev_point(locations = as.data.frame(d), src = "aws", prj = ll_prj)
  All_sp <- cbind(All_sp, as.data.frame(a$elevation))
 
  
  
  All_sp <- cbind(All_sp, tempv)
  
  write.csv(All_sp, "1_All_sp.csv")
  
  All_sp <- read.csv("1_All_sp.csv", header = T)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
   
   




