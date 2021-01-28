######SDM's Projection to different times

library(maptools)
library(ENMeval)
library(rgeos)
library(spocc)
library(spThin)
library(dismo)
library(dplyr)
library(rworldmap)
library(beepr)
source(system.file('shiny/funcs', 'functions.R', package = 'wallace'))
#library(rJava) # required to run Maxent within dismo


setwd("~/Dropbox/1_pos_DOC/desenvolvimento/1_paper_22sp/1_niche_models_chelsa/")



raster_data <- list.files(path = "CHELSA_cur_V1_2B_r2_5m/", pattern = "tif", full.names = T)

predictors <- stack(raster_data)

#plot(predictors)

projCoords <- data.frame(x = c(-65, -65, -39, -39, -65), y = c(-15, -40, -40, -15, -15))
projPoly <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(projCoords)), ID=1)))
predsProj <- raster::crop(predictors, projPoly)
predsProj <- raster::mask(predsProj, projPoly)


pnewmap <- getMap(resolution = "coarse")  # different resolutions available


files <- list.dirs("../2_conStruct/1_no_spatial/", full.names=F, recursive=FALSE)
files <- files[21:22]
file <- files



for(file in files){ 
  
  

#getting occurrence coordinates
occurrences <- read.csv("1_All_sp.csv",header=T)[1:3]
occs <- occurrences[grep(file, occurrences$name), ][1:3]
#occs <- occs[1:180,]
#Spatial thinning selected. Thin distance selected is 20 km.
output <- spThin::thin(occs, 'latitude', 'longitude', 'name', thin.par = 8, reps = 10, locs.thinned.list.return = TRUE, write.files = FALSE, verbose = FALSE)
#Since spThin did 100 iterations, there are 100 different variations of how it thinned your occurrence localities. As there is a stochastic element in the algorithm, some iterations may include more localities than the others, and we need to make sure we maximize the number of localities we proceed with.
# find the iteration that returns the max number of occurrences
maxThin <- which(sapply(output, nrow) == max(sapply(output, nrow)))
# if there's more than one max, pick the first one
maxThin <- output[[ifelse(length(maxThin) > 1, maxThin[1], maxThin)]]  
# subset occs to match only thinned occs
occs <- occs[as.numeric(rownames(maxThin)),] 

if (nrow(occs) > 180) {occs <- occs[sample(nrow(occs), 180), ]

} else {occs <- occs

}

occur_sp <- occs[2:3]
points <- SpatialPoints(occur_sp)



#checking if any point is out of the raster grid
locs.vals <- raster::extract(predsProj[[1]], occur_sp)
# remove occs without environmental values
occur_sp <- occur_sp[!is.na(locs.vals), ]  

rm(occurrences, occs, maxThin)
#####################################
###     ENM DATA PREPARATION     ####
#####################################

occs.xy <- occur_sp[c('longitude', 'latitude')]
sp::coordinates(occs.xy) <- ~ longitude + latitude
crs(occs.xy) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
plot(predsProj[[1]])
plot(occs.xy, add=T)
plot(pnewmap, add=T)

bgExt <- mcp(occs.xy)
bgExt <- rgeos::gBuffer(bgExt, width = 0.5)
plot(bgExt, add=T)
#crop the environmental rasters by the background extent shape
envsBgCrop <- raster::crop(predsProj, bgExt)
#mask the background extent shape from the cropped raster
envsBgMsk <- raster::mask(envsBgCrop, bgExt)
crs(envsBgMsk) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
# sample random background points
bg.xy <- dismo::randomPoints(envsBgMsk, 10000)
# convert matrix output to data frame
bg.xy <- as.data.frame(bg.xy)  
bg.xy <- SpatialPoints(bg.xy)
crs(bg.xy) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

rm(envsBgCrop)
#plot(envsBgMsk)
#jar <- paste(system.file(package="dismo"), "/Java/maxent.jar", sep='')

#######################
#     ENMevaluate     #
#######################
#Jacknife can bu used when less than 20 samples are available
#group.data <- ENMeval::get.jackknife(occ=occs.xy, bg.coords=bg.xy)
group.data <- ENMeval::get.randomkfold(occ=occs.xy, bg.coords=bg.xy, kfolds=4)
# pull out the occurrence and background partition group numbers from the list
occs.grp <- group.data[[1]]
bg.grp <- group.data[[2]]
rms <- seq(1, 4, .5)
#beep()

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#source(system.file('shiny/funcs', 'functions.R', package = 'wallace'))

###primeiro testar qual o melhor modelo 
# iterate model building over all chosen parameter settings
rm(res)
res <- ENMeval::ENMevaluate(occs.xy, envsBgMsk, bg.coords = bg.xy, RMvalues = rms, fc = c('L', 'LQ', 'LQH'), 
                            method = 'user', occs.grp, bg.grp, clamp = T, algorithm = "maxnet", parallel = F) #, numCores = 2)
if(exists("res")) {beep()
} else  {beep(9)
    }
#beep()
#save(res, file= paste(file, "/Maxent_model_full_data_teste_", file, ".rda", sep=""))


evalTbl <- res@results
evalMods <- res@models
names(evalMods) <- res@results$settings
#evalPreds <- res@predictions
write.csv(evalTbl, file= paste(file,"/",file,"_ENMEval.csv", sep=""))
occs <- res@occ.pts
abs <- res@bg.pts
rm(res)
beep()
#}

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################


best_AUC <- evalTbl[order(evalTbl$avg.test.AUC, decreasing = T),][1:5,]
best_AIC <- as.character(best_AUC[order(best_AUC$delta.AICc, decreasing = F),][1,1])


mod <- evalMods[[best_AIC]]

pw <- ENMeval::maxnet.predictRaster(mod, envsBgMsk, type = 'logistic', clamp = TRUE)
plot(pw)

#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
##Projectiong to other times
#time<- "21k_2_5"
for (time in list.files(path='1_projection_past', full.names=F)) {
  

files_past <- list.files(path= paste('1_projection_past/', time, sep=""), pattern='tif', full.names=TRUE)
predictors_past <- stack(files_past)

#names(predictors_past) <- names(predictors)
predspast_Proj <- raster::crop(predictors_past, projPoly)
predspast_Proj <- raster::mask(predspast_Proj, projPoly)

proj <- ENMeval::maxnet.predictRaster(mod, predspast_Proj, type = 'logistic', clamp = T)
plot(proj)
plot(points, add=T)
plot(pnewmap, add=T)

writeRaster(proj, filename= paste(file,"/",file,"_projection_", time, "_AF_data.asc", sep=""),overwrite=TRUE)
pdf(paste(file,"/",file,"_projection_", time, "_AF_data", sep=""))
plot(proj)
plot(pnewmap, add=T)
dev.off()
# get predicted values for occurrence grid cells
occPredVals_proj <- raster::extract(proj, occs.xy)
# define minimum training presence threshold
source(system.file('shiny/funcs', 'functions.R', package = 'wallace'))
thr_proj <- thresh(occPredVals_proj, "p10")
# threshold model prediction
proj_p10 <- proj > thr_proj
# plot the model prediction
pdf(paste(file,"/",file, time,"_projection_AF_threshold_data", sep=""))
plot(proj_p10)
plot(pnewmap, add=T)
dev.off()
writeRaster(proj_p10, filename=paste(file,"/",file, time, "_projection_AF_threshold_data.asc", sep=""),overwrite=TRUE)
#beep()



############MESS####################

pts <- SpatialPoints(rbind(occs,abs))
occEnvVals <- raster::extract(predspast_Proj, pts)

# Threshold mess maps
###projection in space
#reference_points <- extract(envsBgMsk, occEnvVals)
mss <- mess(x=predspast_Proj, v=occEnvVals)
plot(mss)

mess_thresh <- mss
thresh <- 0
mess_thresh[mess_thresh<thresh] <- NA
mess_thresh[mess_thresh>=thresh] <- 1

mess_resamp <- resample(mess_thresh,predspast_Proj)
plot(mess_resamp)

proj_mess <- overlay(proj, mess_resamp, fun = function(x, y) {
  x[is.na(y[])] <- 0
  return(x)
})
plot(proj_mess)
pdf(paste(file,"/",file,time, "_projection_MESS", sep=""))
plot(proj_mess)
plot(pnewmap, add=T)
dev.off()
writeRaster(proj_mess, filename= paste(file,"/",file, time, "_projection_MESS.asc", sep=""),overwrite=TRUE)

proj_mess_p10 <- proj_mess > thr_proj
plot(proj_mess_p10)
plot(pnewmap, add=T)
pdf(paste(file,"/",file, time, "_projection_MESS_p10", sep=""))
plot(proj_mess_p10)
plot(pnewmap, add=T)
dev.off()
writeRaster(proj_mess_p10, filename= paste(file,"/",file, time, "_projection_MESS_p10.asc", sep=""),overwrite=TRUE)
#beep()

}
}


