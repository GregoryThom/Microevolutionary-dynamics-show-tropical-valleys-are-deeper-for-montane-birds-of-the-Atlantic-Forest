#############ResistenceGA#########################
library(ResistanceGA)
library(gdistance)
library(raster)
library(spatialEco)
library(wallace)
source(system.file('shiny/funcs', 'functions.R', package = 'wallace'))

setwd("../inputs_and_data_sets/7_ResistanceGA/")

files <- read.csv("../species.txt", header = F)
#file <- files
for(file in files){ 
  
  ##produce .asc for resistance (current niche model), seasonality and stability(mean probability of occurance between all time slices after aplying threshold)
  dgen <- as.dist(read.csv(paste(file, ".gen_dist", sep = ""), sep = " "))
  res <- raster(paste(file, "/resistance.asc", sep = ""))
 # saz <- raster(paste(file, "/sazonality.asc", sep = ""))
  r1 <- raster::getData(name = "worldclim", var = "bio", res = 2.5)
  saz <- r1[[4]]/100
  sta <- raster(paste(file, "/stability.asc", sep = ""))

  dir.create(file.path(file, file))
  dir.create(file.path(file, "all.comb"))
  write.dir <- paste("7_ResistanceGA/", file,"/",file, "/",sep = "")

  coords <- read.csv(paste(file,"/", file, ".coord",  sep = ""), header = F, sep = " ")
  write.table(coords,file=paste0(write.dir,"../",file,"_samples.txt"),sep="\t",col.names=F,row.names=F)
  
 sample.locales <- SpatialPoints(coords[,2:3])
#plot(res)
#plot(sample.locales, add=T)
#plot(saz)
#plot(sample.locales, add=T)

bgExt <- mcp(sample.locales)

bgExt <- rgeos::gBuffer(bgExt, width = 1.0)
#plot(bgExt, add=T)
#crop the environmental rasters by the background extent shape
envsBgCrop <- raster::crop(res, bgExt)
envsBgCrop_lr <- aggregate(envsBgCrop, fact=10)
crs(envsBgCrop_lr) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

envsBgCrop2 <- raster::crop(saz, envsBgCrop_lr, snap="near")
envsBgCrop2_lr <- aggregate(envsBgCrop2, fact=5)
#envsBgCrop2_lr <- raster.transformation(envsBgCrop2_lr, trans = "norm", smin = 0, smax = 1)
#projCoords <- data.frame(x = c(-51.61681,-51.61681, -43.36681, -43.36681,-51.61681), y = c(-30.49181, -30.49181,-21.57514,-21.57514,-30.49181))
envsBgCrop3 <- raster::crop(sta, envsBgCrop_lr, snap="near")
envsBgCrop3_lr <- aggregate(envsBgCrop3, fact=5)
crs(envsBgCrop3_lr) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"



envsBgCrop2_lr <- projectRaster(envsBgCrop2_lr, envsBgCrop_lr)
envsBgCrop3_lr <- projectRaster(envsBgCrop3_lr, envsBgCrop_lr)


writeRaster(envsBgCrop_lr,
            filename = paste0(write.dir, file, "_res.asc"),
            overwrite = TRUE)

writeRaster(envsBgCrop2_lr,
            filename = paste0(write.dir, file, "_asaz.asc"),
            overwrite = TRUE)

writeRaster(envsBgCrop3_lr,
            filename = paste0(write.dir, file, "_sta.asc"),
            overwrite = TRUE)

gdist.inputs <- gdist.prep(n.Pops = length(sample.locales), 
                           response = as.vector(dgen), 
                           samples = sample.locales, 
                           directions = 8, 
                           longlat = T,
                           method = 'costDistance')
 
GA.inputs <- GA.prep(ASCII.dir = write.dir,
                     method = "LL",
                     Results.dir = "all.comb",
                     max.cat = 500,
                     max.cont = 2500,
                     run = 25,
                     select.trans = list("M", "M", "M"),
                     scale = T,
                     maxiter = 1000,
                     seed = 123,
                     parallel = 6)



#SS_RESULTS <- SS_optim(gdist.inputs=gdist.inputs,
#                       GA.inputs=GA.inputs)

ac <- all_comb(gdist.inputs=gdist.inputs,
               GA.inputs=GA.inputs,
               max.combination = 4,
               results.dir = paste(file, "/all.comb/",sep = ""),
               iters = 1000,
               replicate = 1,
               sample.prop = 0.9,
               nlm = F,
               dist_mod = TRUE,
               null_mod = T)


}
