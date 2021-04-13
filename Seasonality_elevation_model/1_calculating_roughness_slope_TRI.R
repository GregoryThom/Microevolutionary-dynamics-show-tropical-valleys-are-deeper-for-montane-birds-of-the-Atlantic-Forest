library(raster)
library(gdalUtils)
library(rgdal)
source(system.file('shiny/funcs', 'functions.R', package = 'wallace'))
library(rworldmap)
setwd("~/Dropbox/1_pos_DOC/desenvolvimento/1_paper_22sp/Seasonality_elevation_model/")

#a <- grep(".tif", list.files("./"), value = T)

# e <- extent(-50, -45, -25, -21)
# template <- raster(e)
# projection(template) <- '+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
# writeRaster(template, file="MyBigNastyRasty.tif", format="GTiff")
# mosaic_rasters(gdalfile=a,dst_dataset="MyBigNastyRasty.tif",of="GTiff")
# gdalinfo("MyBigNastyRasty.tif")

#x <- raster("MyBigNastyRasty.tif")
x <- raster("SRTM_W_250m_TIF/SRTM_W_250m.tif")
#plot(x)

occurrences <- read.csv("1_coordinates_gbif.csv",header=T)[1:3]
occs <- occurrences

occur_sp <- occs[2:3]
points <- SpatialPoints(occur_sp)

occs.xy <- occur_sp[c('Longitude', 'Latitude')]
sp::coordinates(occs.xy) <- ~ Longitude + Latitude
crs(occs.xy) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
# plot(predsProj[[1]])
# plot(occs.xy, add=T)
#plot(pnewmap, add=T)
bgExt <- mcp(occs.xy)
bgExt <- rgeos::gBuffer(bgExt, width = 0.0)
predsProj <- raster::crop(x, bgExt)
predsProj <- raster::mask(predsProj, bgExt)

pnewmap <- getMap(resolution = "coarse")  # different resolutions available

pdf("1_elevation.pdf")
plot(predsProj, yaxt="n", xaxt="n" )
axis(2,cex.axis=2)
axis(1,cex.axis=2)

plot(pnewmap, add=T)
dev.off()

y <- terrain(predsProj, opt = "slope")
pdf("1_slope.pdf")
plot(y, yaxt="n", xaxt="n" )
axis(2,cex.axis=2)
axis(1,cex.axis=2)
plot(pnewmap, add=T)
dev.off()

w <- terrain(predsProj, opt = "roughness")
pdf("1_roughness.pdf")
plot(w, yaxt="n", xaxt="n" )
axis(2,cex.axis=2)
axis(1,cex.axis=2)
plot(pnewmap, add=T)
dev.off()

t <- terrain(predsProj, opt = "TRI")
pdf("1_TRI.pdf")
plot(t, yaxt="n", xaxt="n" )
axis(2,cex.axis=2)
axis(1,cex.axis=2)
plot(pnewmap, add=T)
dev.off()



seas <- raster("CHELSA_bio10_04.tiff")
seas <- raster::crop(seas, bgExt)
seas <- raster::mask(seas, bgExt)

#plot(seas)

occur_sp$elev <- extract(predsProj, occur_sp[,c("Longitude", "Latitude")])
occur_sp$slop <- extract(y, occur_sp[,c("Longitude", "Latitude")])
occur_sp$rough <- extract(w, occur_sp[,c("Longitude", "Latitude")])
occur_sp$TRI <- extract(t, occur_sp[,c("Longitude", "Latitude")])
occur_sp$seas <- extract(seas, occur_sp[,c("Longitude", "Latitude")])


m1 <- lm(rough ~ elev * Latitude, data = occur_sp)
m1_pred <- predict(m1, newdata = occur_sp, interval = "confidence")
occur_sp$y_hat <- m1_pred[,1]
occur_sp$y_lwr <- m1_pred[,2]
occur_sp$y_upr <- m1_pred[,3]

# plot
ggplot(data=occur_sp, aes(seas, rough)) + 
  geom_point(col="steelblue", size=2) + 
  geom_line(aes(seas, y_hat), col="red") +
  geom_ribbon(aes(ymin=y_lwr, ymax=y_upr), fill="magenta", alpha=.25) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


