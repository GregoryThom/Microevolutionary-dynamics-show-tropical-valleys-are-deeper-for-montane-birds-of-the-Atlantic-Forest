library(raster)
library(sp)
library(wordcloud)
library(geosphere)
library(hypervolume)
library(rworldmap)
library(maps)
library(BBmisc)
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

setwd("../inputs_and_data_sets/")


#First download the data from chelsa data base on the resolution you like "https://chelsa-climate.org/bioclim/"
raster_data <- list.files(path = "Chelsa_30s", pattern = "CHELSA", full.names = T)

predictors <- stack(raster_data)


#Here I'm cropping bioclim layers to the distribution of the Atlantic Forest
projCoords <- data.frame(x = c(-65, -65, -39, -39, -65), y = c(-15, -40, -40, -15, -15))
projPoly <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(projCoords)), ID=1)))
predsProj <- raster::crop(predictors, projPoly)
predsProj <- raster::mask(predsProj, projPoly)


pnewmap <- getMap(resolution = "coarse")  # different resolutions available



#r <- getData("worldclim",var="bio",res=0.5, lon=-54, lat=-18)
#r2 <- getData("worldclim",var="bio",res=0.5, lon=-54, lat=-35)
#r3 <- merge(r, r2)
data_total <- as.data.frame(na.omit(read.csv("1_coordinates_gbif.csv", sep = ",", header = T)[,c('name','longitude', 'latitude', 'pop', 'elevation')]))

climatelayers <- predsProj
climatelayers_n <- raster::scale(climatelayers)

# z-transform climate layers to make axes comparable
   #climatelayers_ss = climatelayers[[c(4,8:11,15:19)]]
#   for (i in 1:nlayers(climatelayers))
#   {
#   climatelayers_ss[[i]] <- (climatelayers_ss[[i]] - cellStats(climatelayers_ss[[i]], 'mean')) / cellStats(climatelayers_ss[[i]], 'sd') 
#   }
   
#   climatelayers_ss_cropped = crop(climatelayers_ss, extent(-60,-40,-37,-18))

   files <- unique(data_total$name)
   files <- files[-grep("Lepidocolaptes", files)]
   vol_all_sp <- matrix(nrow = 0, ncol = 3)
   pca_vari <- matrix(nrow = 1, ncol = 0)
   niche_overlap <- matrix(nrow = 6, ncol = 0)
   
   #file <- files[1]
for(file in files){
   
   data <- na.omit(subset(data_total, name==file))
   norte <- SpatialPoints(subset(data, pop=="north")[,c("longitude","latitude")])
   sul <- SpatialPoints(subset(data, pop=="South")[,c("longitude","latitude")])
   climate_norte = as.data.frame(na.omit(extract(climatelayers_n, norte)))
   climate_sul = as.data.frame(na.omit(extract(climatelayers_n, sul)))
   
   
   
#climate_norte <- extract(r3, norte)[,1:19]
#climate_sul <- extract(r3, sul)[,1:19]
PCA_norte <- prcomp(na.omit(climate_norte))
sum <- summary(PCA_norte) 
sum<- as.data.frame(sum$importance[3,4])
colnames(sum)<- paste(file,"N")
pca_vari <- cbind(pca_vari, sum)

PCA_sul <- prcomp(na.omit(climate_sul))
sum <- summary(PCA_sul) 
sum<- as.data.frame(sum$importance[3,4])
colnames(sum)<- paste(file,"S")
pca_vari <- cbind(pca_vari, sum)


PCA_plot <- as.data.frame(PCA_norte$x)
PCA_plot$group <- row.names(PCA_norte$x)
ggplot(PCA_plot, aes(x=PCA_plot$PC1, y=PCA_plot$PC2)) + geom_point(shape=21)

hv_norte <- hypervolume_gaussian(PCA_norte$x[,1:4], name = 'norte', samples.per.point = 5)
hv_sul <- hypervolume_gaussian(PCA_sul$x[,1:4], name = 'sul', samples.per.point = 5)

#hv_norte <- hypervolume_gaussian(climate_norte[,1:4], name = 'norte', samples.per.point = 5)
#hv_sul <- hypervolume_gaussian(climate_sul[,1:4], name = 'sul', samples.per.point = 5)


hv_set = hypervolume_set(hv_norte, hv_sul, check.memory=FALSE)
volumes <- get_volume(hv_set)
vol <- as.data.frame(volumes)
colnames(vol) <- file
niche_overlap <-cbind(niche_overlap, vol)
vol_taxa <- as.data.frame(as.matrix(volumes)[1:2,])
vol_taxa$sp <- file
vol_taxa$pop <- rownames(vol_taxa)
vol_all_sp <- rbind(vol_all_sp, vol_taxa)




}
   
   pca_vari
   niche_overlap
   vol_all_sp
   
write.csv(niche_overlap, "Niche_hypervolume_all_info.csv")
write.csv(pca_vari, "PCA_Niche_hypervolume.csv")

#vol_all_sp2 <- vol_all_sp[-grep("Phylloscartes difficilis", vol_all_sp$sp),]
#vol_all_sp2$`as.matrix(volumes)[1:2, ]` <- log(vol_all_sp2$`as.matrix(volumes)[1:2, ]`)
vol_all_sp2 <- vol_all_sp
vol_all_sp2$`as.matrix(volumes)[1:2, ]` <- log(vol_all_sp2$`as.matrix(volumes)[1:2, ]`)
write.csv(vol_all_sp2, "log_Niche_hypervolume.txt")
   pdf("1_hypervolume_niche.pdf")
   ggplot(vol_all_sp2, aes(x = vol_all_sp2$pop, y = vol_all_sp2$`as.matrix(volumes)[1:2, ]`)) + geom_boxplot(fill = "gray") + geom_point(size = 1.5) +
     geom_line(aes(group = sp ), alpha = 0.4, colour = "black", linetype = "3313", data = vol_all_sp2)+ 
     labs(title="", 
          x="",
          y="log Niche hypervolume") + theme(legend.title = element_blank()) +
     theme(text = element_text(size=15), axis.text.x = element_text(angle=0, hjust=1)) +
     stat_compare_means(label.y = 2.5, size = 4)
   dev.off()    

op=par(mar=c(3,10,1,1))
barplot(volumes,horiz=TRUE,las=2,main="Hypervolume",cex.names=0.5,col='lightblue')

par(op)
plot(hv_set[[c(3,5,6)]])

#########################################################################################################

norte_map = hypervolume_project(hv_norte, climatelayers[,1,4],reduction.factor=0.1)
sul_map = hypervolume_project(hv_sul, climatelayers,reduction.factor=0.1)


#hyper <- hypervolume_gaussian(PCA_plot[,1:4])
#plot.Hypervolume(hyper, show.3d=T)

plot(norte_map,col=colorRampPalette(c(rgb(1,1,1),rgb(1,0,0)))(100),legend=T,main='norte')
map('world',add=TRUE)
points(latitude~longitude,data=norte,pch=16,cex=0.4)

plot(sul_map,col=colorRampPalette(c(rgb(1,1,1),rgb(0,0,1)))(100),legend=T, main='sul')
map('world',add=TRUE)
points(latitude~longitude,data=sul,pch=16,cex=0.4)
dev.off()

