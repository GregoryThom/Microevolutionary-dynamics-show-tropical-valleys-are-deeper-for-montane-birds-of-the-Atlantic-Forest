library(ape)
library(geiger)
library(caper)
library(Hmisc)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(reshape2)
library(BBmisc)
library(harrietr)
library(caret)
library(broom)
library(motmot)
library(adiv)


blomberg.k <- function(phy, y) {
  bm.phy <- transformPhylo.ML(phy, model="bm", y=y )
  ntip <- Ntip(phy)
  nodetimes <- nodeTimes(phy)
  vcv.int <- vcv(phy)
  vcv.phy <- sum(solve(vcv.int))
  diag.vcv <- sum(diag(vcv.int))
  mean.sq.error.0  <- bm.phy$brownianVariance[[1]]
  mat.state <- y[,1] - bm.phy$root.state
  transpose.mat <- t(mat.state)
  mean.sq.error <- transpose.mat %*% mat.state / (ntip - 1)
  mean.sq.error.div <- 1 / (ntip - 1) * ((diag.vcv - ntip / vcv.phy))
  
  obs <- mean.sq.error / mean.sq.error.0
  exp <- mean.sq.error.div
  k.out <- obs / exp
  return(k.out)
}

# Function to extract the overall ANOVA p-value out of a linear model object
lmp <- function (modelobject) {
  if (class(modelobject) != "pgls") stop("Not an object of class 'pgls' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}



setwd("~/Dropbox/1_pos_DOC/desenvolvimento/1_paper_22sp/8_PGLS/")
datos<-read.csv("1_PGLS_table.csv",header=TRUE)
arbol<-read.tree("1_42pop_tree.tre")
#obj<-name.check(arbol,datos)
#obj
datos_trans <- datos

datos_trans$kipps <- log(datos_trans$kipps)
datos_trans$Effective.pop.size <- log(datos_trans$Effective.pop.size)
datos_trans$Time.size.change <- log(datos_trans$Time.size.change)
datos_trans$Niche_volume <- log(datos_trans$Niche_volume)
datos_trans$Geographic.distance <- sqrt(datos_trans$Geographic.distance)



for(var in c(5, 7:ncol(datos_trans))){
datos_trans[[var]] <- normalize(datos_trans[[var]], method = "standardize", range = c(0, 1), margin = 2)
}


####Blomberg's K
data_blom <- datos_trans[5:ncol(datos_trans)]
#data_blom <- datos_trans[5:8]
rownames(data_blom) <- datos_trans$taxon
data_blom <- na.omit(data_blom)
y=as.matrix(data_blom)
phy=arbol
obj<-name.check(phy,y)
phy <- drop.tip(phy, obj$tree_not_data)

k = data.frame(matrix(vector(), 0, 2,
                       dimnames=list(c(), c("Blomberg's_K", "p-value"))),
                stringsAsFactors=F)




#i=1
for (i in 1:ncol(y)) {
  yb <- as.matrix(y[,i])
  blomK <- adiv::K(phy, yb, nrep = 999, alter = c("greater", "less", "two-sided"))
k2 <- cbind(as.data.frame(blomK$obs), as.data.frame(blomK$pvalue))
  colnames(k2) <- c("Blomberg's_K", "p-value")
  row.names(k2) <- colnames(y)[i]
  k <- rbind(k, k2)
  
}

k 

write.csv(k, "1_BlombergK.csv")

#ggqqplot(datos_trans$TDIV)
#ggdensity(datos_trans$TDIV, 
#          main = "Density plot of tooth length",
#          xlab = "Tooth length")
#shapiro.test(datos_trans$TDIV)


#write.csv(datos_trans, "1_PGLS_table_sqrt.csv")
#obj<-name.check(arbol[,1],datos)
#obj
###Spearmans correlation index
datos_res <- as.matrix(datos_trans[5:ncol(datos_trans)])
sper_cor <- rcorr(datos_res, type = "spearman")
sper_cor

library(lattice)
library(RColorBrewer)
breaks <- seq(-1, 1.2, by=0.2)
#breaks2 <- seq(-1, 102, by=5.0)
myColours2 <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
myColours3 <- colorRampPalette(brewer.pal(11,'Spectral'))

data0.1 <- round(t(as.matrix(sper_cor$r)), 2)
myPanel0.1 <- function(x, y, z, ...) {
  panel.levelplot(x,y,z,...)
  panel.text(x, y, data0.1[cbind(x,y)], cex=0.55) ## use handy matrix indexing
}

data0.2 <- round(t(as.matrix(sper_cor$P)), 2)
myPanel0.2 <- function(x, y, z, ...) {
  panel.levelplot(x,y,z,...)
  panel.text(x, y, data0.2[cbind(x,y)], cex=0.55) ## use handy matrix indexing
}

ht <- levelplot(t(as.matrix(sper_cor$r)), col.regions = myColours2, at=breaks, 
                panel=myPanel0.1, colorkey=list(labels=list(cex=1)), scales=list(x=list(cex=0),y=list(cex=1)),
                xlab=list("", cex=2), ylab=list("Predictors", cex=1.5), cex.lab=1)

htp <- levelplot(t(as.matrix(sper_cor$P)), col.regions = myColours3, at=c(0, 0.01, 0.05,1), 
                panel=myPanel0.2, colorkey=F, scales=list(x=list(cex=0),y=list(cex=1)),
                xlab=list("", cex=2), ylab=list("Predictors", cex=1.5), cex.lab=1)
ht
htp
pdf("../1_Figures/1_heatmap_spearman_correlation_PGLS_predictors.pdf")
ht
dev.off()
pdf("../1_Figures/1_heatmap_spearman_correlation_PGLS_predictors_Pvalue.pdf")
htp
dev.off()

###Testing if response variables are associated to mountain block (North vs South)

comp.data<-comparative.data(arbol, datos_trans, names.col="taxon", vcv.dim=2, warn.dropped=F)

table <- vector(mode = "list", length = ncol(datos_trans[8:19]))
table2 <- vector(mode = "list", length = ncol(datos_trans[8:19]))
names(table) <- colnames(datos_trans)[8:19]
names(table2) <- colnames(datos_trans)[8:19]
#table <-data.frame(matrix(ncol = 0, nrow = 5))
table3_all <- matrix(nrow = 1, ncol = 0)
for(vari in colnames(datos_trans)[8:19]) {
  model<-pgls(as.formula(paste(vari, "~ Loc.1 + kipps + TDIV", sep = "")), data=comp.data, lambda="ML")
  nam <- paste("var_", vari, sep = "")
  assign(nam, model)
  a <- summary(model)
  b <- as.data.frame(a$coefficients[2,])
  b <- rbind(b, a$param[2])
  colnames(b) <- vari
  row.names(b)[5] <- "lambda"
  table[[vari]] <- a
  table2[[vari]] <- a$adj.r.squared
  
  dd <-  as.data.frame(lmp(model))
  colnames(dd) <- vari
  table3_all <- cbind(table3_all, dd)
}

  write.csv(t(table3_all), "1_main_results_PGLS_pvalues.txt", row.names = T)
  
  
table_all <- matrix(nrow = 0, ncol = 3)
n <- names(table[1])
for (n in names(table)) {
aa <- table[[n]]
aa <- as.data.frame(aa$coefficients[2:4,3:4])
aa$va <- n
table_all <- rbind(table_all, aa)
}
write.matrix(table_all, "1_main_results_PGLS.txt")
  
table2_all <- matrix(nrow = 1, ncol = 0)
n <- names(table2[1])
for (n in names(table)) {
  aa <- as.data.frame(table2[[n]])
  colnames(aa) <- n
  table2_all <- cbind(table2_all, aa)
}
write.csv(t(table2_all), "1_main_results_PGLS_r2adj.txt")  



  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ####################################
  ####################################
  ###############################
  ###Testing PGLS using the 42 populations controling for mountain block FULL
  
  comp.data<-comparative.data(arbol, datos_trans, names.col="taxon", vcv.dim=2, warn.dropped=F)
  full_models_42 <- vector(mode = "list", length = ncol(datos_trans[8:15]))
  names(full_models_42) <- colnames(datos_trans)[8:15]
  
  aic_models <- vector(mode = "list", length = ncol(datos_trans[8:15]))
  names(aic_models) <- colnames(datos_trans)[8:15]
  
  aic_wth_all <- data.frame(matrix(ncol = 0, nrow = 4))
  
  table <-data.frame(matrix(ncol = 0, nrow = 5))
  for(vari in colnames(datos_trans)[8:15]) {
    #model1<-pgls(as.formula(paste(vari, "~log.Total_niche_volume.*Geographic.distance*Residuals_resistance_geo*Loc.1", sep = "")), data=comp.data, lambda="ML")
    model2<-pgls(as.formula(paste(vari, "~Geographic.distance+Residuals_resistance_geo+Loc.1", sep = "")), data=comp.data, lambda="ML")
    model3<-pgls(as.formula(paste(vari, "~log.Total_niche_volume.+Residuals_resistance_geo+Loc.1", sep = "")), data=comp.data, lambda="ML")
    model4<-pgls(as.formula(paste(vari, "~log.Total_niche_volume.+Geographic.distance+Loc.1", sep = "")), data=comp.data, lambda="ML")
    model5<-pgls(as.formula(paste(vari, "~log.Total_niche_volume.+Geographic.distance+Residuals_resistance_geo", sep = "")), data=comp.data, lambda="ML")
    model6<-pgls(as.formula(paste(vari, "~Residuals_resistance_geo+Loc.1", sep = "")), data=comp.data, lambda="ML")
    # model<-pgls(as.formula(paste(vari, "~Loc.1", sep = "")), data=comp.data, lambda="ML")
   aic_wth <- aicw(AIC(model2, model3, model4, model5)[,2])
   colnames(aic_wth)[3] <- vari
   aic_wth_all <- cbind(aic_wth_all, aic_wth[3])
    nam <- paste("var_", vari, sep = "")
    assign(nam, model)
    
    a <- summary(model6)
    b <- as.data.frame(a$coefficients)
    
   full_models_42[[vari]] <- a
   #aic_models[[vari]] <- c
  }
  
  full_models_42
  #aic_models
  aic_wth_all
  row.names(aic_wth_all) <- c("- Niche", "- Geo dist", "- resistance", "- locality")
  write.csv(aic_wth_all, "1_aic_wht.csv")
  aic_wth_all
  
  aic_wth_all2 <- as.data.frame(read.csv("1_aic_wht.csv", header = T))
  
  aic_plot <-  qplot(x = var, y = val, fill = model, data = aic_wth_all2, geom = "col") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  plot(aic_plot)
  ##########################################
  ##########################################
  ##########################################
  ##########################################
  ###Testing PGLS using the 42 populations controling for mountain block FULL environment as response
  
  comp.data<-comparative.data(arbol, datos, names.col="taxon", vcv.dim=2, warn.dropped=TRUE)
  full_models_env_42 <- vector(mode = "list", length = ncol(datos[17:20]))
  names(full_models_env_42) <- colnames(datos)[17:20]
  
  table <-data.frame(matrix(ncol = 0, nrow = 5))
  for(vari in colnames(datos)[17:20]) {
    model<-pgls(as.formula(paste(vari, "~Loc.1", sep = "")), data=comp.data, lambda="ML")
    # model<-pgls(as.formula(paste(vari, "~Loc.1", sep = "")), data=comp.data, lambda="ML")
    
    
    nam <- paste("var_", vari, sep = "")
    assign(nam, model)
    a <- summary(model)
    b <- as.data.frame(a$coefficients)
    
    full_models_env_42[[vari]] <- b
    
  }
  full_models_env_42
  
  
  
  
  
  ####################################
  ####################################
  ###############################
  ###Testing PGLS using the diferences between north and south for the variables (table with 21 species)
  
  datos21<-read.csv("1_PGLS_table_21.csv",header=TRUE)
  arbol21<-read.tree("1_21sp_tree.tre")
  
  
    
  comp.data<-comparative.data(arbol21, datos21, names.col="taxon", vcv.dim=2, warn.dropped=TRUE)
  full_models_21 <- vector(mode = "list", length = ncol(datos21[7:15]))
  names(full_models_21) <- colnames(datos21)[7:15]
  
  for(vari in colnames(datos21)[7:15]) {
    model<-pgls(as.formula(paste(vari, "~niche.temp+Niche.prec+Geographic.distance+Resistance.normalized", sep = "")), data=comp.data, lambda="ML")
    # model<-pgls(as.formula(paste(vari, "~Loc.1", sep = "")), data=comp.data, lambda="ML")
    
    
    nam <- paste("var_", vari, sep = "")
    assign(nam, model)
    a <- summary(model)
    b <- as.data.frame(a$coefficients)
    
    full_models_21[[vari]] <- b
    
  }
  full_models_21
  
  #model<-pgls(as.formula(paste(vari, "~niche.temp*Niche.prec*Geographic.distance*Resistance.normalized", sep = "")), data=comp.data, lambda="ML")
 ## }
  #model<-pgls(as.formula(paste(vari, "~Loc.1", sep = "")), data=comp.data, lambda="ML")
  
  
  
  
  ##################################################################################
  ###############################################################################
  ########Testing PGLS for the North
  
  datos_N_21<-datos[1:21,]
  datos_N_21$taxon <- datos21$taxon
  arbol21<-read.tree("1_21sp_tree.tre")
  
  comp.data<-comparative.data(arbol21, datos_N_21, names.col="taxon", vcv.dim=2, warn.dropped=TRUE)
  full_models_N_21 <- vector(mode = "list", length = ncol(datos_N_21[8:16]))
  names(full_models_N_21) <- colnames(datos_N_21)[8:16]
  
  for(vari in colnames(datos_N_21)[8:16]) {
    model<-pgls(as.formula(paste(vari, "~niche.temp+Niche.prec+Geographic.distance+Resistance.normalized", sep = "")), data=comp.data, lambda="ML")
    # model<-pgls(as.formula(paste(vari, "~Loc.1", sep = "")), data=comp.data, lambda="ML")
    
    
    nam <- paste("var_", vari, sep = "")
    assign(nam, model)
    a <- summary(model)
    b <- as.data.frame(a$coefficients)
    
    full_models_N_21[[vari]] <- b
    
  }
  full_models_N_21
  
  
  
  ##################################################################################
  ###############################################################################
  ########Testing PGLS for the South
  
  datos_S_21<-datos[22:42,]
  datos_S_21$taxon <- datos21$taxon
  arbol21<-read.tree("1_21sp_tree.tre")
  
  comp.data<-comparative.data(arbol21, datos_S_21, names.col="taxon", vcv.dim=2, warn.dropped=TRUE)
  full_models_S_21 <- vector(mode = "list", length = ncol(datos_S_21[8:16]))
  names(full_models_S_21) <- colnames(datos_S_21)[8:16]
  
  for(vari in colnames(datos_S_21)[8:16]) {
    model<-pgls(as.formula(paste(vari, "~niche.temp+Niche.prec+Geographic.distance+Resistance.normalized", sep = "")), data=comp.data, lambda="ML")
    # model<-pgls(as.formula(paste(vari, "~Loc.1", sep = "")), data=comp.data, lambda="ML")
    
    
    nam <- paste("var_", vari, sep = "")
    assign(nam, model)
    a <- summary(model)
    b <- as.data.frame(a$coefficients)
    
    full_models_S_21[[vari]] <- b
    
  }
  full_models_S_21
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
 