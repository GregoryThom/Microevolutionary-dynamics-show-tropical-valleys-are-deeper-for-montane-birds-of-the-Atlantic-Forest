#####################LMEM###############
library(lme4)
library(MuMIn)
library(ggplot2)
library(Rcpp)
library(corMLPE)
library(nlme)
library(geiger)
library(caper)
library(Hmisc)
library(ggplot2)
library(dplyr)
library(sjstats)
library(mgcv)
library(ggpubr)
#source("Data_Script/script/corMLPE.R")
orderPops <- function(pops){
  popn <- matrix(as.numeric(factor(as.vector(pops))), nrow=nrow(pops), ncol=ncol(pops))
  popn <- t(apply(popn, 1, function(x) x[order(x)]))
  attr(popn, "order") <- order(popn[,1], popn[,2])
  popn
}

reorderData <- function(data, pop.name, group.name){
  require(plyr)
  ddply(data, group.name, function(z){
    mm <- cbind(as.character(z[,pop.name[[1]]]), as.character(z[,pop.name[[2]]]) )
    np <- orderPops(mm)
    nn <- paste0(pop.name, "_num")
    colnames(np) <- nn
    z <- data.frame(z, np)
    z[attr(np,"order"), ]
  })
}


setwd("~/Dropbox/1_pos_DOC/desenvolvimento/1_paper_22sp/1_comparative_analyses/")

##########ALL_individuals##########
#DMA <- read.csv("1_dist_all_sp_N_and_S_and_NS_normalized.csv", header=T)
#DMA <- read.csv("1_dist_all_sp_N_and_S_and_NS_no_outlier.csv", header=T)
DMA <- read.csv("1_dist_all_sp_N_and_S_and_NS.csv", header=T) ###banco de dados com mesoleuca removida
## Transformations
#logit<-function(x) {x=x+0.01; log(x/(1-x))}
#DMA$dist <- DMA$dist^2
#DMA$Dgeo_m.dist <- sqrt(DMA$Dgeo_m.dist)
#DMA$Dgeo_m.dist <- scale(DMA$Dgeo_m.dist, scale=T)

DMA$resist_dist.dist <- sqrt(DMA$resist_dist.dist)
DMA$resist_dist.dist <- scale(DMA$resist_dist.dist, scale=T)

DMA$saz_dist.dist <- sqrt(DMA$saz_dist.dist)
DMA$saz_dist.dist <- scale(DMA$saz_dist.dist, scale=T)

DMA$sta_dist.dist <- sqrt(DMA$sta_dist.dist)
DMA$sta_dist.dist <- scale(DMA$sta_dist.dist, scale=T)

DMA$re_sa_dist.dist <- sqrt(DMA$re_sa_dist.dist)
DMA$re_sa_dist.dist <- scale(DMA$re_sa_dist.dist, scale=T)

DMA$re_st_dist.dist <- sqrt(DMA$re_st_dist.dist)
DMA$re_st_dist.dist <- scale(DMA$re_st_dist.dist, scale=T)

DMA$sa_st_dist.dist <- sqrt(DMA$sa_st_dist.dist)
DMA$sa_st_dist.dist <- scale(DMA$sa_st_dist.dist, scale=T)

DMA$sa_re_st_dist.dist <- sqrt(DMA$sa_re_st_dist.dist)
DMA$sa_re_st_dist.dist <- scale(DMA$sa_re_st_dist.dist, scale=T)

DMA$resist_resi <- scale(DMA$resist_resi, scale=T)

DMA$saz_resi <- scale(DMA$saz_resi, scale=T)

DMA$sta_resi <- scale(DMA$sta_resi, scale=T)

DMA$TDIV <- log(DMA$TDIV)
DMA$TDIV <- scale(DMA$TDIV, scale=T)

DMA$Kipps <- log(DMA$Kipps)
DMA$Kipps <- scale(DMA$Kipps, scale=T)


#ggqqplot(DMA$re_sa_dist.dist)
#ggdensity(DMA$re_sa_dist.dist, 
#          main = "Density plot of tooth length",
#          xlab = "Tooth length")
#shapiro.test(DMA$re_sa_dist.dist)




RDMA <- reorderData(DMA, c("iso1", "iso2"), "sp")

cont <- lmeControl(maxIter=1000, msMaxIter=3000, msMaxEval=5000, msVerbose=FALSE, opt="optim")








newBM_geo <- lme(fixed=dist ~ 1 + Dgeo_m.dist + TDIV + Kipps + env + Kipps + env,
             random=~Dgeo_m.dist|sp/family,
             correlation=corMLPE(form=~iso1_num+iso2_num|sp/family),
             data=RDMA, method = "ML", control=cont)
#############teste
#newBM_geo <- mgcv::gamm(dist ~ 1 + s(Dgeo_m.dist) + s(resist_dist.dist) + ti(resist_dist.dist,Dgeo_m.dist) + s(saz_dist.dist) + s(saz_dist.dist,Dgeo_m.dist) + s(TDIV),
#                 random= list(family=~1, sp=~1),
#                 correlation=corMLPE(form=~iso1_num+iso2_num|sp/family),
#                 data=RDMA, control=cont, select=T)
#############


newBM_resis <- lme(fixed=dist ~ 1 + resist_dist.dist + TDIV + Kipps + env,
             random=~Dgeo_m.dist|sp/family,
             correlation=corMLPE(form=~iso1_num+iso2_num|sp/family),
             data=RDMA, method = "ML", control=cont)

newBM_saz <- lme(fixed=dist ~ 1 + saz_dist.dist + TDIV + Kipps + env,
              random=~Dgeo_m.dist|sp/family,
              correlation=corMLPE(form=~iso1_num+iso2_num|sp/family),
              data=RDMA, method = "ML", control=cont)

newBM_sta <- lme(fixed=dist ~ 1 + sta_dist.dist + TDIV + Kipps + env,
                 random=~Dgeo_m.dist|sp/family,
                 correlation=corMLPE(form=~iso1_num+iso2_num|sp/family),
                 data=RDMA, method = "ML", control=cont)
#########

#newBM_re_sa <- lme(fixed=dist ~ 1 + re_sa_dist.dist + TDIV + Kipps + env,
#                 random=~Dgeo_m.dist|sp/family,
#                 correlation=corMLPE(form=~iso1_num+iso2_num|sp/family),
#                 data=RDMA, method = "ML", control=cont)

#newBM_re_st <- lme(fixed=dist ~ 1 + re_st_dist.dist + TDIV + Kipps + env,
#                 random=~Dgeo_m.dist|sp/family,
#                 correlation=corMLPE(form=~iso1_num+iso2_num|sp/family),
#                 data=RDMA, method = "ML", control=cont)

#newBM_sa_st <- lme(fixed=dist ~ 1 + sa_st_dist.dist + TDIV + Kipps + env,
#                 random=~Dgeo_m.dist|sp/family,
#                 correlation=corMLPE(form=~iso1_num+iso2_num|sp/family),
#                 data=RDMA, method = "ML", control=cont)

#newBM_sa_re_st <- lme(fixed=dist ~ 1 + sa_re_st_dist.dist + TDIV + Kipps + env,
#                 random=~Dgeo_m.dist|sp/family,
#                 correlation=corMLPE(form=~iso1_num+iso2_num|sp/family),
#                 data=RDMA, method = "ML", control=cont)

#######################




summary(newBM_geo)
summary(newBM_resis)
summary(newBM_saz)
summary(newBM_sta)
#summary(newBM_re_sa)
#summary(newBM_re_st)
#summary(newBM_sa_st)
#summary(newBM_sa_re_st)

r.squaredGLMM(newBM_geo)
r.squaredGLMM(newBM_resis)
r.squaredGLMM(newBM_saz)
r.squaredGLMM(newBM_sta)
#r.squaredGLMM(newBM_re_sa)
#r.squaredGLMM(newBM_re_st)
#r.squaredGLMM(newBM_sa_st)
#r.squaredGLMM(newBM_sa_re_st)

aicw(AIC(newBM_geo, newBM_resis, newBM_saz, newBM_sta)[,2])
anova(newBM_geo, newBM_resis, newBM_saz, newBM_sta)

#########################

newBMfull <- lme(fixed=dist ~ 1 + Dgeo_m.dist + resist_resi + saz_resi + sta_resi +
                   TDIV + Kipps + env,
              random=~Dgeo_m.dist|sp/family,
              correlation=corMLPE(form=~iso1_num+iso2_num|sp/family),
              data=RDMA, method = "ML", control=cont)

aicw(AIC(newBM_geo, newBM_resis, newBM_saz, newBM_sta,newBMfull)[,2])
anova(newBM_geo, newBM_resis, newBM_saz, newBM_sta,newBMfull)
summary(newBMfull)
r.squaredGLMM(newBMfull)
#######################################################
######################NORTH############################
######################NORTH############################
######################NORTH############################
######################NORTH############################
#######################################################

##########NORTH_individuals##########
DMA <- read.csv("1_dist_all_sp_N_and_S_and_NS.csv", header=T)
DMA <- DMA[grep("north", DMA$region),]
#DMA <- read.csv("1_dist_N.csv", header=T)
## Transformations
#logit<-function(x) {x=x+0.01; log(x/(1-x))}
#DMA$dist <- DMA$dist^2
DMA$Dgeo_m.dist <- sqrt(DMA$Dgeo_m.dist)
DMA$Dgeo_m.dist <- scale(DMA$Dgeo_m.dist, scale=T)

DMA$resist_dist.dist <- sqrt(DMA$resist_dist.dist)
DMA$resist_dist.dist <- scale(DMA$resist_dist.dist, scale=T)

DMA$saz_dist.dist <- sqrt(DMA$saz_dist.dist)
DMA$saz_dist.dist <- scale(DMA$saz_dist.dist, scale=T)

DMA$sta_dist.dist <- sqrt(DMA$sta_dist.dist)
DMA$sta_dist.dist <- scale(DMA$sta_dist.dist, scale=T)

DMA$re_sa_dist.dist <- sqrt(DMA$re_sa_dist.dist)
DMA$re_sa_dist.dist <- scale(DMA$re_sa_dist.dist, scale=T)

DMA$re_st_dist.dist <- sqrt(DMA$re_st_dist.dist)
DMA$re_st_dist.dist <- scale(DMA$re_st_dist.dist, scale=T)

DMA$sa_st_dist.dist <- sqrt(DMA$sa_st_dist.dist)
DMA$sa_st_dist.dist <- scale(DMA$sa_st_dist.dist, scale=T)

DMA$sa_re_st_dist.dist <- sqrt(DMA$sa_re_st_dist.dist)
DMA$sa_re_st_dist.dist <- scale(DMA$sa_re_st_dist.dist, scale=T)

DMA$resist_resi <- scale(DMA$resist_resi, scale=T)

DMA$saz_resi <- scale(DMA$saz_resi, scale=T)

DMA$sta_resi <- scale(DMA$sta_resi, scale=T)

DMA$TDIV <- log(DMA$TDIV)
DMA$TDIV <- scale(DMA$TDIV, scale=T)

DMA$Kipps <- log(DMA$Kipps)
DMA$Kipps <- scale(DMA$Kipps, scale=T)


RDMA <- reorderData(DMA, c("iso1", "iso2"), "sp")

cont <- lmeControl(maxIter=1000, msMaxIter=3000, msMaxEval=5000, msVerbose=FALSE, opt="optim")

north_geo <- lme(fixed=dist ~ 1 + Dgeo_m.dist + TDIV + Kipps + env,
                 random=~Dgeo_m.dist|sp/family,
                 correlation=corMLPE(form=~iso1_num+iso2_num|sp/family),
                 data=RDMA, method = "ML", control=cont)

north_resis <- lme(fixed=dist ~ 1 + resist_dist.dist + TDIV + Kipps + env,
                   random=~Dgeo_m.dist|sp/family,
                   correlation=corMLPE(form=~iso1_num+iso2_num|sp/family),
                   data=RDMA, method = "ML", control=cont)

north_saz <- lme(fixed=dist ~ 1 + saz_dist.dist + TDIV + Kipps + env,
                 random=~Dgeo_m.dist|sp/family,
                 correlation=corMLPE(form=~iso1_num+iso2_num|sp/family),
                 data=RDMA, method = "ML", control=cont)

north_sta <- lme(fixed=dist ~ 1 + sta_dist.dist + TDIV + Kipps + env,
                 random=~Dgeo_m.dist|sp/family,
                 correlation=corMLPE(form=~iso1_num+iso2_num|sp/family),
                 data=RDMA, method = "ML", control=cont)
#########

#north_re_sa <- lme(fixed=dist ~ 1 + re_sa_dist.dist + TDIV + Kipps + env,
#                   random=~Dgeo_m.dist|sp/family,
#                   correlation=corMLPE(form=~iso1_num+iso2_num|sp/family),
#                   data=RDMA, method = "ML", control=cont)

#north_re_st <- lme(fixed=dist ~ 1 + re_st_dist.dist + TDIV + Kipps + env,
#                   random=~Dgeo_m.dist|sp/family,
#                   correlation=corMLPE(form=~iso1_num+iso2_num|sp/family),
#                   data=RDMA, method = "ML", control=cont)

#north_sa_st <- lme(fixed=dist ~ 1 + sa_st_dist.dist + TDIV + Kipps + env,
#                   random=~Dgeo_m.dist|sp/family,
#                   correlation=corMLPE(form=~iso1_num+iso2_num|sp/family),
#                   data=RDMA, method = "ML", control=cont)

#north_sa_re_st <- lme(fixed=dist ~ 1 + sa_re_st_dist.dist + TDIV + Kipps + env,
#                      random=~Dgeo_m.dist|sp/family,
#                      correlation=corMLPE(form=~iso1_num+iso2_num|sp/family),
#                      data=RDMA, method = "ML", control=cont)

#######################

summary(north_geo)
summary(north_resis)
summary(north_saz)
summary(north_sta)
#summary(north_re_sa)
#summary(north_re_st)
#summary(north_sa_st)
#summary(north_sa_re_st)

r.squaredGLMM(north_geo)
r.squaredGLMM(north_resis)
r.squaredGLMM(north_saz)
r.squaredGLMM(north_sta)
#r.squaredGLMM(north_re_sa)
#r.squaredGLMM(north_re_st)
#r.squaredGLMM(north_sa_st)
#r.squaredGLMM(north_sa_re_st)


aicw(AIC(north_geo, north_resis, north_saz, north_sta)[,2])
anova(north_geo, north_resis, north_saz, north_sta)
summary(north_resis)

#########################

northfull <- lme(fixed=dist ~ 1 + Dgeo_m.dist + resist_resi + saz_resi 
                 + sta_resi + TDIV + Kipps + env,
                 random=~Dgeo_m.dist|sp/family,
                 correlation=corMLPE(form=~iso1_num+iso2_num|sp/family),
                 data=RDMA, method = "ML", control=cont)


aicw(AIC(north_geo, north_resis, north_saz, north_sta,northfull)[,2])
anova(north_geo, north_resis, north_saz, north_sta,northfull)
summary(northfull)
r.squaredGLMM(northfull)

##Special_model North wo/ stability
#nort_special <- lme(fixed=dist ~ 1 + Dgeo_m.dist + resist_dist.dist + saz_dist.dist 
#                    + TDIV + Kipps + env,
#                 random=~Dgeo_m.dist|sp/family,
#                 correlation=corMLPE(form=~iso1_num+iso2_num|sp/family),
#                 data=RDMA, method = "ML", control=cont)

#aicw(AIC(north_geo, north_resis, north_saz, north_sta,northfull, nort_special)[,2])
#anova(north_geo, north_resis, north_saz, north_sta,northfull, nort_special)
#summary(nort_special)
#r.squaredGLMM(nort_special)

#######################################################
######################SOUTH############################
######################SOUTH############################
######################SOUTH############################
######################SOUTH############################
#######################################################


##########SOUTH_individuals##########
DMA <- read.csv("1_dist_all_sp_N_and_S_and_NS.csv", header=T)
DMA <- DMA[grep("south", DMA$region),]
#DMA <- read.csv("1_dist_N.csv", header=T)
## Transformations
#logit<-function(x) {x=x+0.01; log(x/(1-x))}
#DMA$dist <- DMA$dist^2
DMA$Dgeo_m.dist <- sqrt(DMA$Dgeo_m.dist)
DMA$Dgeo_m.dist <- scale(DMA$Dgeo_m.dist, scale=T)

DMA$resist_dist.dist <- sqrt(DMA$resist_dist.dist)
DMA$resist_dist.dist <- scale(DMA$resist_dist.dist, scale=T)

DMA$saz_dist.dist <- sqrt(DMA$saz_dist.dist)
DMA$saz_dist.dist <- scale(DMA$saz_dist.dist, scale=T)

DMA$sta_dist.dist <- sqrt(DMA$sta_dist.dist)
DMA$sta_dist.dist <- scale(DMA$sta_dist.dist, scale=T)

DMA$re_sa_dist.dist <- sqrt(DMA$re_sa_dist.dist)
DMA$re_sa_dist.dist <- scale(DMA$re_sa_dist.dist, scale=T)

DMA$re_st_dist.dist <- sqrt(DMA$re_st_dist.dist)
DMA$re_st_dist.dist <- scale(DMA$re_st_dist.dist, scale=T)

DMA$sa_st_dist.dist <- sqrt(DMA$sa_st_dist.dist)
DMA$sa_st_dist.dist <- scale(DMA$sa_st_dist.dist, scale=T)

DMA$sa_re_st_dist.dist <- sqrt(DMA$sa_re_st_dist.dist)
DMA$sa_re_st_dist.dist <- scale(DMA$sa_re_st_dist.dist, scale=T)

DMA$resist_resi <- scale(DMA$resist_resi, scale=T)

DMA$saz_resi <- scale(DMA$saz_resi, scale=T)

DMA$sta_resi <- scale(DMA$sta_resi, scale=T)

#ggqqplot(DMA$sta_resi)
#ggdensity(DMA$sta_resi, 
#          main = "Density plot of tooth length",
#          xlab = "Tooth length")
#shapiro.test(DMA$sta_resi)


DMA$TDIV <- log(DMA$TDIV)
DMA$TDIV <- scale(DMA$TDIV, scale=T)

DMA$Kipps <- log(DMA$Kipps)
DMA$Kipps <- scale(DMA$Kipps, scale=T)

RDMA <- reorderData(DMA, c("iso1", "iso2"), "sp")

cont <- lmeControl(maxIter=1000, msMaxIter=3000, msMaxEval=5000, msVerbose=FALSE, opt="optim")

south_geo <- lme(fixed=dist ~ 1 + Dgeo_m.dist + TDIV + Kipps + env,
                 random=~Dgeo_m.dist|sp/family,
                 correlation=corMLPE(form=~iso1_num+iso2_num|sp/family),
                 data=RDMA, method = "ML", control=cont)

south_resis <- lme(fixed=dist ~ 1 + resist_dist.dist + TDIV + Kipps + env,
                   random=~Dgeo_m.dist|sp/family,
                   correlation=corMLPE(form=~iso1_num+iso2_num|sp/family),
                   data=RDMA, method = "ML", control=cont)

south_saz <- lme(fixed=dist ~ 1 + saz_dist.dist + TDIV + Kipps + env,
                 random=~Dgeo_m.dist|sp/family,
                 correlation=corMLPE(form=~iso1_num+iso2_num|sp/family),
                 data=RDMA, method = "ML", control=cont)

south_sta <- lme(fixed=dist ~ 1 + sta_dist.dist + TDIV + Kipps + env,
                 random=~Dgeo_m.dist|sp/family,
                 correlation=corMLPE(form=~iso1_num+iso2_num|sp/family),
                 data=RDMA, method = "ML", control=cont)
#########

#south_re_sa <- lme(fixed=dist ~ 1 + re_sa_dist.dist + TDIV + Kipps + env,
#                   random=~Dgeo_m.dist|sp/family,
#                   correlation=corMLPE(form=~iso1_num+iso2_num|sp/family),
#                   data=RDMA, method = "ML", control=cont)

#south_re_st <- lme(fixed=dist ~ 1 + re_st_dist.dist + TDIV + Kipps + env,
#                   random=~Dgeo_m.dist|sp/family,
#                   correlation=corMLPE(form=~iso1_num+iso2_num|sp/family),
#                   data=RDMA, method = "ML", control=cont)

#south_sa_st <- lme(fixed=dist ~ 1 + sa_st_dist.dist + TDIV + Kipps + env,
#                   random=~Dgeo_m.dist|sp/family,
#                   correlation=corMLPE(form=~iso1_num+iso2_num|sp/family),
#                   data=RDMA, method = "ML", control=cont)

#south_sa_re_st <- lme(fixed=dist ~ 1 + sa_re_st_dist.dist + TDIV + Kipps + env,
#                      random=~Dgeo_m.dist|sp/family,
#                      correlation=corMLPE(form=~iso1_num+iso2_num|sp/family),
#                      data=RDMA, method = "ML", control=cont)

#######################

summary(south_geo)
summary(south_resis)
summary(south_saz)
summary(south_sta)
#summary(south_re_sa)
#summary(south_re_st)
#summary(south_sa_st)
#summary(south_sa_re_st)

r.squaredGLMM(south_geo)
r.squaredGLMM(south_resis)
r.squaredGLMM(south_saz)
r.squaredGLMM(south_sta)
#r.squaredGLMM(south_re_sa)
#r.squaredGLMM(south_re_st)
#r.squaredGLMM(south_sa_st)
#r.squaredGLMM(south_sa_re_st)

aicw(AIC(south_geo, south_resis, south_saz, south_sta)[,2])
anova(south_geo, south_resis, south_saz, south_sta)
summary(south_geo)

#########################

southfull <- lme(fixed=dist ~ 1 + Dgeo_m.dist + resist_resi + saz_resi + sta_resi +
                   TDIV + Kipps + env,
                 random=~Dgeo_m.dist|sp/family,
                 correlation=corMLPE(form=~iso1_num+iso2_num|sp/family),
                 data=RDMA, method = "ML", control=cont)

#south <- lme(fixed=dist ~ 1 + Dgeo_m.dist + resist_dist.dist + saz_dist.dist,

summary(southfull)
r.squaredGLMM(southfull)
aicw(AIC(south_geo, south_resis, south_saz, south_sta,southfull)[,2])
anova(south_geo, south_resis, south_saz, south_sta,southfull)


###########FINAL RESULTS###################################!
############################################################
a <- aicw(AIC(newBM_geo, newBM_resis, newBM_saz, newBM_sta,newBMfull)[,2])
b <- anova(newBM_geo, newBM_resis, newBM_saz, newBM_sta,newBMfull)
c <- summary(newBMfull)
d <- r.squaredGLMM(newBMfull)

write.csv(a, "final_results/1_all_FA_aicw.csv")
write.csv(c$tTable, "final_results/1_all_FA_summary_full_model.csv")
write.csv(d, "final_results/1_all_FA_r.squared_full_model.csv")

a <- aicw(AIC(north_geo, north_resis, north_saz, north_sta,northfull)[,2])
b <- anova(north_geo, north_resis, north_saz, north_sta,northfull)
c <- summary(northfull)
d <- r.squaredGLMM(northfull)

write.csv(a, "final_results/1_north_aicw.csv")
write.csv(c$tTable, "final_results/1_north_summary_full_model.csv")
write.csv(d, "final_results/1_north_r.squared_full_model.csv")

a <- aicw(AIC(south_geo, south_resis, south_saz, south_sta, southfull)[,2])
b <- anova(south_geo, south_resis, south_saz,south_sta, southfull)
c <- summary(south_geo)
d <- r.squaredGLMM(south_geo)

write.csv(a, "final_results/1_south_aicw.csv")
write.csv(c$tTable, "final_results/1_south_summary_geo_model.csv")
write.csv(d, "final_results/1_south_r.squared_geo_model.csv")


aicw(AIC(newBM_geo, newBM_resis, newBM_saz, newBM_sta,newBMfull)[,2])
aicw(AIC(north_geo, north_resis, north_saz, north_sta,northfull)[,2])
aicw(AIC(south_geo, south_resis, south_saz, south_sta, southfull)[,2])




newBM <- lme(fixed= dist ~ 1,
             random=~1|sp,
             correlation=corMLPE(form=~iso1_num+iso2_num|sp),
             data=RDMA, method = "ML", control=cont)
summary(newBM)             
r.squaredGLMM(newBM)             
x <- as.data.frame(newBM$residuals)  
colnames(x) <- c("fixed", "random")


RDMA2 <- cbind(RDMA, x)


fit1 <- lm(random ~ Dgeo_m.dist, data = RDMA2)
summ1 <- summary(fit1)
R2geo_ind <- as.data.frame(cbind(summ1$adj.r.squared, summ1$coef[2,4]))
#row.names(R2geo_ind) <- file
#R2geo <- rbind(R2geo, R2geo_ind)



pdf("../01a_manuscrito/1_Figures/Regression_geo_gen.pdf")
ggplot(RDMA2, aes(x = Dgeo_m.dist, y = random)) + 
  stat_density_2d(geom = "polygon", aes(alpha = (..level..)^2, fill = region), show.legend = FALSE) +
  geom_point(show.legend = FALSE, size=.7, aes(col=region))+
  stat_smooth(method = "lm", col = "black") + theme_classic() + 
  scale_fill_manual(values = c("#00AFBB", "#00CC00", "#FC4E07")) + 
  scale_color_manual(values = c("#3333FF", "#339900", "#FC4E07")) +
  theme(legend.title = element_blank())
dev.off()
