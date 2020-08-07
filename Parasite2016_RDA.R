library(ggvegan)
library(ggthemes)
library(raster)
library(ggplot2)
library(readr)
library(gridExtra)
library(plyr)
library(afex)
library(MASS)
library(dplyr) # building data matrix
library(glmmTMB) # GLMMs
library(car) # ANOVA
library(vegan) # for MEMs

# Set working directory
setwd('C:/Users/pascalh/Documents/GitHub/Stickleback-parasites-2016')

# Load data
data_2016 <- read.csv("data_2016.csv", sep=';') #field and parasite data
env <- read.csv("Environment_R.csv", sep=',') #all environmental data
#env_av <- read.csv("Env_av.csv", sep=';') #environmental variables (average values)
#env_max <- read.csv("env_max.csv", sep=';') #environmental variables (max. values)
spavar <- read.csv("space2.csv", sep=';') #spatial variables: network centrality and upstream distance
distance_matrix <- read.csv("distance_matrix.csv", sep=';') #spatial variables: distance matrix

# Calculation of MEMs  ####
KautoDist <- data.frame(cmdscale(distance_matrix),rownames(distance_matrix));colnames(distance_matrix)<- c("X1","X2","LocationID")
spa <- KautoDist[,c(1,2)]
library(raster)
dist <- pointDistance(spa, allpairs = TRUE, lonlat = FALSE)
spa.dist <- as.dist(dist)

# Test and forward selectio of PCNM variables
spanning <- spantree(spa.dist)
dmin <- max(spanning$dist)
spa.dist[spa.dist>dmin] <- 4*dmin
xy.PCoA <- cmdscale(spa.dist,k=nrow(spa)-1,eig=TRUE)
nb.ev<-length(which(xy.PCoA$eig > 0.0001))
xy.PCoA$eig

# Broken stick model
ev <- xy.PCoA$eig
ev <- subset(ev, xy.PCoA$eig > 0.0001)
n <- length(ev)
bsm <- data.frame(j=seq(1:n), p=0)
bsm$p[1] <- 1/n
for (i in 2:n) {
  bsm$p[i] = bsm$p[i-1] + (1/(n+1-i))
}
bsm$p <- 100*bsm$p/n
bsm
par(mfrow=c(1,1))
barplot(ev, main="Eigenvalues", col="bisque", las=2)
abline(h=mean(ev), col="red")
legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")
barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=TRUE, main="% variance", col=c("bisque", 2), las=2)
nb.ev <- length(which(xy.PCoA$eig > mean(ev)))

spa.PCNM <- as.data.frame(xy.PCoA$points[1:nrow(spa),1:3])

cor(spavar$updist, spa.PCNM$V1)
plot(spavar$updist, spa.PCNM$V1) #first MEM corresponds to distance from Demer-Dijle confluence
cor(spavar$netcen, spa.PCNM$V2)
plot(spavar$netcen, spa.PCNM$V2) #second MEM corresponds to network centrality

# remove correlated variables
#write.csv(cor(cbind(env[,-1], spavar), use="pairwise.complete.obs"), file="collinearity.csv")
# uncorrelated variables: T_max, T_av, T_min, con_min, COD_min, KjN_min, NH4_min, NO3_min, SM_min, T, con, COD, NO3, SM
# correlation within variable type: pH (all 4 measures, keep pH_av), O2_av and O2_min (keep O2_av), O2_sat_av and O2_sat_min (keep 02_sat_av), Cl (all except Cl, which is correlated with Cl_max -> keep only Cl_av); COD_av and COD_max (keep COD_av); 
# correlation between variable types: pH_av and con_av (keep con_av), O2_av and O2_sat_av (keep O2_sat_av), O2_max and O2_sat_max (keep O2_sat_max), O2 and O2_sat (keep O2_sat)
# BOC, NO2 oPO4 -> keep only NO2_av
# NO3 and Nt -> keep only NO3_av
# KjN, NH4 and Pt -> keep only NH4_av
# pH and conductivity -> keep only con_av
# SO4: all measures correlated -> keep only SO4_av
# SM: all measures correlated -> keep only SM_av
# Cl: all measures correlated -> keep only Cl_av
# COD: all measures somewhat correlated -> keep only COD_av
# temp: all measures somewhat correlated -> keep only T_av
# O2 and O2_sat: all measures somewhat correlated -> keep only O2_sat_av
# speciesrichness is marginally related to network centrality (corr.coef. -0.60) -> keep only network centrality
# updist and updist2 are correlated -> keep only updist
# updist3 is not related to any other parameter -> keep

#write.csv(cor(cbind(env[,c("T_av","con_av","O2_sat_av","Cl_av","COD_av","NH4_av","NO3_av","NO2_av","SO4_av")], spavar), use="pairwise.complete.obs"), file="collinearity_selected.csv")

# make a new data frame by combinding parasite, environmental and space data
env_exp <- env %>% slice(rep(1:n(), table(as.factor(data_2016$site))))
#env_av_exp <- env_av %>% slice(rep(1:n(), table(as.factor(data_2016$site))))
#env_max_exp <- env_max %>% slice(rep(1:n(), table(as.factor(data_2016$site))))
spavar_exp <- spavar %>% slice(rep(1:n(), table(as.factor(data_2016$site))))
spa.PCNM_exp <- spa.PCNM %>% slice(rep(1:n(), table(as.factor(data_2016$site))))

data0 <- cbind(data_2016, env_exp[,-1], spavar_exp, spa.PCNM_exp)
NAs <- 1>rowSums(is.na(data0[,c("weight", "length", "Sex")])) # identify fish with any of the following data missing: length, weight, sex
table(NAs)
data <- data0[NAs,] # remove fish with missing data
data$site <- as.factor(data$site)
data$fish <- as.factor(data$fish)
data$length <- as.numeric(data$length)


#### CALCULATE PARAMETERS ####
names(data)
parsum = aggregate(data[,c(23:25,27:33)], by = list(data[,1]), function(x){sum(x, na.rm = T)}) 
avin = aggregate(data[,c(23:25,27:33)], by = list(data[,1]), function(x){mean(x[x >0], na.rm = T)}); avin[is.na(avin)] <- 0
avab = aggregate(data[,c(23:25,27:33)], by = list(data[,1]), function(x){mean(x, na.rm =T)})
prev = aggregate(data[,c(23:25,27:33)], by = list(data[,1]), function(x){sum(x >0, na.rm = T)/length(x)})
medin = aggregate(data[,c(23:25,27:33)], by = list(data[,1]), function(x){median(x[x >0], na.rm = T)}) ; medin[is.na(medin)] <- 0

#parasite data is overdispersed (mostly so for Trichodina), if using average abundance data, species matrix needs to be transformed
datao <- na.omit(data[,c(1,23:25,27:33)])
ddata <- dispweight(datao[,-1])
avab <- aggregate(ddata, by = list(datao[,1]), function(x){mean(x, na.rm =T)})

# Hellinger transformation of the species dataset ####
# Component communties: Average abundance is used to compensate for differences in sampling size between sites
spe.hel.avab <- vegdist(decostand(cbind(avab[,-1]), "hellinger"), method="bray")
spe.hel.medin <- vegdist(decostand(medin[,-1],"hellinger"), method="euc")

# Infracommunities: Dissimilarities are calculated at the individual level and then averaged within site; a dummy parasite species is added to avoid similarties due to zero-parasites
data_infra <- na.omit(data[,c(1,23:25,27:33)])
data_infra_disp <- dispweight(datao[,-1])
braycurtis <- vegdist(decostand(cbind(data_infra_disp,rep(1,nrow(data_infra))), na.rm=T, method="hellinger"), method="bray", na.rm=T)
meandist <- meandist(braycurtis, data_infra[,1])

env_select <- env[,c("T_av","con_av","O2_sat_av","Cl_av","COD_av","NH4_av","NO3_av","NO2_av")]
#env_select <- env[,c("T","con","O2_sat","Cl","COD","NH4","NO3","NO2")]
#env_select <- env[,c("T_max","con_max","O2_sat_max","Cl_max","COD_max","NH4_max","NO3_max","NO2_max")]

plot(env_select$NO2_av)

spe.rda <- dbrda(meandist ~ T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av, env_select)
#spe.rda <- dbrda(meandist ~ Cl_av, env_select)
spe.rda <- dbrda(spe.hel.avab ~ T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av, env_select)
spe.hel <- decostand(cbind(avab[,-1]), "hellinger")
spe.rda <- rda(spe.hel ~ T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av, env_select)

plot(spe.rda, scaling = 1)
summary(spe.rda)
anova(spe.rda)
anova(spe.rda, by="term")
anova.cca(spe.rda, step=1000);
anova.cca(spe.rda, step=1000, by="term");
RsquareAdj(spe.rda)$adj.r.squared;
RsquareAdj(spe.rda)$r.squared

mod0 <- dbrda(meandist ~ 1, env_select)  # Model with intercept only  #edit_PH
mod1 <- dbrda(meandist ~ ., env_select)  # Model with all explanatory variables  #edit_PH
step.res <- ordiR2step(mod0, mod1, direction = "both",perm.max = 200)
step.res$anova  # Summary table

# Check whether exclusion of site 12 changes results
mod0 <- dbrda(meandist[-c(10),-c(10)] ~ 1, env_select[-c(10),])  # Model with intercept only  #edit_PH
mod1 <- dbrda(meandist[-c(10),-c(10)] ~ ., env_select[-c(10),])  # Model with all explanatory variables  #edit_PH
step.res <- ordiR2step(mod0, mod1, direction = "both",perm.max = 200)
step.res$anova  # Summary table

g <- autoplot(spe.rda, arrows = FALSE, geom = c("point", "text")) + geom_text(size = 10) + theme_few() + xlim(c(-0.7,1.2)) + scale_color_manual(values = c("black","#009999","#006666"))
g + geom_hline(yintercept = 0, linetype="dotted") + geom_vline(xintercept = 0, linetype="dotted")
ggsave("environmentspecies.png", units="in", width=10, height=10, dpi=300)

# Spatial effects
spe.rda <- dbrda(meandist ~ spavar$netcen + spavar$updist)

plot(spe.rda, scaling = 1)
summary(spe.rda)
anova(spe.rda)
anova.cca(spe.rda, step=1000);
RsquareAdj(spe.rda)$adj.r.squared;
RsquareAdj(spe.rda)$r.squared


spe.rda <- dbrda(spe.hel.avab ~ spavar$netcen + spavar$updist)

plot(spe.rda, scaling = 1)
summary(spe.rda)
anova(spe.rda)
anova.cca(spe.rda, step=1000);
RsquareAdj(spe.rda)$adj.r.squared;
RsquareAdj(spe.rda)$r.squared

#variance partitioning
spe.varpart1 <- varpart(meandist, cbind(spavar[,2:3],spa.PCNM), env_select)
spe.varpart1 <- varpart(meandist, cbind(spa.PCNM[,1:2]), env_select)
spe.varpart1 <- varpart(meandist, cbind(spavar[,2:3]), env_select)
par(mfrow=c(1,2))
showvarparts(2)
plot(spe.varpart1,digits=2)
spe.varpart1

spe.varpart1 <- varpart(spe.hel.avab, cbind(spavar[,2:3],spa.PCNM), env_select)
spe.varpart1 <- varpart(spe.hel.avab, cbind(spa.PCNM[,1:2]), env_select)
spe.varpart1 <- varpart(spe.hel.avab, cbind(spavar[,2:3]), env_select)
par(mfrow=c(1,2))
showvarparts(2)
plot(spe.varpart1,digits=2)
spe.varpart1


#cor(cbind(spavar[,2:3],spa.PCNM))


plot(avab[,2], env_select$T_max)
avin

prev = aggregate(data[,c(23:25,27:33)], by = list(data[,1]), function(x){sum(x >0, na.rm = T)/length(x)})

# Hellinger transformation of the species dataset ####
spe.hel <- vegdist(decostand(avab[,-1],"hellinger"), method="bray") # Bray-Curtis distances

env_select <- env[,c("T_max","con_av","O2_sat_av","Cl_av","COD_av","NH4_av","NO3_av","NO2_av")]

spe.rda <- dbrda(spe.hel ~ T_max + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av, env_select)
plot(spe.rda, scaling = 1)
summary(spe.rda)
anova(spe.rda)
anova.cca(spe.rda, step=1000);
RsquareAdj(spe.rda)$adj.r.squared;
RsquareAdj(spe.rda)$r.squared




















# Set working directory
setwd('C:/Users/pascalh/Documents/GitHub/Stickleback-parasites-2016')

# Load data
data_2016 <- read.csv("data_2016.csv", sep=';') #field and parasite data
env <- read.csv("Environment_R.csv", sep=',') #all environmental data
spavar <- read.csv("space2.csv", sep=';') #spatial variables: network centrality and upstream distance

# make a new data frame by combinding parasite, environmental and space data
env_exp <- env %>% slice(rep(1:n(), table(as.factor(data_2016$site))))
spavar_exp <- spavar %>% slice(rep(1:n(), table(as.factor(data_2016$site))))

data <- cbind(data_2016, env_exp[,-1], spavar_exp)
data$site <- as.factor(data$site)
data$fish <- as.factor(data$fish)
data$length <- as.numeric(data$length)

data0 <- cbind(data_2016, env_exp[,-1], spavar_exp, spa.PCNM_exp)
NAs <- 1>rowSums(is.na(data0[,c(23:25,27:33)]))
table(NAs)
data <- data0[NAs,] # remove fish with missing data

spe.hel <- vegdist(decostand(cbind(data[,c(23:25,27:33)],rep(1,nrow(data))),"hellinger", na.rm=T), "bray", na.rm=T)

hist(as.matrix(data2[,c(23,24,25,27,28,29,30,31,32,33)]), breaks=2000)

dummy_env <- matrix(0, nrow(data2), length(levels(data2$site)))
for(i in 1:length(levels(data2$site))){
  dummy_env[,i] <- ifelse(data2$site == levels(data2$site)[i], 1, 0)
}

env_select <- cbind(data[,c("T_av","con_av","O2_sat_av","Cl_av","COD_av","NH4_av","NO3_av","NO2_av")])

#### Environmental variables ####
#spe.rda <- rda(spe.hel, env_select)
spe.rda <- dbrda(spe.hel ~ T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av, env_select)
plot(spe.rda, scaling = 1)
summary(spe.rda)
anova(spe.rda)
anova.cca(spe.rda, step=1000);
RsquareAdj(spe.rda)$adj.r.squared;
RsquareAdj(spe.rda)$r.squared


# Spatial variables  ####
spa.rda <- dbrda(spe.hel ~ data$netcen + data$updist)
summary(spa.rda)
plot(spa.rda)
anova(spa.rda)
anova.cca(spa.rda, step=1000)
RsquareAdj(spa.rda)$r.squared;
RsquareAdj(spa.rda)$adj.r.squared;

#variance partitioning
spe.varpart1 <- varpart(spe.hel, cbind(data$netcen + data$updist), env_select)
par(mfrow=c(1,2))
showvarparts(2)
plot(spe.varpart1,digits=2)
spe.varpart1




#setwd('C:/Users/u0113095/Google Drive/PhD/2 Parasite/2016/Parasite2016_analysis') #edit_PH

ab <- read.csv("abundance.csv", sep=";")
ab <- ab[,-1]
pre <- read.csv("prevalence.csv", sep=';')

spavar <- read.csv("space2.csv", sep=';') #spatial variables
spavar <- spavar[,-1]
spa <- read.csv("space.csv", sep=';')
plot(spa$long ~ spa$lat)
dist <- pointDistance(spa, allpairs = TRUE, lonlat = FALSE)

env_max <- read.csv("env_max.csv", sep=';')
env_av <- read.csv("env_av.csv", sep=';') #edit_PH

env_all <- cbind(env_max, env_av) #edit_PH

cor(env_all) #edit_PH

env_select <- (env_all[,c(4,6,8,10:12)]) #edit_PH

env <- read.csv("Environment_R.csv", sep=',') #all environmental data
env_select <- (env[,c("T_av","con_av","O2_sat_av","Cl_av","COD_av","NH4_av","NO3_av","NO2_av")]) #edit_PH

# Hellinger transformation of the species dataset ####
spe.hel <- decostand(pre,"hellinger") #prevalence or abundance

#### Environmental variables ####
spe.rda <- rda(spe.hel, env_max)
plot(spe.rda, scaling = 1)
summary(spe.rda)
anova(spe.rda)
anova.cca(spe.rda, step=1000);
RsquareAdj(spe.rda)$adj.r.squared;
RsquareAdj(spe.rda)$r.squared

mod0 <- rda(spe.hel ~ 1, env_select)  # Model with intercept only  #edit_PH
mod1 <- rda(spe.hel ~ ., env_select)  # Model with all explanatory variables  #edit_PH
step.res <- ordiR2step(mod0, mod1, direction = "backward",perm.max = 200)
step.res$anova  # Summary table

g <- autoplot(spe.rda, arrows = FALSE, geom = c("point", "text")) + geom_text(size = 10) + theme_few() + xlim(c(-0.7,1.2)) + scale_color_manual(values = c("black","#009999","#006666"))
g + geom_hline(yintercept = 0, linetype="dotted") + geom_vline(xintercept = 0, linetype="dotted")
ggsave("environmentspecies.png", units="in", width=10, height=10, dpi=300)

# Spatial variables  ####
distance_matrix <- read.csv("distance_matrix.csv", sep=';')
KautoDist <- data.frame(cmdscale(distance_matrix),rownames(distance_matrix));colnames(distance_matrix)<- c("X1","X2","LocationID")
spa <- KautoDist[,c(1,2)]
library(raster)
dist <- pointDistance(spa, allpairs = TRUE, lonlat = FALSE)
spa.dist <- as.dist(dist)

# Test and forward selectio of PCNM variables
spanning <- spantree(spa.dist)
dmin <- max(spanning$dist)
spa.dist[spa.dist>dmin] <- 4*dmin
xy.PCoA <- cmdscale(spa.dist,k=nrow(spa)-1,eig=TRUE)
nb.ev<-length(which(xy.PCoA$eig > 0.0001))
xy.PCoA$eig

# Broken stick model
ev <- xy.PCoA$eig
ev <- subset(ev, xy.PCoA$eig > 0.0001)
n <- length(ev)
bsm <- data.frame(j=seq(1:n), p=0)
bsm$p[1] <- 1/n
for (i in 2:n) {
  bsm$p[i] = bsm$p[i-1] + (1/(n+1-i))
}
bsm$p <- 100*bsm$p/n
bsm
par(mfrow=c(1,1))
barplot(ev, main="Eigenvalues", col="bisque", las=2)
abline(h=mean(ev), col="red")
legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")
barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=TRUE, main="% variance", col=c("bisque", 2), las=2)
nb.ev <- length(which(xy.PCoA$eig > mean(ev)))

spa.PCNM <- as.data.frame(xy.PCoA$points[1:nrow(spa),1:2])

# Distance matrix  ####
spa.rda <- rda(spe.hel,spa.PCNM)
summary(spa.rda)
plot(spa.rda)
anova(spa.rda)
anova.cca(spa.rda, step=1000)
R2 <- RsquareAdj(spa.rda)$r.squared;
R2
R2 <- RsquareAdj(spa.rda)$adj.r.squared;
R2

mod0 <- rda(spe.hel ~ 1, spa.PCNM)  # Model with intercept only  #edit_PH
mod1 <- rda(spe.hel ~ ., spa.PCNM)  # Model with all explanatory variables  #edit_PH
step.res <- ordiR2step(mod0, mod1, direction = "backward",perm.max = 200)
step.res$anova  # Summary table


g <- autoplot(spa.rda, arrows = FALSE, geom = c("point", "text")) + geom_text(size = 16) + theme_few() + scale_color_manual(values = c("black","darkgrey","black"))
g + geom_hline(yintercept = 0, linetype="dotted") + geom_vline(xintercept = 0, linetype="dotted") 
ggsave("environmentspace2.png", units="in", width=10, height=10, dpi=300)

# Spatial variables  ####
spa.rda <- rda(env_select,spa.PCNM)  #edit_PH
summary(spa.rda)
plot(spa.rda)
anova(spa.rda)
anova.cca(spa.rda, step=1000)
R2 <- RsquareAdj(spa.rda)$r.squared;
R2 <- RsquareAdj(spa.rda)$adj.r.squared;
R2

mod0 <- rda(spe.hel ~ 1, space)  # Model with intercept only
mod1 <- rda(spe.hel ~ ., space)  # Model with all explanatory variables
step.res <- ordiR2step(mod0, mod1, direction = "backward",perm.max = 200)
step.res$anova  # Summary table


g <- autoplot(spa.rda, arrows = FALSE, geom = c("point", "text")) + geom_text(size = 16) + theme_few() + scale_color_manual(values = c("black","darkgrey","black"))
g + geom_hline(yintercept = 0, linetype="dotted") + geom_vline(xintercept = 0, linetype="dotted") 
ggsave("environmentspace2.png", units="in", width=10, height=10, dpi=300)

#variance partitioning
spe.varpart1 <- varpart(spe.hel,space, env)
par(mfrow=c(1,2))
showvarparts(2)
plot(spe.varpart1,digits=2)
spe.varpart1




