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
setwd('C:/Users/u0113095/Documents/GitHub/Stickleback-parasites-2016')
# Load data
data_2016 <- read.csv("data_2016_1211.csv", sep=';') #field and parasite data
env <- read.csv("Environment_RDA.csv", sep=';') #all environmental data
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


#### CALCULATE INFECTION PARAMETERS ####
names(data)
parsum = aggregate(data[,c(21:23,25:31)], by = list(data[,1]), function(x){sum(x, na.rm = T)}) 
avin = aggregate(data[,c(21:23,25:31)], by = list(data[,1]), function(x){mean(x[x >0], na.rm = T)}); avin[is.na(avin)] <- 0
avab = aggregate(data[,c(21:23,25:31)], by = list(data[,1]), function(x){mean(x, na.rm =T)})
prev = aggregate(data[,c(21:23,25:31)], by = list(data[,1]), function(x){sum(x >0, na.rm = T)/length(x)})
medin = aggregate(data[,c(21:23,25:31)], by = list(data[,1]), function(x){median(x[x >0], na.rm = T)}) ; medin[is.na(medin)] <- 0

#parasite data is overdispersed (mostly so for Trichodina), if using average abundance data, species matrix needs to be transformed
datao <- na.omit(data[,c(21:23,25:31)])
ddata <- dispweight(datao[,-1])
avab <- aggregate(ddata, by = list(datao[,1]), function(x){mean(x, na.rm =T)})

###########################
# RDA on infracommunities #
###########################

# Infracommunities: Bray-Curtis dissimilarities are calculated at the individual host level Hellinger-transformed parasite data and then averaged within site
# A dummy parasite species is added to avoid problems with non-infected fishes
data_infra <- na.omit(data[,c(1,21:23,25:31)])
data_infra_disp <- dispweight(data_infra[,-1])
braycurtis <- vegdist(decostand(cbind(data_infra_disp,rep(1,nrow(data_infra))), na.rm=T, method="hellinger"), method="bray", na.rm=T)
meandist_bray <- meandist(braycurtis, data_infra[,1])

# Check whether Euclidean and Bray-Curtis distances are comparable
braycurtis <- vegdist(decostand(cbind(data_infra_disp,rep(1,nrow(data_infra))), na.rm=T, method="hellinger"), method="bray", na.rm=T)
meandist_bray <- meandist(braycurtis, data_infra[,1])
euc <- vegdist(decostand(cbind(data_infra_disp,rep(1,nrow(data_infra))), na.rm=T, method="hellinger"), method="euc", na.rm=T)
meandist_euc <- meandist(euc, data_infra[,1])
plot(meandist_bray[1:37,1:37], meandist_euc[1:37,1:37])
mantel(meandist_bray[1:37,1:37], meandist_euc[1:37,1:37])


# environmental variables
env_select <- env[,c("Temperature","Conductivity","Oxygen","COD","NH4","Nt","Meander", "Poolriffle")]
env_select$Poolriffle <- as.factor(env_select$Poolriffle)
env_select$Meander <- as.factor(env_select$Meander)

# Assess the effect of environmental variables on parasite infracommunity dissimilarities using distance based RDA
spe.rda <- dbrda(meandist_bray ~ Temperature + Conductivity + Oxygen + COD + NH4 + Nt + Meander +
                   Poolriffle, env_select)
mod0 <- dbrda(meandist_bray ~ 1, env_select)  # Model with intercept only  #edit_PH
mod1 <- dbrda(meandist_bray ~ ., env_select)  # Model with all explanatory variables  #edit_PH
step.res <- ordiR2step(mod0, mod1, direction = "both",perm.max = 200)
step.res$anova
step.res$anova  # Summary table

spe.rda <- dbrda(meandist_bray ~  COD + NH4 + Meander, env_select)
plot(spe.rda, scaling = 1) # it is for technical reasons not possible to plot both site and species scores
summary(spe.rda)
anova(spe.rda)
anova(spe.rda, by="term")
anova.cca(spe.rda, step=1000);
anova.cca(spe.rda, step=1000, by="term");
RsquareAdj(spe.rda)$adj.r.squared;
RsquareAdj(spe.rda)$r.squared


# Check whether exclusion of site 12 changes results
spe.rda <- dbrda(meandist_bray[-10,-10] ~ Temperature + Conductivity + Oxygen + COD + NH4 + Nt + Meander +
                   Poolriffle, env_select[-10,])

plot(spe.rda, scaling = 1) # it is for technical reasons not possible to plot both site and species scores
summary(spe.rda)
anova(spe.rda)
anova(spe.rda, by="term")
anova.cca(spe.rda, step=1000);
anova.cca(spe.rda, step=1000, by="term");
RsquareAdj(spe.rda)$adj.r.squared;
RsquareAdj(spe.rda)$r.squared

mod0 <- dbrda(meandist_bray[-10,-10] ~ 1, env_select[-10,])  # Model with intercept only  #edit_PH
mod1 <- dbrda(meandist_bray[-10,-10] ~ ., env_select[-10,])  # Model with all explanatory variables  #edit_PH
step.res <- ordiR2step(mod0, mod1, direction = "both",perm.max = 200)
step.res$anova  # Summary table

png(file="Results/Infra_env.png", res=600, width=3000, height=3000)
pRDAplot <- plot(spe.rda, choices = c(1, 2), type="n", cex.lab=1, xlab="dbRDA 1 (13.20%) ", ylab="dbRDA2 (3.31%)")
with(env_select, points(spe.rda, display = "sites", cex=0.7, pch = 19, col = 'black'))
text(spe.rda, "bp",choices = c(1, 2), col="black", cex=0.6)
dev.off()

# Spatial effects
space <- cbind(spavar[,c(2,3)], spa.PCNM[,3])
space$Upstream <- space$updist
space$NetCen <- space$netcen
space$PCNM3 <- space$`spa.PCNM[, 3]`
space <- space[,c(4:6)]

spe.rda <- dbrda(meandist_bray~Upstream + NetCen + PCNM3, space)

plot(spe.rda, scaling = 1)
summary(spe.rda)
anova(spe.rda)
anova.cca(spe.rda, step=1000);
anova.cca(spe.rda, step=1000, by="term");
RsquareAdj(spe.rda)$adj.r.squared;
RsquareAdj(spe.rda)$r.squared

png(file="Results/Infra_spa.png", res=600, width=3000, height=3000)
pRDAplot <- plot(spe.rda, choices = c(1, 2), type="n", cex.lab=1, xlab="dbRDA 1 (4.15%) ", ylab="dbRDA2 (3.70%)")
with(space, points(spe.rda, display = "sites", cex=0.7, pch = 19, col = 'black'))
text(spe.rda, "bp",choices = c(1, 2), col="black", cex=0.6)
dev.off()


#Variation partitioning
spe.varpart1 <- varpart(meandist_bray, env_select, space)
plot(spe.varpart1,digits=2)
spe.varpart1

anova.cca(dbrda(meandist_bray ~ space[,1] + space[,2] + space[,3] +
                  Condition(Temperature + Conductivity + Oxygen + COD +
                              NH4 + Nt + Meander + Poolriffle), data=env_select), step=1000)

anova.cca(dbrda(meandist_bray ~ Temperature + Conductivity + Oxygen + COD +
                  NH4 + Nt + Meander + Poolriffle + Condition(space[,1] + space[,2] + space[,3]),
                data=env_select), step=1000)

# PermANOVA (adonis)
adonis(meandist_bray~ . + space[,1] + space[,2] + space[,3] , data=env_select)

#Removing location 12
#Variation partitioning
spe.varpart1 <- varpart(meandist_bray[-10,-10], cbind(spavar[-10,2:3]), env_select[-10,])
par(mfrow=c(1,2))
showvarparts(2)
plot(spe.varpart1,digits=2)
spe.varpart1
# Unique space fraction
anova_space <- anova.cca(dbrda(meandist_bray[-10,-10] ~ spavar[-10,2] + spavar[-10,3] + Condition(T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av), data=env_select[-10,]), step=1000)
anova_space
# Unique environment fraction
anova_env <- anova.cca(dbrda(meandist_bray[-10,-10] ~ . + Condition(spavar[-10,2] + spavar[-10,3]), data=env_select[-10,]), step=1000)
anova_env

# PermANOVA (adonis)
adonis(meandist_bray[-10,-10] ~ . + spavar[-10,2] + spavar[-10,3], data=env_select[-10,])
adonis(meandist_bray[-10,-10] ~ con_av + spavar[-10,"updist"], data=env_select[-10,])




################################
# RDA on component communities #
################################

# Component communities: Bray-Curtis dissimilarities based on Hellinger transformed average abundance data
#parasite data is overdispersed (mostly so for Trichodina), if using average abundance data, species matrix needs to be transformed
datao <- na.omit(data[,c(1,21:23,25:31)])
ddata <- dispweight(datao[,-1])
avab <- aggregate(ddata, by = list(datao[,1]), function(x){mean(x, na.rm =T)})

spe.hel_bray <- vegdist(decostand(avab[,-1], na.rm=T, method="hellinger"), method="bray", na.rm=T)

# Check whether Euclidean and Bray-Curtis distances are comparable
spe.hel_euc <- vegdist(decostand(avab[,-1], na.rm=T, method="hellinger"), method="euc", na.rm=T)
plot(spe.hel_bray, spe.hel_euc)
mantel(spe.hel_bray, spe.hel_euc)

# environmental variables
env_select <- env[,c("T_av","con_av","O2_sat_av","Cl_av","COD_av","NH4_av","NO3_av","NO2_av")]

# Assess the effect of environmental variables on parasite infracommunity dissimilarities using distance based RDA
spe.rda <- dbrda(spe.hel_bray ~ Temperature + Conductivity + Nt + Meander +
                   Poolriffle, env_select)
sppscores(spe.rda) <- decostand(avab[,-1], na.rm=T, method="hellinger") 

plot(spe.rda, scaling = 1) # it is for technical reasons not possible to plot both site and species scores
summary(spe.rda)
anova(spe.rda)
anova(spe.rda, by="term")
anova.cca(spe.rda, step=1000);
anova.cca(spe.rda, step=1000, by="term");
RsquareAdj(spe.rda)$adj.r.squared;
RsquareAdj(spe.rda)$r.squared

mod0 <- dbrda(spe.hel_bray ~ 1, env_select)  # Model with intercept only  #edit_PH
mod1 <- dbrda(spe.hel_bray ~ Temperature + Conductivity + Nt + Meander +
                Poolriffle, env_select)  # Model with all explanatory variables  #edit_PH
step.res <- ordiR2step(mod0, mod1, direction = "both",perm.max = 200)
step.res$anova  # Summary table


png(file="Results/Comp_env.png", res=600, width=3000, height=3000)
pRDAplot <- plot(spe.rda, choices = c(1, 2), type="n", cex.lab=1, xlab="dbRDA 1 (9.96%) ", ylab="dbRDA2 (6.65%)")
with(env_select, text(spe.rda, display = "species", cex=0.7, pch = 19, col = 'grey'))
with(env_select, points(spe.rda, display = "sites", cex=0.7, pch = 19, col = 'black'))
text(spe.rda, "bp",choices = c(1, 2), col="black", cex=0.6)
dev.off()

# Spatial effects
spe.rda <- dbrda(spe.hel_bray ~ Upstream + NetCen + PCNM3, space)
sppscores(spe.rda) <- decostand(avab[,-1], na.rm=T, method="hellinger")
plot(spe.rda, scaling = 2)
summary(spe.rda)
anova(spe.rda)
anova.cca(spe.rda, step=1000);
anova.cca(spe.rda, step=1000, by='term');
RsquareAdj(spe.rda)$adj.r.squared;
RsquareAdj(spe.rda)$r.squared

png(file="Results/Comp_spa.png", res=600, width=3000, height=3000)
pRDAplot <- plot(spe.rda, choices = c(1, 2), type="n", cex.lab=1, xlab="dbRDA 1 (9.96%) ", ylab="dbRDA2 (6.65%)")
with(space, text(spe.rda, display = "species", cex=0.7, pch = 19, col = 'grey'))
with(space, points(spe.rda, display = "sites", cex=0.7, pch = 19, col = 'black'))
text(spe.rda, "bp",choices = c(1, 2), col="black", cex=0.6)
dev.off()

#Variation partitioning
# Temperature + Conductivity + Nt + Meander +Poolriffle
spe.varpart1 <- varpart(spe.hel_bray, space, env_select[,c(1,2,6,7,8)])
plot(spe.varpart1,digits=2)
spe.varpart1
anova.cca(dbrda(spe.hel_bray ~ space[,1] + space[,2] + space[,3] +
                  Condition(Temperature + Conductivity + Oxygen + Nt + Meander 
                            + Poolriffle), data=env_select), step=1000)
anova.cca(dbrda(spe.hel_bray ~ Temperature + Conductivity + Oxygen +
                  Nt + Meander + Poolriffle + Condition(space[,1] + space[,2] + space[,3]),
                data=env_select), step=1000)


# PermANOVA (adonis)
adonis(spe.hel_bray ~ Temperature + Conductivity + Nt + Meander +Poolriffle + spavar[,1]+  spavar[,2] + spavar[,3], data=env_select)





install.packages("lmPerm")
library(lmPerm)


summary(model <- aovp(diag(meandist_bray)~.,env_select))
step.model <- stepAIC(model,
                      direction = "both", 
                      trace = FALSE)
step.model
summary(model <- aovp(diag(meandist_bray)~T_av+con_av+COD_av+NH4_av+NO2_av,env_select))

library(car)
scatterplotMatrix(new_env_select<-as.data.frame(cbind(env_select$T_av, (env_select$con_av)^2, (env_select$O2_sat_av)^3, (env_select$Cl_av)^(1/3), env_select$COD_av, (env_select$NH4_av)^(1/3), env_select$NO3_av^(1/3))))
cor(new_env_select)
plot(new_env_select)
summary(model <- lm(diag(meandist_bray)~.,new_env_select))
plot(model)
step.model <- stepAIC(model,
                      direction = "both", 
                      trace = FALSE)
step.model
summary(model <- lm(diag(meandist_bray)~V1+V4,new_env_select))
plot(model)

?scatterplotMatrix
summary(model <- lm(sqrt(avab$gyro)~.,env_select))
step.model <- stepAIC(model,
                      direction = "both", 
                      trace = FALSE)
step.model
summary(model <- lm(avab$gyro~Cl_av,env_select))
plot(model)
summary(model <- lm(avab$gyro[-10]~Cl_av,env_select[-10,]))
plot(model)


summary(model <- lm(sqrt(avab$tricho)~.,env_select))
step.model <- stepAIC(model,
                      direction = "both", 
                      trace = FALSE)
step.model
summary(model <- lm(avab$gyro~Cl_av+COD_av,env_select))
plot(model)
summary(model <- lm(avab$gyro[-10]~Cl_av+COD_av,env_select[-10,]))
plot(model)

summary(model <- lm(sqrt(avab$glugea)~.,env_select))
step.model <- stepAIC(model,
                      direction = "both", 
                      trace = FALSE)
step.model
summary(model <- lm(avab$glugea~O2_sat_av+NH4_av+NO3_av,env_select))
plot(model)

summary(model <- lm(sqrt(avab$contra)~.,env_select))
step.model <- stepAIC(model,
                      direction = "both", 
                      trace = FALSE)
step.model

summary(model <- lm(sqrt(avab$angui)~.,env_select))
step.model <- stepAIC(model,
                      direction = "both", 
                      trace = FALSE)
step.model
summary(model <- lm(avab$glugea~T_av+O2_sat_av+NH4_av+NO2_av,env_select))
plot(model)




summary(model <- lm(avab$gyro~O2_sat_av+Cl_av+COD_av+NO3_av+NO2_av,env_select))

plot(model)

?lmp
res <- resid(lmp(avab$gyro~O2_sat_av+COD_av+NO3_av+NO2_av,env_select))
plot(res~Cl_av, env_select)


summary(model <- lm(diag(meandist_bray)~.,env_select))
summary(model <- rlm(diag(meandist_bray)~.,env_select))
library(MASS)
#install.packages("MuMIn")
library(MuMIn)
dredge(model, rank = "AIC")
#install.packages("sfsmisc")
library(sfsmisc)
f.robftest(model, var = "T_av")

summary(model <- lm(avab$gyro~.,env_select))
summary(model <- rlm(avab$gyro~.,env_select))
dredge(model, rank = "BIC")
summary(model <- rlm(avab$gyro~Cl_av+con_av+NO3_av,env_select))
f.robftest(model, var = "Cl_av")
f.robftest(model, var = "con_av")
f.robftest(model, var = "NO3_av")

plot(model)

step.model <- stepAIC(lm(diag(meandist_bray)~.,env_select),
                      direction = "both", 
                      trace = FALSE)
step.model
summary(lm(diag(meandist_bray)~T_av+con_av+COD_av+NH4_av+NO2_av,env_select))

res_T_av <- resid(lm(diag(meandist_bray)~con_av+COD_av+NH4_av+NO2_av,env_select))
plot(res_T_av~env_select$T_av)
lines(lowess(res_T_av~env_select$T_av), col=3)

res_COD_av <- resid(lm(diag(meandist_bray)~T_av+con_av+NH4_av+NO2_av,env_select))
plot(res_COD_av~env_select$COD_av)
lines(lowess(res_COD_av~env_select$COD_av), col=3)

library(caret)
# Define training control
set.seed(123)
train.control <- trainControl(method = "repeatedcv", 
                              number = 3, repeats = 10)
x <- diag(meandist_bray)
new_data <- cbind(env_select, x)
colnames(new_data)
# Train the model
model <- train(x~T_av+con_av+O2_sat_av+Cl_av+COD_av+NH4_av+NO3_av+NO2_av, data=new_data, method = "lmStepAIC",
               trControl = train.control)
summary(model)
model2 <- lm(x~T_av+con_av+COD_av+NH4_av+NO2_av, data=new_data)
summary(model2)

step.model <- stepAIC(lm(x~T_av+con_av+COD_av+NH4_av+NO2_av, data=new_data),
                      direction = "both", 
                      trace = FALSE)
step.model



#Step Backward and remove one variable at a time
tctrl <- trainControl(method = "cv",number=3,
                      repeats=10)
#Declare exactly which parameters you want to test
rpart_opts <- expand.grid(cp = seq(0.0,0.01, by = 0.001))

rpart_model <- train(x~., data = new_data, method="rpart",
                     metric = "RMSE", trControl = tctrl,
                     tuneGrid = rpart_opts)

summary(rpart_model)

?confusionMatrix
# Summarize the results
print(model)
#plot(model)
?caret::train

install.packages("mlbench")
library(mlbench)
data(BostonHousing)

lmFit <- train(medv ~ . + rm:lstat,
               data = BostonHousing,
               method = "lm",
               class
               metric = "ROC")
lmFit$finalModel
ggplot(lmFit)

###### Old code


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




