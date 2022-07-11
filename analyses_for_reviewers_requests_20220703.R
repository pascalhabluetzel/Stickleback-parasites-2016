# The script runs on R v. 4.1.2

# Empty environment
rm(list=ls())

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # requires installation of package "rstudioapi"

# This script codes the Bayesian Model Averaging for the stickleback parasite study
# It also includes figures

# Load libraries
library(ggplot2) # for making figures
library(ggeffects) # for making figures
library(vegan) # for function dispweight
library(MASS) # for step-wise model selection
library(car) # for VIF
library(BAS) # for Bayesian models
library(dplyr) # utility package
library(nlme) # for linear mixed effect models
library(piecewiseSEM) # for SEM
library(boral) # boral depends on jags for installation, which needs to be installed manually on linux ('sudo apt install jags')

##########################
#### 0. PRELIMINARIES ####
##########################

#### READ AND PREPARE DATA ####
# Parasite data
data <- read.csv("data_2016_2303.csv", sep=';')
data$site <- as.factor(data$site)

# Environmental data (VMM)
environment <- read.csv("Environment_update.csv", sep=';')
spavar <- read.csv("space2.csv", sep=';') #spatial variables: network centrality and upstream distance
plot(spavar$netcen); plot(density(spavar$netcen))
plot(spavar$updist); plot(density(spavar$updist))
#environmental data (from field observations)
field_data <- read.csv("field_data.csv", sep=',')
environment2 <- cbind(environment[,c(1,49,52:53,55,57,60,63)], field_data[-c(8,10,25,27,37),33:34], spavar[,c(2,3)])
environment2$pool_riffle <- as.factor(environment2$pool_riffle)
environment2$meander <- as.factor(environment2$meander)

#### Matrix for PIP (Posterior Inclusion Probability) ####
PIP <- matrix(nrow=12, ncol=14)
rownames(PIP) <- c("Condition", "Length", "Temperature", "Oxygen saturation", "Conductivity", "COD", "Ammonium", "Total nitrogen", "Pool riffle pattern", "Meander", "Network centrality", "Upstream distance")
colnames(PIP) <- c("Condition", "Length", "Gyrodactylus abundance", "Gyrodactylus prevalence", "Gyrodactylus infection intensity", "Trichodina abundance", "Trichodina prevalence", "Trichodina infection intensity", "Glugea", "Contracaecum", "Aguillicola",
                   "PI", "PI_ecto", "PI_endo")  

#### CALCULATE PARASITE PARAMETERS ####
names(data)
#parasite data is overdispersed (mostly so for Trichodina), if using average abundance data, species matrix needs to be transformed
datao <- na.omit(data[,c(1,22:24,26:32)]) #remove fish without parasite counts
ddata <- dispweight(datao[,-1]) #correct for overdispersion of the parasite count data
avab <- aggregate(ddata, by = list(datao[,1]), function(x){mean(x, na.rm =T)})
prev = aggregate(data[,c(22:24,26:32)], by = list(data[,1]), function(x){sum(x >0, na.rm = T)/length(x)})
medin = aggregate(data[,c(22:24,26:32)], by = list(data[,1]), function(x){median(x[x >0], na.rm = T)}) 
pa = aggregate(data[,c(22:24,26:32)], by = list(data[,1]), function(x){ifelse(mean(x, na.rm =T)>0,1,0)}) 

avab[is.na(avab)] <- 0
prev[is.na(prev)] <- 0
medin[is.na(medin)] <- 0

#### PARASITE INDEX ####
PI <- 1:nrow(data)
for(j in 1:nrow(data)){
  PI[j] <- 
    max(data$Gyr[j]/sd(data$Gyr, na.rm=T))*(data$Gyr[j]/sd(data$Gyr, na.rm=T))
  + max(data$Tri[j]/sd(data$Tri, na.rm=T))*(data$Tri[j]/sd(data$Tri, na.rm=T))
  + max(data$Glu[j]/sd(data$Glu, na.rm=T))*(data$Glu[j]/sd(data$Glu, na.rm=T))
  + max(data$Con[j]/sd(data$Con, na.rm=T))*(data$Con[j]/sd(data$Con, na.rm=T))
  + max(data$CysL[j]/sd(data$CysL, na.rm=T))*(data$CysL[j]/sd(data$CysL, na.rm=T))
  + max(data$Pro[j]/sd(data$Pro, na.rm=T))*(data$Pro[j]/sd(data$Pro, na.rm=T))
  + max(data$Aca[j]/sd(data$Aca, na.rm=T))*(data$Aca[j]/sd(data$Aca, na.rm=T))
  + max(data$Cam[j]/sd(data$Cam, na.rm=T))*(data$Cam[j]/sd(data$Cam, na.rm=T))
  + max(data$Ang[j]/sd(data$Ang, na.rm=T))*(data$Ang[j]/sd(data$Ang, na.rm=T))
  + max(data$CisI[j]/sd(data$CisI, na.rm=T))*(data$CisI[j]/sd(data$CysI, na.rm=T))
}

PI_ecto <- 1:nrow(data)
for(j in 1:nrow(data)){
  PI_ecto[j] <- 
    max(data$Gyr[j]/sd(data$Gyr, na.rm=T))*(data$Gyr[j]/sd(data$Gyr, na.rm=T))
  + max(data$Tri[j]/sd(data$Tri, na.rm=T))*(data$Tri[j]/sd(data$Tri, na.rm=T))
  + max(data$Glu[j]/sd(data$Glu, na.rm=T))*(data$Glu[j]/sd(data$Glu, na.rm=T))
}

# BORAL analysis (model-based analysis of multivariate abundance data using Bayesian Markov chain Monte Carlo methods for parameter estimation)
# EXAMPLE to test the code
install.packages("mvabund")
library(mvabund) ## Load a dataset from the mvabund package
data(spider)
y <- spider$abun
X <- scale(spider$x)
n <- nrow(y)
p <- ncol(y)
example_mcmc_control <- list(n.burnin = 10, n.iteration = 100, n.thin = 1)
testpath <- file.path(tempdir(), "jagsboralmodel.txt")
spiderfit_nb <- boral(y, X = X, family = "negative.binomial",
                     ranef.ids = data.frame(region = rep(1:7,each=4)),
                     mcmc.control = example_mcmc_control, model.name = testpath)
spiderfit_nb$ranef.coefs.median
spiderfit_nb$ranef.sigma.median

ranefsplot(object = spiderfit_nb)
coefsplot(covname = "soil.dry", object = spiderfit_nb)

spiderfit_nb <- boral(y, X, family = "negative.binomial", ranef.ids = data.frame(region = rep(1:7,each=4)), lv.control = list(num.lv = 2, type = "independent"), save.model = TRUE)

summary(spiderfit_nb)
plot(spiderfit_nb)
envcors <- get.enviro.cor(spiderfit_nb)
rescors <- get.residual.cor(spiderfit_nb)
library(corrplot)
corrplot(envcors$sig.cor, type = "lower", diag = FALSE, title = "Correlations due to covariates", mar = c(3,0.5,2,1), tl.srt = 45)
corrplot(rescors$sig.cor, type = "lower", diag = FALSE, title = "Residual correlations", mar = c(3,0.5,2,1), tl.srt = 45)

# Real analysis
data$site <- as.factor(data$site)
levels(data$site) <- levels(as.factor(environment2$site))
data_m <- merge(data, environment2, by = "site")
data_all <- na.omit(data_m) 
names(data_all)

avcondition <- aggregate(data$SMI, by = list(data[,1]), function(x){mean(x, na.rm =T)})[,2]
avlength <- aggregate(data$length, by = list(data[,1]), function(x){mean(x, na.rm =T)})[,2]

y <- round(cbind(avab$Gyr, avab$Tri, avab$Glu, avab$Con, avab$Ang))
y <- round(cbind(medin$Gyr, medin$Tri, medin$Glu, medin$Con, medin$Ang))
X <- cbind(avcondition,
           avlength,
           environment2$T_av,
           environment2$O2_sat_av,
           environment2$Con_av,
           environment2$COD_av, 
           environment2$NH4._av, 
           environment2$Nt_av, 
           environment2$netcen, 
           environment2$updist, 
           as.numeric(environment2$pool_riffle), 
           as.numeric(environment2$meander))
colnames(X) <- c("avcondition", "avlength", "T", "O2", "Con", "COD", "NH4", "Nt", "netcen", "updist", "pool_riffle", "meander")

example_mcmc_control <- list(n.burnin = 1000, n.iteration = 10000, n.thin = 1)
testpath <- file.path(tempdir(), "jagsboralmodel.txt")
spiderfit_nb <- boral(y, X = X,
                      family = "negative.binomial",
                      mcmc.control = example_mcmc_control,
                      model.name = testpath,
                      lv.control = list(num.lv = 2, type = "independent"),
                      save.model = TRUE)
spiderfit_nb$ranef.coefs.median
spiderfit_nb$ranef.sigma.median

ranefsplot(object = spiderfit_nb)
plot(spiderfit_nb)
par(mfrow=c(3,4), ask=F)
coefsplot(covname = "avcondition", object = spiderfit_nb) #Condition
coefsplot(covname = "avlength", object = spiderfit_nb) #Length
coefsplot(covname = "T", object = spiderfit_nb) #Temperature
coefsplot(covname = "O2", object = spiderfit_nb) #Oxygen
coefsplot(covname = "Con", object = spiderfit_nb) #Conductivity
coefsplot(covname = "COD", object = spiderfit_nb) #COD
coefsplot(covname = "NH4", object = spiderfit_nb) #NH4
coefsplot(covname = "Nt", object = spiderfit_nb) #Nt
coefsplot(covname = "netcen", object = spiderfit_nb) #netcen
coefsplot(covname = "updist", object = spiderfit_nb) #updist
coefsplot(covname = "pool_riffle", object = spiderfit_nb) #poolriffle
coefsplot(covname = "meander", object = spiderfit_nb) #meander

envcors <- get.enviro.cor(spiderfit_nb)
rescors <- get.residual.cor(spiderfit_nb)
library(corrplot)
corrplot(envcors$sig.cor, type = "lower", diag = FALSE, title = "Correlations due to covariates", mar = c(3,0.5,2,1), tl.srt = 45)
corrplot(rescors$sig.cor, type = "lower", diag = FALSE, title = "Residual correlations", mar = c(3,0.5,2,1), tl.srt = 45)

# Test whether parasite infection depends on the effect of the environment on condition



pca_env <- prcomp(all_data[,c("T_av", "O2_sat_av", "Con_av", "COD_av", "NH4._av", "Nt_av")])
pcaenv <- pca_env$x[, 1]
all_data$pcaenv <- pcaenv

model <- lme(PI ~ pcaenv*SMI, random=~1|site, data=all_data, na.action = na.omit)
model
summary(model)


all_data1 <- na.exclude(all_data[,c("PI", "SMI", "site", "netcen.x", "pcaenv")])
all_data$netcen
colnames(all_data)

psem1 <- psem(
  lme(PI ~ pcaenv + SMI, random=~1|site, data=all_data1),
  lme(SMI ~ pcaenv, random=~1|site, data = all_data1)
)

basisSet(psem1)
summary(psem1, aicc=TRUE)
AIC(psem1, aicc=TRUE)

all_data2 <- na.exclude(all_data[,c("PI_ecto", "Con_av", "length", "SMI", "site", "netcen.x", "updist.x", "pcaenv")])
all_data$netcen
colnames(all_data)

psem2 <- psem(
  lme(PI_ecto ~ pcaenv + SMI, random=~1|site, data=all_data2),
  lme(SMI ~ pcaenv, random=~1|site, data = all_data2)
)

basisSet(psem2)
plot(all_data2$netcen.x, all_data2$pcaenv)
summary(psem2, aicc=TRUE)
AIC(psem2, aicc=TRUE)

all_data3 <- na.exclude(all_data[,c("PI_endo", "SMI", "site", "netcen.x", "pcaenv")])
all_data$netcen
colnames(all_data)

psem3 <- psem(
  lme(PI_endo ~ pcaenv + SMI, random=~1|site, data=all_data3),
  lme(SMI ~ pcaenv, random=~1|site, data = all_data3)
)

summary(psem3, aicc=TRUE)
AIC(psem3, aicc=TRUE)

######################
#### 1. Condition ####
######################

# Effect of environment (average) on host condition
avcondition <- aggregate(data$SMI, by = list(data[,1]), function(x){mean(x, na.rm =T)})[,2]
step.model <- stepAIC(lm(avcondition ~ T_av + O2_sat_av + Con_av^2 + COD_av + 
                           log(NH4._av) + log(Nt_av) + pool_riffle + meander + 
                           spavar$netcen + spavar$updist, data=environment2),
                      direction = "both", 
                      trace = FALSE)
step.model

summary(model <- lm(avcondition ~ T_av + pool_riffle + meander, data=environment2))

plot(model)

vif(model)
res <- resid(lm(avcondition ~ pool_riffle, data=environment2))
boxplot(res ~ environment2$meander)
res <- resid(lm(avcondition ~ meander, data=environment2))
boxplot(res ~ environment2$meander)

# Bayesian Model Averaging
bas.model <- bas.lm(avcondition ~  T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(Nt_av) + 
                      pool_riffle + meander + spavar$netcen + spavar$updist, 
                    data=environment2, prior="JZS")
summary(bas.model)
image(bas.model, rotate=F)
coef.model <- coef(bas.model)
abs(coef.model$postmean)-2*coef.model$postsd > 0
plot(confint(coef.model, parm = 2:11))
confint <- confint(coef.model, parm = 2:11)
write.table(confint, 'condition.txt', sep="\t")
pip <- summary(bas.model)
PIP[c(3:12),1] <- pip[2:11,1]*sign(coef.model$postmean[2:11])
coef.model$postmean[2:11]


#### Length ####
# Effect of environment (average) on host length
avlength <- aggregate(data$length, by = list(data[,1]), function(x){mean(x, na.rm =T)})[,2]
summary(avlength); plot(density(avlength))

model <- lm(avlength ~ T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + 
              log(Nt_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2)
summary(model)
plot(model)
vif(model)

stepAIC(lm(avlength ~ T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + 
             log(Nt_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2),
        direction = "both", 
        trace = FALSE)
summary(model <- lm(avlength ~ Con_av + log(NH4._av) + 
                      log(Nt_av) + spavar$netcen, 
                    data=environment2))
plot(model)
vif(model)
ggplot(environment2, aes(avlength, Con_av) ) +
  geom_point() +
  stat_smooth(method = lm, formula = y ~ poly(x, 2, raw = TRUE))


# Bayesian Model Averaging
bas.model <- bas.lm(avlength ~ T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + 
                      log(Nt_av) + pool_riffle + meander + spavar$netcen + spavar$updist, 
                    data=environment2, prior="JZS")
summary(bas.model)
image(bas.model, rotate=F)
coef.model <- coef(bas.model)
abs(coef.model$postmean)-2*coef.model$postsd > 0
plot(confint(coef.model, parm = 2:11))
confint <- confint(coef.model, parm = 2:11)
write.table(confint, 'length.txt', sep="\t")

pip <- summary(bas.model)
PIP[c(3:12),2] <- pip[2:11,1]*sign(coef.model$postmean[2:11])

#Figure
BPM <- predict(bas.model, estimator = "BPM")
variable.names(BPM)
netcen <- spavar$netcen
updist <- spavar$updist
fit <- lm(avlength ~ O2_sat_av + Con_av^2 + netcen + updist, data=environment2)
predict <- ggpredict(fit, terms = "Con_av")

png(file="figure.png", res=600, width=3000, height=3000)
par(mfrow=c(3,3))
ggplot(predict, aes(x, predicted)) +
  theme_bw() +
  geom_line(color="red", size=1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
  geom_point(data = environment2, aes(x=environment2$Con_av, y=avlength)) +
  labs(x=expression("Conductivity ["*mu*"S/cm]"), y=expression("Average length [mm SL]")) +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12)) +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

predict <- ggpredict(fit, terms = "netcen")
png(file="Results/Length_netcen.png", res=600, width=3000, height=3000)
ggplot(predict, aes(x, predicted)) +
  theme_bw() +
  geom_line(color="red") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
  geom_point(data = environment2, aes(x=netcen, y=avlength)) +
  labs(x=expression("Network centrality"), y=expression("Average length [mm SL]")) + 
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12)) +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()


#### Gyrodactylus ####
#### Average abundance ####
plot(avab$Gyr); plot(density(avab$Gyr))
plot((avab$Gyr)^(1/3)); plot(density((avab$Gyr)^(1/3)))
model <- lm((avab$Gyr)^(1/3) ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av +
              log(NH4._av) + log(Nt_av) + pool_riffle + meander + spavar$netcen + spavar$updist,
            data=environment2)
summary(model)
plot(model)                      
step.model <- stepAIC(lm((avab$Gyr)^(1/3) ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av +
                           log(NH4._av) + log(Nt_av) + pool_riffle + meander + spavar$netcen + spavar$updist
                         , data=environment2),
                      direction = "both", 
                      trace = FALSE)
step.model

model <- lm((avab$Gyr)^(1/3) ~ avlength + T_av + COD_av + log(Nt_av) + pool_riffle, data=environment2)
summary(model)
plot(model) 

# Bayesian Model Averaging
bas.model <- bas.glm(avab$Gyr ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av 
                     + log(NH4._av) + log(Nt_av) + pool_riffle + meander + spavar$netcen + 
                       spavar$updist, data=environment2, prior="JZS")
summary(bas.model)
image(bas.model, rotate=F)
coef.model <- coef(bas.model)
abs(coef.model$postmean)-2*coef.model$postsd > 0
confint(coef.model)
plot(confint(coef.model, parm = 2:13))
confint <- confint(coef.model, parm = 2:13)
write.table(confint, 'GyroAA.txt', sep="\t")
pip <- summary(bas.model)
PIP[c(1:12),3] <- pip[2:13,1]*sign(coef.model$postmean[2:13])

#Figure
BPM <- predict(bas.model, estimator = "BPM")
variable.names(BPM)
netcen <- spavar$netcen
updist <- spavar$updist
fit <- lm(avab$gyro ~ avlength + O2_sat_av + Con_av^2 + COD_av + log1p(NH4._av) + netcen + updist, data=environment2)
predict <- ggpredict(fit, terms = "Con_av")

png(file="Results/GyroAA_conductivity.png", res=600, width=3000, height=3000)
ggplot(predict, aes(x, predicted)) +
  theme_bw() +
  geom_line(color="red") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
  geom_point(data = environment2, aes(x=environment2$Con_av, y=avab$gyro)) +
  labs(x=expression("Conductivity ["*mu*"S/cm]"), y=expression(italic(Gyrodactylus)~"sp. [average abundance]"))+
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12)) +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

predict <- ggpredict(fit, terms = "COD_av")
png(file="Results/GyroAA_COD.png", res=600, width=3000, height=3000)
ggplot(predict, aes(x, predicted)) +
  theme_bw() +
  geom_line(color="red") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
  geom_point(data = environment2, aes(x=environment2$COD_av, y=avab$gyro)) +
  labs(x=expression("Chemical oxygen demand [mg/L]"), y=expression(italic(Gyrodactylus)~"sp. [average abundance]")) +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12)) +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

#### Prevalence ####
plot(prev$Gyr); plot(density(prev$Gyr))

model <- lm(prev$Gyr ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + 
              log(NH4._av) + log(Nt_av) + pool_riffle + meander + spavar$netcen + spavar$updist, 
            data=environment2)
summary(model)
plot(model)                      
step.model <- stepAIC(lm(prev$Gyr ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + 
                           log(NH4._av) + log(Nt_av) + pool_riffle + meander + spavar$netcen + spavar$updist, 
                         data=environment2),
                      direction = "both", 
                      trace = FALSE)
step.model

model <- lm(prev$Gyr ~ Con_av + COD_av + pool_riffle + meander, data=environment2)
summary(model)
plot(model)  

# Have a look at a Bayesian approach
bas.model <- bas.lm(prev$Gyr ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + 
                      log(NH4._av) + log(Nt_av) + pool_riffle + meander + spavar$netcen + spavar$updist, 
                    data=environment2, prior="JZS")
summary(bas.model)
image(bas.model, rotate=F)
coef.model <- coef(bas.model)
abs(coef.model$postmean)-2*coef.model$postsd > 0
confint(coef.model)
plot(confint(coef.model, parm = 2:13))
confint <- confint(coef.model, parm = 2:13)
write.table(confint, 'GyroPrev.txt', sep="\t")
pip <- summary(bas.model)
PIP[c(1:12),4] <- pip[2:13,1]*sign(coef.model$postmean[2:13])

#Figure
BPM <- predict(bas.model, estimator = "BPM")
variable.names(BPM)
netcen <- spavar$netcen
updist <- spavar$updist
fit <- lm(prev$gyro ~ Con_av^2 + log1p(NH4._av) + log(Nt_av) + updist, data=environment2)
predict <- ggpredict(fit, terms = "Con_av")

png(file="Results/GyroP_conductivity.png", res=600, width=3000, height=3000)
ggplot(predict, aes(x, predicted)) +
  theme_bw() +
  geom_line(color="red") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
  geom_point(data = environment2, aes(x=environment2$Con_av, y=prev$gyro)) +
  labs(x=expression("Conductivity ["*mu*"S/cm]"), y=expression(italic(Gyrodactylus)~"sp. [prevalence]")) + 
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12)) +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

dev.off()

#### Median abundance ####
plot(medin$Gyr); plot(density(medin$Gyr))
plot((medin$Gyr)^(1/3)); plot(density((medin$Gyr)^(1/3)))

model <- lm((medin$Gyr)^(1/3) ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + 
              log(NH4._av) + log(Nt_av) + pool_riffle + meander + spavar$netcen + spavar$updist,
            data=environment2)
summary(model)
#plot(model)                      

step.model <- stepAIC(lm((medin$Gyr)^(1/3) ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + 
                           COD_av + log(NH4._av) + log(Nt_av) + log(SM_av) + pool_riffle + meander + 
                           spavar$netcen + spavar$updist, data=environment2),
                      direction = "both", 
                      trace = FALSE)
step.model

model <- lm((medin$Gyr)^(1/3) ~ avlength + T_av + COD_av + pool_riffle, data=environment2)
summary(model)
#plot(model)                      

# Bayesian Model Averaging
bas.model <- bas.lm((medin$Gyr)^(1/3) ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + 
                      COD_av + log(NH4._av) + log(Nt_av) + log(SM_av) + pool_riffle + meander + 
                      spavar$netcen + spavar$updist, data=environment2, prior="JZS")
summary(bas.model)
image(bas.model, rotate=F)
coef.model <- coef(bas.model)
abs(coef.model$postmean)-2*coef.model$postsd > 0
confint(coef.model)
plot(confint(coef.model, parm = 2:13))
confint <- confint(coef.model, parm = 2:13)
write.table(confint, 'GyroMedIn.txt', sep="\t")

pip <- summary(bas.model)
PIP[c(1:12),5] <- pip[2:13,1]*sign(coef.model$postmean[2:13])


#### Trichodina ####
#### Average abundance ####
plot(avab$Tri); plot(density(avab$Tri))
plot((avab$Tri)^(1/3)); plot(density((avab$Tri)^(1/3)))

model <- lm((avab$Tri)^(1/3) ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + 
              COD_av + log(NH4._av) + log(Nt_av) + log(SM_av) + pool_riffle + meander + 
              spavar$netcen + spavar$updist, data=environment2)
summary(model)
#plot(model)                      

step.model <- stepAIC(lm((avab$Tri)^(1/3) ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + 
                           COD_av + log(NH4._av) + log(Nt_av) + log(SM_av) + pool_riffle + meander + 
                           spavar$netcen + spavar$updist, data=environment2),
                      direction = "both", 
                      trace = FALSE)
step.model
model <- lm((avab$Tri)^(1/3) ~ Con_av + COD_av + pool_riffle + meander + spavar$netcen,
            data=environment2)
summary(model)
plot(model) 

res <- resid(lm((avab$tricho)^(1/3) ~ COD_av + pool_riffle + meander + spavar$netcen, data=environment2))
plot(res ~ environment2$Con_av)
lines(lowess(res ~ environment2$Con_av), col=3)
cor.test(res, environment2$Con_av, method="spearman")
res <- resid(lm((avab$tricho)^(1/3) ~ Con_av + pool_riffle + meander + spavar$netcen, data=environment2))
plot(res ~ environment2$COD_av)
lines(lowess(res ~ environment2$COD_av), col=3)
cor.test(res, environment2$COD_av, method="spearman")
res <- resid(lm((avab$tricho)^(1/3) ~ Con_av + COD_av + pool_riffle + meander, data=environment2))
plot(res ~ spavar$netcen)
lines(lowess(res ~ spavar$netcen), col=3)
cor.test(res, spavar$netcen, method="spearman")

# Bayesian Model Averaging
bas.model <- bas.lm((avab$Tri)^(1/3) ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + 
                      COD_av + log(NH4._av) + log(Nt_av) + log(SM_av) + pool_riffle + meander + 
                      spavar$netcen + spavar$updist, data=environment2, prior="JZS")
summary(bas.model)
image(bas.model, rotate=F)
coef.model <- coef(bas.model)
abs(coef.model$postmean)-2*coef.model$postsd > 0
confint(coef.model)
plot(confint(coef.model, parm = 2:13))
confint <- confint(coef.model, parm = 2:13)
write.table(confint, 'TrichoAvAb.txt', sep="\t")

pip <- summary(bas.model)
PIP[c(1:12),6] <- pip[2:13,1]*sign(coef.model$postmean[2:13])

#### Prevalence ####
plot(prev$Tri); plot(density(prev$Tri))

model <- lm(prev$Tri ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + 
              COD_av + log(NH4._av) + log(Nt_av) + log(SM_av) + pool_riffle + meander + 
              spavar$netcen + spavar$updist, data=environment2)
summary(model)
plot(model)                      
step.model <- stepAIC(lm(prev$Tri ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + 
                           COD_av + log(NH4._av) + log(Nt_av) + log(SM_av) + pool_riffle +
                           meander + spavar$netcen + spavar$updist, data=environment2),
                      direction = "both", 
                      trace = FALSE)
step.model

model <- lm(prev$Tri ~ Con_av + log(NH4._av), data=environment2)
summary(model)
plot(model)                      
res <- resid(lm(prev$tricho ~ log(NH4._av), data=environment2))
plot(res ~ environment2$Con_av)
lines(lowess(res ~ environment2$Con_av), col=3)
cor.test(res, environment2$Con_av, method="spearman")
res <- resid(lm(prev$tricho ~ Con_av, data=environment2))
plot(res ~ log(environment2$NH4._av))
lines(lowess(res ~ log(environment2$NH4._av)), col=3)
cor.test(res, environment2$NH4._av, method="spearman")

# Bayesian Model Averaging
bas.model <- bas.lm(prev$Tri ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + 
                      COD_av + log(NH4._av) + log(Nt_av) + log(SM_av) + pool_riffle +
                      meander + spavar$netcen + spavar$updist, data=environment2, prior="JZS")
summary(bas.model)
image(bas.model, rotate=F)
coef.model <- coef(bas.model)
abs(coef.model$postmean)-2*coef.model$postsd > 0
confint(coef.model)
plot(confint(coef.model, parm = 2:12))
confint <- confint(coef.model, parm = 2:13)
write.table(confint, 'TrichoPrev.txt', sep="\t")

pip <- summary(bas.model)
PIP[c(1:12),7] <- pip[2:13,1]*sign(coef.model$postmean[2:13])


#### Median abundance ####
plot(medin$Tri); plot(density(medin$Tri))
plot((medin$Tri)^(1/3)); plot(density((medin$Tri)^(1/3)))

model <- lm((medin$Tri)^(1/3) ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + 
              log(NH4._av) + log(Nt_av) + pool_riffle + meander + spavar$netcen + spavar$updist, 
            data=environment2)
summary(model)
plot(model)                      
step.model <- stepAIC(lm((medin$Tri)^(1/3) ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 
                         + COD_av + log(NH4._av) + log(Nt_av) + pool_riffle + meander + spavar$netcen
                         + spavar$updist, data=environment2),
                      direction = "both", 
                      trace = FALSE)
step.model

model <- lm((medin$Tri)^(1/3) ~ Con_av^2 + COD_av + log(NH4._av) + pool_riffle + meander + 
              spavar$netcen, data=environment2)
summary(model)
#plot(model)                      
res <- resid(lm(prev$tricho ~ COD_av + log(NH4._av) + pool_riffle + meander + spavar$netcen, data=environment2))
plot(res ~ environment2$Con_av)
lines(lowess(res ~ environment2$Con_av), col=3)
cor.test(res, environment2$Con_av^2, method="spearman")

# Bayesian Model Averaging
bas.model <- bas.lm((medin$Tri)^(1/3) ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + 
                      COD_av + log(NH4._av) + log(Nt_av) + pool_riffle + meander + spavar$netcen + 
                      spavar$updist, data=environment2, prior="JZS")
summary(bas.model)
image(bas.model, rotate=F)
coef.model <- coef(bas.model)
abs(coef.model$postmean)-2*coef.model$postsd > 0
confint(coef.model)
plot(confint(coef.model, parm = 2:12))
confint <- confint(coef.model, parm = 2:13)
write.table(confint, 'TrichoMedin.txt', sep="\t")

pip <- summary(bas.model)
PIP[c(1:12),8] <- pip[2:13,1]*sign(coef.model$postmean[2:13])

#Figure
BPM <- predict(bas.model, estimator = "BPM")
variable.names(BPM)
updist <- spavar$updist
fit <- lm((medin$tricho)^(1/3) ~ T_av + Con_av^2 + COD_av + log(NH4._av) + updist, data=environment2)
predict <- ggpredict(fit, terms = "Con_av")

png(file="Results/Tricho_conductivity.png", res=600, width=3000, height=3000)
ggplot(predict, aes(x, predicted^3)) +
  theme_bw() +
  geom_line(color="red") +
  geom_ribbon(aes(ymin = conf.low^3, ymax = conf.high^3), alpha = .1) +
  geom_point(data = environment2, aes(x=environment2$Con_av, y=medin$tricho)) +
  labs(x=expression("Conductivity ["*mu*"S/cm]"), y=expression(italic(Trichodina)~"sp. [median infection intensity]")) +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12)) +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

png(file="Results/Tricho_NH4.png", res=600, width=3000, height=3000)
updist <- spavar$updist
fit <- lm((medin$tricho)^(1/3) ~ T_av + Con_av^2 + COD_av + log1p(NH4._av) + updist, data=environment2) # add a small number to the ammonium value to avoid infinite values in subsequent calculations
predict <- ggpredict(fit, terms = "NH4._av")
ggplot(predict, aes(x, predicted^3)) +
  theme_bw() +
  geom_line(color="red") +
  geom_ribbon(aes(ymin = conf.low^3, ymax = conf.high^3), alpha = .1) +
  geom_point(data = environment2, aes(x=environment2$NH4._av, y=medin$tricho)) +
  labs(x=expression("NH4+"), y=expression(italic(Trichodina)~"sp. [median infection intensity]")) +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12)) +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

#### Glugea ####
model <- glm(pa$Glu ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + 
               log(NH4._av) + log(Nt_av) + log(SM_av) + meander, data=environment2, 
             family=binomial(link="logit"))
summary(model)
#plot(model)   
step.model <- stepAIC(glm(pa$Glu ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + 
                            log(NH4._av) + log(Nt_av) + log(SM_av) + meander, data=environment2, 
                          family=binomial(link="logit")),
                      direction = "both", 
                      trace = FALSE)
step.model

model <- glm(pa$Glu ~ avlength + T_av + O2_sat_av + Con_av^2 + meander, data=environment2, family=binomial(link="logit"))
summary(model)
plot(model) 

plot(as.factor(pa$glugea), avlength)
plot(as.factor(pa$glugea), environment2$Con_av)
table(as.factor(pa$glugea), environment2$meander)

# Bayesian Model Averaging
pa_glugea <- pa$Glu
netcen <- spavar$netcen
updist <- spavar$updist
bas.model <- bas.glm(pa_glugea ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + 
                       log(NH4._av) + log(Nt_av) + pool_riffle + meander + netcen + updist, 
                     data=environment2, family = binomial(link = "logit"), betaprior=bic.prior(), modelprior=uniform())
summary(bas.model)
image(bas.model, rotate=F)
coef.model <- coef(bas.model)
abs(coef.model$postmean)-2*coef.model$postsd > 0
confint(coef.model)
plot(confint(coef.model, parm = 2:12))

confint <- confint(coef.model, parm = 2:13)
write.table(confint, 'Glugea.txt', sep="\t")

pip <- summary(bas.model)
PIP[c(1:12),9] <- pip[2:13,1]*sign(coef.model$postmean[2:13])

#Figure
BPM <- predict(bas.model, estimator = "BPM")
variable.names(BPM)
updist <- spavar$updist
fit <- glm(pa_glugea ~ avlength + O2_sat_av + Con_av^2 + meander + updist, data=environment2, family = binomial(link = "logit"))
predict <- ggpredict(fit, terms = "avlength")
ggpredict(fit,se=TRUE,interactive=TRUE,digits=3)
png(file="Glu_conductivity.png", res=600, width=3000, height=3000)
ggplot(predict, aes(x, predicted)) +
  theme_bw() +
  geom_line(color="red") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
  geom_point(data = environment2, aes(x=avlength, y=pa_glugea)) +
  labs(x=expression("Conductivity ["*mu*"S/cm]"), y=expression(italic(Trichodina)~"sp. [median infection intensity]"))+
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12))
dev.off()

#### Contracaecum ####
model <- glm(pa$Con ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + 
               log(NH4._av) + log(Nt_av) + log(SM_av) + pool_riffle + meander + spavar$netcen + 
               spavar$updist, data=environment2, family=binomial(link="logit"))
summary(model)
plot(model)   
step.model <- stepAIC(glm(pa$Con ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + 
                            log(NH4._av) + log(Nt_av) + log(SM_av) + pool_riffle + meander + spavar$netcen + 
                            spavar$updist, data=environment2, family=binomial(link="logit")),
                      direction = "both", 
                      trace = FALSE)
step.model

model <- glm(pa$Con ~ avlength + avcondition + O2_sat_av + Con_av^2 + log(NH4._av) + log(Nt_av)+ spavar$netcen, 
             data=environment2, family=binomial(link="logit"))
summary(model)

bas.model <- bas.lm(pa$Con ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + 
                      log(NH4._av) + log(Nt_av) + log(SM_av) + pool_riffle + meander + spavar$netcen + 
                      spavar$updist, data=environment2, prior="JZS")
summary(bas.model)
image(bas.model, rotate=F)
coef.model <- coef(bas.model)
abs(coef.model$postmean)-2*coef.model$postsd > 0
confint(coef.model)
plot(confint(coef.model, parm = 2:12))

confint <- confint(coef.model, parm = 2:13)
write.table(confint, 'Contra.txt', sep="\t")

pip <- summary(bas.model)
PIP[c(1:12),10] <- pip[2:13,1]*sign(coef.model$postmean[2:13])

#### Anguillicola ####
model <- glm(pa$Ang ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av)
             + log(Nt_av) + log(SM_av) + pool_riffle + meander + spavar$netcen + spavar$updist, 
             data=environment2, family=binomial(link="logit"))
summary(model)
#plot(model)   

# Model does not converge. Remove one predictor predictors.
model <- glm(pa$Ang ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + 
               log(Nt_av) + log(SM_av) + pool_riffle + meander, data=environment2, family=binomial(link="logit"))
summary(model)
#plot(model) 

stepAIC(glm(pa$Ang ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + 
              log(NH4._av) + log(Nt_av) + log(SM_av) + pool_riffle + meander + 
              spavar$netcen + spavar$updist, data=environment2, family=binomial(link="logit")),
        direction = "both", 
        trace = FALSE)

bas.model <- bas.glm(pa$Ang ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + 
                       log(NH4._av) + log(Nt_av) + log(SM_av) + pool_riffle + meander + 
                       spavar$netcen + spavar$updist, data=environment2, betaprior=g.prior(100), family=binomial)
summary(bas.model)
image(bas.model, rotate=F)
coef.model <- coef(bas.model)
abs(coef.model$postmean)-2*coef.model$postsd > 0
confint(coef.model)
plot(confint(coef.model, parm = 2:12))
confint <- confint(coef.model, parm = 2:13)
write.table(confint, 'Angui.txt', sep="\t")

pip <- summary(bas.model)
PIP[c(1:12),11] <- pip[2:13,1]*sign(coef.model$postmean[2:13])


# Effect of environment on Parasite Index
avPI <- aggregate(data$PI, by = list(data[,1]), function(x){mean(x, na.rm =T)})[,2]
avPI_ecto <- aggregate(data$PI_ecto, by = list(data[,1]), function(x){mean(x, na.rm =T)})[,2]
avPI_endo <- aggregate(data$PI_endo, by = list(data[,1]), function(x){mean(x, na.rm =T)})[,2]


# Bayesian
bas.model <- bas.lm(avPI ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av)
                    + log(Nt_av) + pool_riffle + meander + spavar$netcen + spavar$updist, 
                    data=environment2, prior="JZS")
summary(bas.model)
image(bas.model, rotate=F)
coef.model <- coef(bas.model)
abs(coef.model$postmean)-2*coef.model$postsd > 0
confint(coef.model)
plot(confint(coef.model, parm = 2:12))
confint <- confint(coef.model)
write.table(confint, 'avPI.txt', sep="\t")
pip <- summary(bas.model)
PIP <- cbind(PIP,pip[2:13,1]*sign(coef.model$postmean[2:13]))


# Bayesian
bas.model <- bas.lm(avPI_ecto ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + 
                      log(NH4._av) + log(Nt_av) + pool_riffle + meander + spavar$netcen + spavar$updist, 
                    data=environment2, prior="JZS")
summary(bas.model)
image(bas.model, rotate=F)
coef.model <- coef(bas.model)
abs(coef.model$postmean)-2*coef.model$postsd > 0
confint(coef.model)
plot(confint(coef.model, parm = 2:12))
confint <- confint(coef.model)
write.table(confint, 'avPI_ecto.txt', sep="\t")

pip <- summary(bas.model)
PIP <- cbind(PIP,pip[2:13,1]*sign(coef.model$postmean[2:13]))

#Figure
BPM <- predict(bas.model, estimator = "BPM")
variable.names(BPM)
updist <- spavar$updist
netcen <- spavar$netcen
fit <- lm(avPI_ecto ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(Nt_av) + pool_riffle + meander, data=environment2)
predict <- ggpredict(fit, terms = "COD_av")

png(file="Results/PIecto_COD.png", res=600, width=3000, height=3000)
ggplot(predict, aes(x, predicted)) +
  theme_bw() +
  geom_line(color="red") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1) +
  geom_point(data = environment2, aes(x=environment2$COD_av, y=avPI_ecto)) +
  labs(x="Chemical oxygen demand", y="Parasitation index (ectoparasites)")+
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12)) +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

# Bayesian
bas.model <- bas.lm(avPI_endo ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + 
                      log(NH4._av) + log(Nt_av) + pool_riffle + meander + spavar$netcen + spavar$updist,
                    data=environment2, prior="JZS")
summary(bas.model)
image(bas.model, rotate=F)
coef.model <- coef(bas.model)
abs(coef.model$postmean)-2*coef.model$postsd > 0
confint(coef.model)
plot(confint(coef.model, parm = 2:12))

confint <- confint(coef.model)
write.table(confint, 'avPI_endo.txt', sep="\t")

pip <- summary(bas.model)
PIP <- cbind(PIP,pip[2:13,1]*sign(coef.model$postmean[2:13]))

write.matrix(PIP, file="PIP.txt", sep="\t")
table <- as.data.frame(PIP)
write.table(table, file='table.txt', sep="\t")









#Structural Equation Model
library(piecewiseSEM)

