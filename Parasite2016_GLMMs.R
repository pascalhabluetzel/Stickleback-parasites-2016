# Load libraries
library(dplyr) # building data matrix
library(glmmTMB) # GLMMs
library(car) # ANOVA
library(performance) # for checking model requirements
library(buildmer) # evaluating and comparing different models (stepwise model selection)
library(lme4) # GLMMs

# Set working directory
setwd('C:/Users/pascalh/Documents/GitHub/Stickleback-parasites-2016')

# Load data
data_2016 <- read.csv("data_2016.csv", sep=';') #field and parasite data
env_av <- read.csv("Env_av.csv", sep=';') #environmental variables (average values)
env_max <- read.csv("env_max.csv", sep=';') #environmental variables (max. values)
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
cor(cbind(env_max, env_av, spavar))
# Nitrogen (max and average) is related to phosphorus (max and average) -> keep only average nitrogen
# max. and average oxygen are highly correlated -> keep only average
# max temperature and average temperature are not related to other parameters
# max and average pH are related to max and average conductivity -> keep only average conductivity
# speciesrichness is marginally related to network centrality (corr.coef. -0.60) -> keep only network centrality
# updist and updist2 are correlated -> keep only updist
# updist3 is not related to any other parameter -> keep

# make a new data frame by combinding parasite, environmental and space data
env_av_exp <- env_av %>% slice(rep(1:n(), table(as.factor(data_2016$site))))
env_max_exp <- env_max %>% slice(rep(1:n(), table(as.factor(data_2016$site))))
spavar_exp <- spavar %>% slice(rep(1:n(), table(as.factor(data_2016$site))))
spa.PCNM_exp <- spa.PCNM %>% slice(rep(1:n(), table(as.factor(data_2016$site))))

data <- cbind(data_2016, env_av_exp, env_max_exp, spavar_exp, spa.PCNM_exp)
data$site <- as.factor(data$site)
data$fish <- as.factor(data$fish)
data$length <- as.numeric(data$length)

# Effect of environment on host condition
confactor <- resid(lm(data$weight~data$length + data$Sex), na.action=na.exclude)
summary(confactor)
select <- buildglmmTMB(confactor ~ sqrt(T_av) + sqrt(Temperature) + sqrt(O2_av) + sqrt(con_av) + sqrt(KjN_av) + sqrt(netcen) + sqrt(updist) + (1|site),
                              data=data,
                              crit="AIC")
select
model <- lme(confactor ~ sqrt(T_av) + sqrt(Temperature) + sqrt(netcen), random=~1|site, data=data, na.action=na.omit)
summary(model)

# Effect of environment on Parasite Index

#calculate parasite infection index
PI_ecto <- 1:nrow(data)
for(j in 1:nrow(data)){
  PI_ecto[j] <- max(data$gyro[j]/sd(data$gyro, na.rm=T))*(data$gyro[j]/sd(data$gyro, na.rm=T)) + max(data$tricho[j]/sd(data$tricho, na.rm=T))*(data$tricho[j]/sd(data$tricho, na.rm=T)) + max(data$glugea[j]/sd(data$glugea, na.rm=T))*(data$glugea[j]/sd(data$glugea, na.rm=T))
}

PI_endo <- 1:nrow(data)
for(j in 1:nrow(data)){
  PI_endo[j] <- 
    max(data$contra[j]/sd(data$contra, na.rm=T))*(data$contra[j]/sd(data$contra, na.rm=T))
  + max(data$cystsliver[j]/sd(data$cystsliver, na.rm=T))*(data$cystsliver[j]/sd(data$cystsliver, na.rm=T))
  + max(data$proteo[j]/sd(data$proteo, na.rm=T))*(data$proteo[j]/sd(data$proteo, na.rm=T))
  + max(data$acantho[j]/sd(data$acantho, na.rm=T))*(data$acantho[j]/sd(data$acantho, na.rm=T))
  + max(data$cama[j]/sd(data$cama, na.rm=T))*(data$cama[j]/sd(data$cama, na.rm=T))
  + max(data$angui[j]/sd(data$angui, na.rm=T))*(data$angui[j]/sd(data$angui, na.rm=T))
  + max(data$cistsintestine[j]/sd(data$cistsintestine, na.rm=T))*(data$cistsintestine[j]/sd(data$cistsintestine, na.rm=T))
}

PI <- 1:nrow(data)
for(j in 1:nrow(data)){
  PI[j] <- 
    max(data$gyro[j]/sd(data$gyro, na.rm=T))*(data$gyro[j]/sd(data$gyro, na.rm=T))
  + max(data$tricho[j]/sd(data$tricho, na.rm=T))*(data$tricho[j]/sd(data$tricho, na.rm=T))
  + max(data$glugea[j]/sd(data$glugea, na.rm=T))*(data$glugea[j]/sd(data$glugea, na.rm=T))
  + max(data$contra[j]/sd(data$contra, na.rm=T))*(data$contra[j]/sd(data$contra, na.rm=T))
  + max(data$cystsliver[j]/sd(data$cystsliver, na.rm=T))*(data$cystsliver[j]/sd(data$cystsliver, na.rm=T))
  + max(data$proteo[j]/sd(data$proteo, na.rm=T))*(data$proteo[j]/sd(data$proteo, na.rm=T))
  + max(data$acantho[j]/sd(data$acantho, na.rm=T))*(data$acantho[j]/sd(data$acantho, na.rm=T))
  + max(data$cama[j]/sd(data$cama, na.rm=T))*(data$cama[j]/sd(data$cama, na.rm=T))
  + max(data$angui[j]/sd(data$angui, na.rm=T))*(data$angui[j]/sd(data$angui, na.rm=T))
  + max(data$cistsintestine[j]/sd(data$cistsintestine, na.rm=T))*(data$cistsintestine[j]/sd(data$cistsintestine, na.rm=T))
}

select <- buildglmmTMB(PI ~ Sex + sqrt(T_av) + confactor + sqrt(Temperature) + sqrt(O2_av) + sqrt(con_av) + sqrt(KjN_av) + sqrt(netcen) + sqrt(updist) + (1|site),
                       data=data,
                       crit="AIC")
select
model <- lme(PI ~ sqrt(T_av) + confactor, random=~1|site, data=data, na.action=na.omit)
summary(model)

select <- buildglmmTMB(PI_ecto ~ Sex + sqrt(T_av) + confactor + sqrt(Temperature) + sqrt(O2_av) + sqrt(con_av) + sqrt(KjN_av) + sqrt(netcen) + sqrt(updist) + (1|site) + (1|ecto_screener),
                       data=data,
                       crit="AIC")
select
model <- lme(PI_ecto ~ confactor + sqrt(con_av), random=~1|site, data=data, na.action=na.omit)
summary(model)

select <- buildglmmTMB(PI_endo ~ Sex + sqrt(T_av) + confactor + sqrt(Temperature) + sqrt(O2_av) + sqrt(con_av) + sqrt(KjN_av) + sqrt(netcen) + sqrt(updist) + (1|site) + (1|endo_screener),
                       data=data,
                       crit="AIC")
select
model <- lme(PI_endo ~ confactor + sqrt(Temperature) + sqrt(netcen), random=~1|site, data=data, na.action=na.omit)
summary(model)

# Effect of environment on individual parasite taxa

#### 1.2 ABUNDANCE ####

# ZIGLMM (glmmTMB package)

#check whether spatial variables are correlated with environmental variables
cor(cbind(data$T_av, data$Temperature, data$O2_av, data$con_av, data$KjN_av, data$netcen, data$updist, data$updist3, data$speciesrichness))

#Gyrodactylus
fit_zipoisson <- buildglmmTMB(gyro ~ Sex + sqrt(length) + confactor + sqrt(T_av) + sqrt(Temperature) + sqrt(O2_av) + sqrt(con_av) + sqrt(KjN_av) + sqrt(netcen) + sqrt(updist) + (1|site) + (1|ecto_screener),
                              data=data,
                              ziformula= ~ Sex + sqrt(length) + sqrt(T_av) + sqrt(Temperature) + sqrt(O2_av) + sqrt(con_av) + sqrt(KjN_av) + sqrt(netcen) + sqrt(updist),
                              family=poisson,
                              crit="AIC")
fit_zipoisson

fit_zipoisson <- glmmTMB(gyro ~ sqrt(Temperature) + (1|site),
                         data=data,
                         ziformula= ~ Sex + sqrt(length) + confactor + sqrt(T_av) + sqrt(Temperature) + sqrt(O2_av) + sqrt(con_av) + sqrt(KjN_av) + sqrt(netcen) + sqrt(updist),
                         family=poisson)
summary(fit_zipoisson)
Anova(fit_zipoisson)

X.cond = model.matrix(lme4::nobars(formula(fit_zipoisson)[-2]), data)
beta.cond = fixef(fit_zipoisson)$cond
pred.cond = X.cond %*% beta.cond
ziformula = fit_zipoisson$modelInfo$allForm$ziformula
X.zi = model.matrix(lme4::nobars(ziformula), data)
beta.zi =fixef(fit_zipoisson)$zi
pred.zi = X.zi %*% beta.zi
pred.ucount =exp(pred.cond)*(1-plogis(pred.zi))

collinearity <- check_collinearity(fit_zipoisson)
plot(collinearity)
normality <- check_normality(fit_zipoisson, effects="random")


#Trichodina (netcen and updist are not in the model because the model did not converge)
fit_zipoisson <- buildglmmTMB(tricho ~ Sex + sqrt(length) + confactor + sqrt(T_av) + sqrt(Temperature) + sqrt(O2_av) + sqrt(con_av) + sqrt(KjN_av) + (1|site) + (1|ecto_screener),
                              data=data,
                              ziformula= ~ Sex + sqrt(length) + sqrt(T_av) + sqrt(Temperature) + sqrt(O2_av) + sqrt(con_av) + sqrt(KjN_av),
                              family=poisson,
                              crit="AIC")
fit_zipoisson

fit_zipoisson <- glmmTMB(tricho ~ Sex + sqrt(length) + confactor + sqrt(con_av) + sqrt(KjN_av) + (1|site),
                         data=data,
                         ziformula= ~ Sex + sqrt(length) + confactor + sqrt(T_av) + sqrt(Temperature) + sqrt(O2_av) + sqrt(con_av) + sqrt(KjN_av),
                         family=poisson)
summary(fit_zipoisson)
Anova(fit_zipoisson)

collinearity <- check_collinearity(fit_zipoisson)
plot(collinearity)
normality <- check_normality(fit_zipoisson, effects="random")


#Glugea
fit_zipoisson <- buildglmmTMB(glugea ~ Sex + sqrt(length) + confactor + sqrt(T_av) + sqrt(Temperature) + sqrt(O2_av) + sqrt(con_av) + sqrt(KjN_av) + sqrt(netcen) + sqrt(updist) + (1|site) + (1|ecto_screener),
                              data=data,
                              ziformula= ~ Sex + sqrt(length) + confactor + sqrt(T_av) + sqrt(Temperature) + sqrt(O2_av) + sqrt(con_av) + sqrt(KjN_av) + sqrt(netcen) + sqrt(updist),
                              family=poisson,
                              crit="AIC")
fit_zipoisson

fit_zipoisson <- glmmTMB(glugea ~ sqrt(length) + sqrt(T_av) + confactor + sqrt(Temperature) + sqrt(con_av) + sqrt(KjN_av) + sqrt(netcen) + sqrt(updist) + (1|site),
                         data=data,
                         ziformula= ~ Sex + sqrt(length) + confactor + sqrt(T_av) + sqrt(Temperature) + sqrt(O2_av) + sqrt(con_av) + sqrt(KjN_av) + sqrt(netcen) + sqrt(updist),
                         family=poisson)
summary(fit_zipoisson)
Anova(fit_zipoisson)

collinearity <- check_collinearity(fit_zipoisson)
plot(collinearity)
normality <- check_normality(fit_zipoisson, effects="random")
plot(normality)


#Contracaecum
fit_zipoisson <- buildglmmTMB(contra ~ Sex + sqrt(length) + confactor + sqrt(T_av) + sqrt(Temperature) + sqrt(O2_av) + sqrt(con_av) + sqrt(KjN_av) + sqrt(netcen) + sqrt(updist) + (1|site) + (1|endo_screener),
                              data=data,
                              ziformula= ~ Sex + sqrt(length) + confactor + sqrt(T_av) + sqrt(Temperature) + sqrt(O2_av) + sqrt(con_av) + sqrt(KjN_av) + sqrt(netcen) + sqrt(updist),
                              family=poisson,
                              crit="AIC")
fit_zipoisson

fit_zipoisson <- glmmTMB(contra ~ sqrt(length) + confactor + sqrt(T_av) + sqrt(netcen),
                         data=data,
                         ziformula= ~ Sex + sqrt(length) + confactor + sqrt(T_av) + sqrt(Temperature) + sqrt(O2_av) + sqrt(con_av) + sqrt(KjN_av) + sqrt(netcen) + sqrt(updist),
                         family=poisson)
summary(fit_zipoisson)# AIC 291.3
fit_zipoisson <- glmmTMB(contra ~ sqrt(length) + sqrt(T_av) + sqrt(netcen) + (1|site),
                         data=data,
                         ziformula= ~ Sex + sqrt(length) + sqrt(T_av) + sqrt(Temperature) + sqrt(O2_av) + sqrt(con_av) + sqrt(KjN_av) + sqrt(netcen) + sqrt(updist),
                         family=poisson)
summary(fit_zipoisson)# AIC 292.1
Anova(fit_zipoisson)

#Anguillicola crassus
fit_zipoisson <- buildglmmTMB(angui ~ Sex + sqrt(length) + confactor + sqrt(T_av) + sqrt(Temperature) + sqrt(O2_av) + sqrt(con_av) + sqrt(KjN_av) + sqrt(netcen) + sqrt(updist) + (1|site) + (1|endo_screener),
                              data=data,
                              ziformula= ~ Sex + sqrt(length) + confactor + sqrt(T_av) + sqrt(Temperature) + sqrt(O2_av) + sqrt(con_av) + sqrt(KjN_av) + sqrt(netcen) + sqrt(updist),
                              family=poisson,
                              crit="AIC")
fit_zipoisson

fit_zipoisson <- glmmTMB(angui ~ Sex + sqrt(length) + sqrt(T_av) + sqrt(O2_av) + sqrt(con_av) + sqrt(KjN_av) + sqrt(updist),
                         data=data,
                         ziformula= ~ Sex + sqrt(length) + confactor + sqrt(T_av) + sqrt(Temperature) + sqrt(O2_av) + sqrt(con_av) + sqrt(KjN_av) + sqrt(netcen) + sqrt(updist),
                         family=poisson)
summary(fit_zipoisson)# AIC 441.4
fit_zipoisson <- glmmTMB(angui ~ sqrt(length) + sqrt(KjN_av) + (1|site),
                         data=data,
                         ziformula= ~ Sex + sqrt(length) + confactor + sqrt(T_av) + sqrt(Temperature) + sqrt(O2_av) + sqrt(con_av) + sqrt(KjN_av) + sqrt(netcen) + sqrt(updist),
                         family=poisson)
summary(fit_zipoisson)# AIC 438.0
Anova(fit_zipoisson)



