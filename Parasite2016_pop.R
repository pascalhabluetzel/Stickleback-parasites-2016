library(ggplot2)
library(readr)
library(ggplot2)
library(gridExtra)
library(plyr)
library(vegan)
library(afex)
library(MASS)
library(effects)
library(lme4)
library(car)
library(corrplot)
library(dplyr) # building data matrix
library(glmmTMB) # GLMMs
library(car) # ANOVA
library(performance) # for checking model requirements
library(buildmer) # evaluating and comparing different models (stepwise model selection)
library(lme4) # GLMMs
library(vegan) # for MEMs
library(nlme) # for lme

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
avin = aggregate(data[,c(23:25,27:33)], by = list(data[,1]), function(x){mean(x[x >0], na.rm = T)}) 
avab = aggregate(data[,c(23:25,27:33)], by = list(data[,1]), function(x){mean(x, na.rm =T)})
prev = aggregate(data[,c(23:25,27:33)], by = list(data[,1]), function(x){sum(x >0, na.rm = T)/length(x)})
medin = aggregate(data[,c(23:25,27:33)], by = list(data[,1]), function(x){median(x[x >0], na.rm = T)}) 

avlength <- aggregate(data$length, by = list(data[,1]), function(x){mean(x, na.rm =T)})[,2]
avcondition <- aggregate(condition, by = list(data[,1]), function(x){mean(x, na.rm =T)})[,2]

avin[is.na(avin)] <- 0
avab[is.na(avab)] <- 0
prev[is.na(prev)] <- 0
medin[is.na(medin)] <- 0

# General structure of the tests for individual parasites:
# 1. Models for environmental effect on host condition and length
# 2. Four measures at populatin level: average abundance, prevalence and median infection intensity
# 3. Three parasite indices: PI all parasites, PI for ectoparasites and PI for endoparasites 
# 4. Five parasite taxa: Gyrodactylus, Trichodina, Glugea, Contracaecum, Anguillicola
# 5. For each response variable:
#    6. Plot of raw data (value vs. site)
#    7. Model with all explanatory variables
#    8. Check for testing requirements
#    9. Model selection
#    10. Check for testing requirements
#    11. Plot of residuals

#### Condition ####
plot(avcondition)

model <- lm(avcondition ~ T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + spavar$netcen + spavar$updist, data=env)
summary(model)
#plot(model)                      

step.model <- stepAIC(lm(avcondition ~ T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + spavar$netcen + spavar$updist, data=env),
                      direction = "both", 
                      trace = FALSE)
step.model

summary(model <- lm(avcondition ~ con_av + NO2_av + spavar$netcen, data=env))
#plot(model)

res_con_av <- resid(lm(avcondition ~ NO2_av + spavar$netcen, data=env))
plot(res_con_av ~ con_av, data=env)
abline(lm(res_con_av ~ con_av, data=env))


#### Length ####
plot(avlength)

model <- lm(avlength ~ T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + spavar$netcen + spavar$updist, data=env)
summary(model)
#plot(model)                      

step.model <- stepAIC(lm(avlength ~ T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + spavar$netcen + spavar$updist, data=env),
                      direction = "both", 
                      trace = FALSE)
step.model

# Check the effect of the outlier (location 12)
step.model <- stepAIC(lm(avlength[-c(10)] ~ T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + spavar$netcen[-c(10)] + spavar$updist[-c(10)], data=env[-c(10),]),
                      direction = "both", 
                      trace = FALSE)
step.model

# When location 12 is excluded, the effect of Cl disappears. Do not use Cl in the final model.

summary(model <- lm(avlength ~ con_av + NO3_av + spavar$netcen, data=env))
#plot(model)

res_netcen <- resid(lm(avlength ~ con_av + NO3_av, data=env))
plot(res_netcen ~ spavar$netcen)
abline(lm(res_netcen ~ spavar$netcen))


#### Gyrodactylus ####
#### Average infection intensity ####
plot(avin$gyro)

model <- lm(avin$gyro ~ avlength + avcondition + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + spavar$netcen + spavar$updist, data=env)
summary(model)
#plot(model)                      

step.model <- stepAIC(lm(avin$gyro ~ avlength + avcondition + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + spavar$netcen + spavar$updist, data=env),
                      direction = "both", 
                      trace = FALSE)
step.model

summary(model <- lm(avin$gyro ~ avlength + T_av + con_av + NO3_av, data=env))
#plot(model)

res_avlength <- resid(lm(avin$gyro ~ T_av + con_av + NO3_av, data=env))
plot(res_avlength ~ avlength, data=env)
abline(lm(res_avlength ~ avlength, data=env))

res_T_av <- resid(lm(avin$gyro ~ avlength + con_av + NO3_av, data=env))
plot(res_T_av ~ T_av, data=env)
abline(lm(res_T_av ~ T_av, data=env))

res_con_av <- resid(lm(avin$gyro ~ avlength + T_av + NO3_av, data=env))
plot(res_con_av ~ con_av, data=env)
abline(lm(res_con_av ~ con_av, data=env))

res_NO3_av <- resid(lm(avin$gyro ~ avlength + T_av + con_av, data=env))
plot(res_NO3_av ~ NO3_av, data=env)
abline(lm(res_NO3_av ~ NO3_av, data=env))


#### Gyrodactylus ####
#### Average abundance ####
plot(avab$gyro)

model <- lm(avab$gyro ~ avlength + avcondition + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + spavar$netcen + spavar$updist, data=env)
summary(model)
#plot(model)                      

step.model <- stepAIC(lm(avab$gyro ~ avlength + avcondition + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + spavar$netcen + spavar$updist, data=env),
                      direction = "both", 
                      trace = FALSE)
step.model

# Check the effect of the outlier (location 12)
step.model <- stepAIC(lm(avab$gyro[-c(10)] ~ avlength[-c(10)] + avcondition[-c(10)] + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + spavar$netcen[-c(10)] + spavar$updist[-c(10)], data=env[-c(10),]),
                      direction = "both", 
                      trace = FALSE)
step.model

# The exclusion of location 12 does change the outcome significantly. Use only the model that was suggested upon the exclusion of location 12.
summary(model <- lm(avab$gyro ~ avlength, data=env))
#plot(model)


#### Gyrodactylus ####
#### Prevalence ####
plot(prev$gyro)

model <- lm(prev$gyro ~ avlength + avcondition + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + spavar$netcen + spavar$updist, data=env)
summary(model)
#plot(model)                      

step.model <- stepAIC(lm(prev$gyro ~ avlength + avcondition + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + spavar$netcen + spavar$updist, data=env),
                      direction = "both", 
                      trace = FALSE)
step.model

# Check the effect of the exclusion of location 12
step.model <- stepAIC(lm(prev$gyro[-c(10)] ~ avlength[-c(10)] + avcondition[-c(10)] + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + spavar$netcen[-c(10)] + spavar$updist[-c(10)], data=env[-c(10),]),
                      direction = "both", 
                      trace = FALSE)
step.model


summary(model <- lm(prev$gyro ~ avcondition + con_av + COD_av, data=env))
#plot(model)

res_avcondition <- resid(lm(prev$gyro ~ con_av + COD_av, data=env))
plot(res_avcondition ~ avcondition, data=env)
abline(lm(res_avcondition ~ avcondition, data=env))

res_con_av <- resid(lm(prev$gyro ~ avcondition + COD_av, data=env))
plot(res_con_av ~ con_av, data=env)
abline(lm(res_con_av ~ con_av, data=env))

res_COD_av <- resid(lm(prev$gyro ~ avcondition + con_av, data=env))
plot(res_COD_av ~ COD_av, data=env)
abline(lm(res_COD_av ~ COD_av, data=env))

#### Gyrodactylus ####
#### Median infection intensity ####
plot(medin$gyro)

model <- lm(medin$gyro ~ avlength + avcondition + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + spavar$netcen + spavar$updist, data=env)
summary(model)
#plot(model)                      

step.model <- stepAIC(lm(medin$gyro ~ avlength + avcondition + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + spavar$netcen + spavar$updist, data=env),
                      direction = "both", 
                      trace = FALSE)
step.model

# Check the effect of the exclusion of location 12
step.model <- stepAIC(lm(medin$gyro[-c(10)] ~ avlength[-c(10)] + avcondition[-c(10)] + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + spavar$netcen[-c(10)] + spavar$updist[-c(10)], data=env[-c(10),]),
                      direction = "both", 
                      trace = FALSE)
step.model

summary(model <- lm(medin$gyro ~ T_av, data=env))
#plot(model)


#### Trichodina ####
#### Average abundance ####
plot(avab$tricho)

model <- lm(avab$tricho ~ avlength + avcondition + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + spavar$netcen + spavar$updist, data=env)
summary(model)
#plot(model)                      

step.model <- stepAIC(lm(avab$tricho ~ avlength + avcondition + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + spavar$netcen + spavar$updist, data=env),
                      direction = "both", 
                      trace = FALSE)
step.model

# Check the effect of location 12
step.model <- stepAIC(lm(avab$tricho[-10] ~ avlength[-10] + avcondition[-10] + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + spavar$netcen[-10] + spavar$updist[-10], data=env[-10,]),
                      direction = "both", 
                      trace = FALSE)
step.model

summary(model <- lm(avab$tricho ~ con_av + COD_av, data=env))
#plot(model)
summary(model <- lm(avab$tricho[-10] ~ con_av + COD_av, data=env[-10,]))

#### Trichodina ####
#### Prevalence ####
plot(prev$tricho)

model <- lm(prev$tricho ~ avlength + avcondition + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + spavar$netcen + spavar$updist, data=env)
summary(model)
#plot(model)                      

step.model <- stepAIC(lm(prev$tricho ~ avlength + avcondition + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + spavar$netcen + spavar$updist, data=env),
                      direction = "both", 
                      trace = FALSE)
step.model

# Check the effect of removal of location 12
step.model <- stepAIC(lm(prev$tricho[-10] ~ avlength[-10] + avcondition[-10] + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + spavar$netcen[-10] + spavar$updist[-10], data=env[-10,]),
                      direction = "both", 
                      trace = FALSE)
step.model

# No difference, continue with all sites

summary(model <- lm(prev$tricho ~ T_av + con_av + O2_sat_av + COD_av + NH4_av + spavar$updist, data=env))
#plot(model)

res_con_av <- resid(lm(prev$tricho ~ T_av + O2_sat_av + COD_av + NH4_av + spavar$updist, data=env))
plot(res_con_av ~ con_av, data=env)
abline(lm(res_con_av ~ con_av, data=env))

res_COD_av <- resid(lm(prev$tricho ~ T_av + con_av + O2_sat_av + NH4_av + spavar$updist, data=env))
plot(res_COD_av ~ COD_av, data=env)
abline(lm(res_COD_av ~ COD_av, data=env))

res_NH4_av <- resid(lm(prev$tricho ~ T_av + con_av + O2_sat_av + COD_av + spavar$updist, data=env))
plot(res_NH4_av ~ NH4_av, data=env)
abline(lm(res_NH4_av ~ NH4_av, data=env))

res_updist <- resid(lm(prev$tricho ~ T_av + con_av + O2_sat_av + COD_av + NH4_av, data=env))
plot(res_updist ~ spavar$updist, data=env)
abline(lm(res_updist ~ spavar$updist, data=env))

#### Trichodina ####
#### Median infection intensity ####
plot(medin$trich)

model <- lm(medin$gyro ~ avlength + avcondition + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + spavar$netcen + spavar$updist, data=env)
summary(model)
#plot(model)                      

step.model <- stepAIC(lm(medin$tricho ~ avlength + avcondition + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + spavar$netcen + spavar$updist, data=env),
                      direction = "both", 
                      trace = FALSE)
step.model

# Check the effect of removal of location 12
step.model <- stepAIC(lm(medin$tricho[-10] ~ avlength[-10] + avcondition[-10] + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + spavar$netcen[-10] + spavar$updist[-10], data=env[-10,]),
                      direction = "both", 
                      trace = FALSE)
step.model

# Use only variables from the model selection without location 12

summary(model <- lm(medin$tricho ~ con_av + COD_av + NH4_av, data=env))
#plot(model)

summary(model <- lm(medin$tricho[-10] ~ con_av + COD_av + NH4_av, data=env[-10,]))
#plot(model)

# Only effect of conductivity remains significant after removal of location 12

res_con_av <- resid(lm(medin$tricho ~ COD_av + NH4_av, data=env))
plot(res_con_av ~ con_av, data=env)
abline(lm(res_con_av ~ con_av, data=env))


#### Glugea ####
plot(avab$glugea)

glugea_pa <- ifelse(avab$glugea>0, 1, 0)

model <- glm(glugea_pa ~ avlength + avcondition + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + spavar$netcen + spavar$updist, data=env, family="binomial")
summary(model)
#plot(model)                      

step.model <- stepAIC(glm(glugea_pa ~ avlength + avcondition + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + spavar$netcen + spavar$updist, data=env, family="binomial"),
                      direction = "both", 
                      trace = FALSE)
step.model

model <- glm(glugea_pa ~ con_av + spavar$updist, data=env, family="binomial")
summary(model)

# plotting residual effects makes little sense with binomial data -> use model predictions
plot(as.factor(glugea_pa), env$con_av)
plot(as.factor(glugea_pa), spavar$updist)


#### Contracaecum ####
plot(avab$contra)

contra_pa <- ifelse(avab$contra>0, 1, 0)

model <- glm(contra_pa ~ avlength + avcondition + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + spavar$netcen + spavar$updist, data=env, family="binomial")
summary(model)
#plot(model)                      

step.model <- stepAIC(glm(contra_pa ~ avlength + avcondition + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + spavar$netcen + spavar$updist, data=env, family="binomial"),
                      direction = "both", 
                      trace = FALSE)
step.model

model <- glm(contra_pa ~ avlength + O2_sat_av + Cl_av + NO3_av + NO2_av, data=env, family="binomial")
summary(model)

# plotting residual effects makes little sense with binomial data -> use model predictions
plot(as.factor(contra_pa), avlength)
plot(as.factor(contra_pa), env$O2_sat_av)


#### Anguillicola ####
plot(avab$angui)

angui_pa <- ifelse(avab$angui>0, 1, 0)

model <- glm(angui_pa ~ avlength + avcondition + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + spavar$netcen + spavar$updist, data=env, family="binomial")
summary(model)
#plot(model)                      

step.model <- stepAIC(glm(angui_pa ~ avlength + avcondition + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + spavar$netcen + spavar$updist, data=env, family="binomial"),
                      direction = "both", 
                      trace = FALSE)
step.model

model <- glm(angui_pa ~ avcondition + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av, data=env, family="binomial")
summary(model)

# plotting residual effects makes little sense with binomial data -> use model predictions
plot(as.factor(angui_pa), env$T_av)



library(piecewiseSEM)

SEM <- as.data.frame(cbind(prev$gyro, avcondition, avlength, env$T_av, env$con_av, env$COD_av, env$NO3_av))
colnames(SEM) <- c("gyro", "avcondition", "avlength", "T_av", "con_av", "COD_av", "NO3_av")
model <- psem(
  lm(gyro ~ avcondition + avlength + T_av + con_av + COD_av + NO3_av, data=SEM),
  lm(avcondition ~ con_av, data=SEM)
)
summary(model)

model <- psem(
  lm(gyro ~ avcondition + avlength + T_av + con_av + COD_av + NO3_av, data=SEM),
  lm(avcondition ~ con_av, data=SEM)
)
summary(model)
