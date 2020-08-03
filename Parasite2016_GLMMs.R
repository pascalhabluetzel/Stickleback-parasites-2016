# Load libraries
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
cor(spavar$updist, spa.PCNM$V2)
plot(spavar$updist, spa.PCNM$V2) #second MEM also corresponds to distance from Demer-Dijle confluence
cor(spavar$netcen, spa.PCNM$V2)
plot(spavar$netcen, spa.PCNM$V2) #second MEM corresponds to network centrality
plot(spavar$updist, spavar$netcen)
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
# temp: measures somewhat correlated -> keep only T_av
# O2 and O2_sat: all measures somewhat correlated -> keep only O2_sat_av
# speciesrichness is marginally related to network centrality (corr.coef. -0.60) -> keep only network centrality
# updist and updist2 are correlated -> keep only updist
# updist3 is not related to any other parameter -> keep?

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

# Effect of environment (average) on host condition
condition <- resid(lm(data$weight~data$length + data$Sex), na.action=na.exclude)
summary(condition)
select <- buildglmmTMB(condition ~ T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + netcen + updist,
                              data=data,
                              include = ~ (1|site),
                              direction = c('order', 'forward'),
                              crit="LRT")
select

model <- lme(condition ~ con_av, random=~1|site, data=data, na.action=na.omit)
summary(model)

# Effect of environment (maximum) on host condition
condition <- resid(lm(data$weight~data$length + data$Sex), na.action=na.exclude)
summary(condition)
select <- buildglmmTMB(condition ~ T_max + con_max + O2_sat_max + Cl_max + COD_max + NH4_max + NO3_max + NO2_max + netcen + updist,
                       data=data,
                       include = ~ (1|site),
                       direction = c('order', 'forward'),
                       crit="LRT")
select
model <- lme(condition ~ T_max, random=~1|site, data=data, na.action=na.omit)
summary(model)

res_T_max <- resid(lme(condition ~ 1, random=~1|site, data=data, na.action=na.omit))
plot(res_T_max ~ data$T_max)
abline(lm(res_T_max ~ data$T_max))

# Effect of environment (point value) on host condition
condition <- resid(lm(data$weight~data$length + data$Sex), na.action=na.exclude)
summary(condition)
select <- buildglmmTMB(condition ~ T + con + O2_sat + Cl + COD + NH4 + NO3 + NO2 + netcen + updist,
                       data=data,
                       include = ~ (1|site),
                       direction = c('order', 'forward'),
                       crit="LRT")
select
model <- lme(condition ~ NO3, random=~1|site, data=data, na.action=na.omit)
summary(model)

res_NO3 <- resid(lme(condition ~ 1, random=~1|site, data=data, na.action=na.omit))
plot(res_NO3 ~ data$NO3)
abline(lm(res_NO3 ~ data$NO3))



# Effect of environment (average) on host condition
condition <- resid(lm(data$weight~data$length + data$Sex), na.action=na.exclude)
summary(condition)
select <- buildglmmTMB(condition ~ T + con + O2_sat + Cl + COD + NH4 + NO3 + NO2 + netcen + updist + site,
                       data=data,
                       include = ~ site,
                       direction = c('order', 'forward'),
                       crit="LRT")
select
model <- lme(condition ~ 1, random=~1|site, data=data, na.action=na.omit)
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

select <- buildglmmTMB(PI ~ Sex + sqrt(length) + condition + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + netcen + updist,
                       data=data,
                       include = ~ (1|site),
                       direction = c('order', 'forward'),
                       crit="LRT")
select
model <- lme(PI ~ sqrt(length) + Cl_av + T_av, random=~1|site, data=data, na.action=na.omit)
summary(model)

collinearity <- check_collinearity(model)
plot(collinearity)

res_length <- resid(lme(PI ~ Cl_av + T_av, random=~1|site, data=data, na.action=na.omit))
plot(res_length ~ data$length)
abline(lm(res_length ~ data$length))

res_Cl <- resid(lme(PI ~ length + T_av, random=~1|site, data=data, na.action=na.omit))
plot(res_Cl ~ data$Cl_av)
abline(lm(res_Cl ~ data$Cl_av))

res_T <- resid(lme(PI ~ sqrt(length) + Cl_av, random=~1|site, data=data, na.action=na.omit))
plot(res_T ~ data$T_av)
abline(lm(res_T ~ data$T_av))

collinearity <- check_collinearity(model)
plot(collinearity)



select <- buildglmmTMB(PI_ecto ~ Sex + sqrt(length) + condition + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + netcen + updist + (1|ecto_screener),
                       data=data,
                       include = ~ (1|site),
                       direction = c('order', 'forward'),
                       crit="LRT")
select
model <- lme(PI_ecto ~ sqrt(length) + Cl_av + con_av + NO3_av, random=~1|site, data=data, na.action=na.omit)
summary(model)

collinearity <- check_collinearity(model)
plot(collinearity)

select <- buildglmmTMB(PI_endo ~ Sex + sqrt(length) + confactor + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + netcen + updist + (1|endo_screener),
                       data=data,
                       include = ~ (1|site),
                       direction = c('order', 'forward'),
                       crit="AIC")
select
model <- lme(PI_endo ~ sqrt(length) + T_av + NH4_av + NO3_av, random=~1|site, data=data, na.action=na.omit)
summary(model)

collinearity <- check_collinearity(model)
plot(collinearity)

# Effect of environment on individual parasite taxa

#### 1.2 ABUNDANCE ####

# ZIGLMM (glmmTMB package)

#Gyrodactylus
fit_zipoisson <- buildglmmTMB(gyro ~ Sex + sqrt(length) + condition + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + netcen + updist + (1|ecto_screener),
                              data=data,
                              ziformula= ~ Sex + sqrt(length) + confactor + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + netcen + updist,
                              family=poisson,
                              include = ~ (1|site),
                              crit="LRT")
fit_zipoisson

fit_zipoisson <- buildglmmTMB(gyro ~ Sex + sqrt(length) + condition + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + netcen + updist + (1|ecto_screener),
                              data=data,
                              ziformula= ~ 1 + (1|site),
                              family=poisson,
                              include = ~ (1|site),
                              direction = c('order', 'forward'),
                              crit="LRT")
fit_zipoisson

fit_zipoisson <- glmmTMB(gyro ~ 1 + (1|site) + (1|ecto_screener),
                         data=data,
                         ziformula= ~ 1 + (1|site),
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
fit_zipoisson <- buildglmmTMB(tricho ~ Sex + length + condition + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + netcen + updist + (1|ecto_screener),
                              data=data,
                              ziformula= ~ 1 + (1|site),
                              family=poisson,
                              include = ~ (1|site),
                              direction = c('order', 'forward'),
                              crit="LRT")
fit_zipoisson

fit_zipoisson <- glmmTMB(tricho ~ Sex + length + condition + Cl_av + (1|ecto_screener) + (1|site),
                         data=data,
                         ziformula= ~ 1 + (1|site),
                         family=poisson)
summary(fit_zipoisson)
Anova(fit_zipoisson)

collinearity <- check_collinearity(fit_zipoisson)
plot(collinearity)
normality <- check_normality(fit_zipoisson, effects="random")

model <- glmmTMB(tricho ~ condition + Sex + Cl_av + (1|ecto_screener) + (1|site),
                         data=data,
                         ziformula= ~ 1 + (1|site),
                         family=poisson)
res_length <- resid(model)
plot(res_length ~ data$length)
abline(lm(res_length ~ data$length))

model <- glmmTMB(tricho ~ sqrt(length) + Cl_av + (1|ecto_screener) + (1|site),
                 data=data,
                 ziformula= ~ 1 + (1|site),
                 family=poisson)
res_sex <- resid(model)
boxplot(res_sex ~ data$Sex)


#Glugea
fit_zipoisson <- buildglmmTMB(glugea ~ Sex + sqrt(length) + confactor + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + netcen + updist + (1|ecto_screener),
                              data=data,
                              ziformula= ~ 1 + (1|site),
                              family=poisson,
                              include = ~ (1|site),
                              direction = c('order', 'forward'),
                              crit="LRT")
fit_zipoisson

fit_zipoisson <- glmmTMB(glugea ~ con_av + (1|site),
                         data=data,
                         ziformula= ~ 1 + (1|site),
                         family=poisson)
summary(fit_zipoisson)
Anova(fit_zipoisson)

collinearity <- check_collinearity(fit_zipoisson)
plot(collinearity)
normality <- check_normality(fit_zipoisson, effects="random")
plot(normality)


#Contracaecum
fit_zipoisson <- buildglmmTMB(contra ~ Sex + sqrt(length) + confactor + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + netcen + updist + (1|endo_screener),
                              data=data,
                              ziformula= ~ 1 + (1|site),
                              family=poisson,
                              include = ~ (1|site),
                              direction = c('order', 'forward'),
                              crit="LRT")
fit_zipoisson

fit_zipoisson <- glmmTMB(contra ~ 1 + Sex + (1|site),
                         data=data,
                         ziformula= ~ 1 + (1|site),
                         family=poisson)
summary(fit_zipoisson)
Anova(fit_zipoisson)

#Anguillicola crassus
fit_zipoisson <- buildglmmTMB(angui ~ Sex + sqrt(length) + confactor + T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + netcen + updist + (1|endo_screener),
                              data=data,
                              ziformula= ~ 1 + (1|site),
                              family=poisson,
                              include = ~ (1|site),
                              direction = c('order', 'forward'),
                              crit="LRT")
fit_zipoisson

fit_zipoisson <- glmmTMB(angui ~ 1 + NH4_av + T_av + (1|site),
                         data=data,
                         ziformula= ~ 1 + (1|site),
                         family=poisson)
summary(fit_zipoisson)
Anova(fit_zipoisson)



