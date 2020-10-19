library(devtools)
install_github("jslefche/piecewiseSEM@devel", build_vignette = TRUE)
library(piecewiseSEM)
library(nlme)
library(lme4)
library(nlme)
library(dplyr) # building data matrix
library(MASS) # for GLMMs
library(vegan)

setwd('C:/Users/pascalh/Documents/GitHub/Stickleback-parasites-2016')
#setwd('C:/Users/u0113095/Google Drive/PhD/2 Parasite/2016/Analysis_2020/data')

# Load data
data_2016 <- read.csv("data_2016.csv", sep=';') #field and parasite data
env <- read.csv("Environment_R.csv", sep=',') #all environmental data
env_av <- read.csv("Env_av.csv", sep=';') #environmental variables (average values)
env_max <- read.csv("env_max.csv", sep=';') #environmental variables (max. values)
spavar <- read.csv("space2.csv", sep=';') #spatial variables: network centrality and upstream distance
distance_matrix <- read.csv("distance_matrix.csv", sep=';') #spatial variables: distance matrix

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

data <- cbind(data_2016, env_av_exp, env_max_exp, spavar_exp)
data$site <- as.factor(data$site)
data$fish <- as.factor(data$fish)
data$length <- as.numeric(data$length)

data0 <- cbind(data_2016, env_av_exp[,-1], spavar_exp)
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
condition <- resid(lm(data$weight~data$length + data$Sex), na.action=na.exclude)
avcondition <- aggregate(condition, by = list(data[,1]), function(x){mean(x, na.rm =T)})[,2]

avin[is.na(avin)] <- 0
avab[is.na(avab)] <- 0
prev[is.na(prev)] <- 0
medin[is.na(medin)] <- 0

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


#add PI and condition index to data frame
data <- cbind(data, PI_ecto, PI_endo, PI)

#remove fish with NAs
data <- data[complete.cases(data[ , c("PI_ecto", "Sex", "length", "weight")]),]

#calculate condition index
confactor <- resid(lm(data$weight~data$length+data$Sex), na.action=na.exclude)
summary(confactor)

#add PI and condition index to data frame
data <- cbind(data, confactor)

avPI <- aggregate(data$PI, by = list(data[,1]), function(x){mean(x, na.rm =T)})[,2]



dataX <- as.data.frame(cbind(avab[,-c(1)], avcondition, avlength, avPI, env[,c(4,7,8,9,11,13,14,15)], spavar[,2:3]))
names(dataX)
plot(avPI)
model <- psem(
  lm(avPI ~ avcondition + avlength + T_av + con_av + O2_sat_av + COD_av + NH4_av + netcen + updist, data=dataX),
  lm(avlength ~ T_av + con_av + O2_sat_av + COD_av + NH4_av + netcen + updist, data=dataX), data=dataX)
summary(model)

model <- psem(
  lme(PI_ecto ~ Sex + confactor + sqrt(length) + sqrt(T_av) + sqrt(Temperature) + sqrt(O2_av) + sqrt(con_av) + sqrt(KjN_av) + sqrt(netcen) + sqrt(updist), random=~1|site, data=data),
  lme(confactor ~ Sex + sqrt(T_av) + sqrt(Temperature) + sqrt(O2_av) + sqrt(con_av) + sqrt(KjN_av) + sqrt(netcen) + sqrt(updist), random=~1|site, data=data),
  lme(sqrt(length) ~ Sex + sqrt(T_av) + sqrt(Temperature) + sqrt(O2_av) + sqrt(con_av) + sqrt(KjN_av) + sqrt(netcen) + sqrt(updist), random=~1|site, data=data)
)
summary(model)

model <- psem(
  lme(PI_endo ~ Sex + confactor + sqrt(length) + sqrt(T_av) + sqrt(Temperature) + sqrt(O2_av) + sqrt(con_av) + sqrt(KjN_av) + sqrt(netcen) + sqrt(updist), random=~1|site, data=data),
  lme(confactor ~ Sex + sqrt(T_av) + sqrt(Temperature) + sqrt(O2_av) + sqrt(con_av) + sqrt(KjN_av) + sqrt(netcen) + sqrt(updist), random=~1|site, data=data),
  lme(sqrt(length) ~ Sex + sqrt(T_av) + sqrt(Temperature) + sqrt(O2_av) + sqrt(con_av) + sqrt(KjN_av) + sqrt(netcen) + sqrt(updist), random=~1|site, data=data)
)
summary(model)

model <- psem(
  lme(PI ~ Sex + confactor + length + T_av + Temperature + O2_av + con_av + KjN_av + netcen + updist, random=~1|site, data=data),
  lme(confactor ~ Sex + T_av + Temperature + O2_av + con_av + KjN_av + netcen + updist, random=~1|site, data=data),
  lme(length ~ Sex + T_av + Temperature + O2_av + con_av + KjN_av + netcen + updist, random=~1|site, data=data)
)
summary(model)

data$length <- sqrt(data$length)
data$T_av <- sqrt(data$T_av)
data$Temperature <- sqrt(data$Temperature)
data$O2_av <- sqrt(data$O2_av)
data$con_av <- sqrt(data$con_av)
data$KjN_av <- sqrt(data$KjN_av)
data$netcen <- sqrt(data$netcen)
data$updist <- sqrt(data$updist)

model <- psem(
  lme(PI ~ Sex + confactor + length + T_av + Temperature + O2_av + con_av + KjN_av + netcen + updist, random=~1|site, data=data),
  lme(confactor ~ T_av + Temperature + O2_av + con_av + KjN_av + netcen + updist, random=~1|site, data=data)
)
summary(model)

model <- psem(
  lme(PI ~ confactor + T_av, random=~1|site, data=data),
  lme(confactor ~ T_av + Temperature + netcen, random=~1|site, data=data)
)
summary(model)

model <- psem(
  lme(PI ~ confactor + T_av, random=~1|site, data=data),
  lme(confactor ~ T_av + Temperature + netcen, random=~1|site, data=data)
)
summary(model)

model <- psem(
  lme(PI_ecto ~ confactor + con_av, random=~1|site, data=data),
  lme(confactor ~ T_av + Temperature + netcen, random=~1|site, data=data)
)
summary(model)

model <- psem(
  lme(PI_endo ~ confactor + Temperature + netcen, random=~1|site, data=data),
  lme(confactor ~ T_av + Temperature + netcen, random=~1|site, data=data)
)
summary(model)





model_PI <- lme(PI ~ Sex + confactor + length + T_av + Temperature + O2_av + con_av + KjN_av + netcen + updist, random=~1|site, data=data)
model_confactor <- lme(confactor ~ T_av + Temperature + O2_av + con_av + KjN_av + netcen + updist, random=~1|site, data=data)
summary(model_PI)
summary(model_confactor)

model <- psem(
  lme(PI_endo ~ Sex + confactor + length + T_av + Temperature + O2_av + con_av + KjN_av + netcen + updist, random=~1|site, data=data),
  lme(confactor ~ Sex + T_av + Temperature + O2_av + con_av + KjN_av + netcen + updist, random=~1|site, data=data),
  lme(length ~ Sex + T_av + Temperature + O2_av + con_av + KjN_av + netcen + updist, random=~1|site, data=data)
)
summary(model)


model <- psem(
  lm(prev$gyro ~ avcondition + avlength + T_av + con_av + COD_av + NO3_av, data=env),
  lm(avcondition ~ con_av, data=env),
)
summary(model)



?glmmPQL
model <- psem(
  lm(PI_ecto ~ Sex + confactor + sqrt(length) + sqrt(T_av) + sqrt(Temperature) + sqrt(O2_av) + sqrt(con_av) + sqrt(KjN_av) + sqrt(netcen) + sqrt(updist), data=data),
  lm(confactor ~ Sex + sqrt(T_av) + sqrt(Temperature) + sqrt(O2_av) + sqrt(con_av) + sqrt(KjN_av) + sqrt(netcen) + sqrt(updist), data=data)
  )
summary(model)

model <- psem(
  lme(PI_ecto ~ Sex + sqrt(updist), random=~1|site, data=data),
  lme(confactor ~ Sex + sqrt(updist), random=~1|site, data=data),
  lme(length ~ Sex + sqrt(updist), random=~1|site, data=data)
)
summary(model)

model <- psem(
  lme(PI_ecto ~ Sex + confactor + length + T_av + Temperature + O2_av + con_av + KjN_av + netcen + updist, random=~1|site, data=data, na.action=na.omit),
  lme(confactor ~ Sex + T_av + Temperature + O2_av + con_av + KjN_av + netcen + updist, random=~1|site, data=data, na.action=na.omit),
  lme(length ~ Sex + T_av + Temperature + O2_av + con_av + KjN_av + netcen + updist, random=~1|site, data=data, na.action=na.omit)
)
summary(model)



?psem
lme(PI_ecto ~  gyro ~ Sex + site, data=data, na.action=na.omit)
?lme
summary(PI_ecto)

data <- as.matrix(data)
summary(data)
?na.action
?psem






datat <- read.csv("data_2016_SEM.csv", sep=';')
#data$site <- as.factor(data$site)
#data$fish <- as.factor(data$fish)
#data$length <- as.numeric(data$length)
#data$updist2 <- as.numeric(data$updist2)

datat$PI_ecto2 = sign(datat$PI_ecto) * abs(datat$PI_ecto)^(1/3)
datat$PI_endo2 = sign(datat$PI_endo) * abs(datat$PI_endo)^(1/3)
datat$PI2 = sign(datat$PI) * abs(datat$PI)^(1/3)

summary(data$PI_ecto)
summary(datat$PI_ecto2)

model = psem(
  lme(PI_ecto2 ~  Nitrogen + pH + Oxygen + Conductivity + Temp + Phosphorus + length + Sex + Updist3 + nr_species, random=~1|site, data = datat),
  lme(confactor2 ~ Updist3 + nr_species + PI_ecto2 + Nitrogen + pH + Oxygen + Conductivity + Temp + Phosphorus + length + Sex, random=~1|site, 
      data = datat), 
  lme(length ~ pH + Nitrogen + Oxygen + Conductivity + Sex + Temp + Phosphorus + updist2 + nr_species, random=~1|site, data = datat))
summary(model)

model = psem(
  lm(PI_endo2 ~  Nitrogen + pH + Oxygen + Conductivity + Temp + Phosphorus + length + Sex, data = data),
  lm(confactor2 ~ PI_endo2 + length + Nitrogen + pH + Oxygen + Conductivity + Temp + Phosphorus + Sex, data = data), 
  lm(length ~ pH + Nitrogen + Oxygen + Conductivity + Sex + Temp + Phosphorus, data = data))
summary(model)$coefficients

model = psem(
  lme(PI_ecto2 ~  Nitrogen + pH + Oxygen + Conductivity + Temp + Phosphorus + length + Sex + Updist3 + nr_species, random=~1|site, data = data),
  lme(confactor2 ~ Updist3 + nr_species + PI_ecto2 + Nitrogen + pH + Oxygen + Conductivity + Temp + Phosphorus + length + Sex, random=~1|site, 
     data = data), 
  lme(length ~ pH + Nitrogen + Oxygen + Conductivity + Sex + Temp + Phosphorus + updist2 + nr_species, random=~1|site, data = data))
summary(model)

model = psem(
  lm(PI_endo2 ~  Nitrogen + pH + Oxygen + Conductivity + Temp + Phosphorus + length + Sex, data = data),
  lm(confactor2 ~ PI_endo2 + length + Nitrogen + pH + Oxygen + Conductivity + Temp + Phosphorus + Sex, data = data), 
  lm(length ~ pH + Nitrogen + Oxygen + Conductivity + Sex + Temp + Phosphorus, data = data))
summary(model)$coefficients


