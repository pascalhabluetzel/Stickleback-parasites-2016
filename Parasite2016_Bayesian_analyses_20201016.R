#The script runs on R v. 3.4.2

#This script codes the Bayesian Model Averging for the stickleback parasite study
#It also includes figures

library(ggplot2) # for making figures
library(ggeffects) # for making figures
#library(readr)
#library(gridExtra)
#library(plyr)
library(vegan) # for function dispweight
library(MASS) # for step-wise model selection
#library(lme4)
library(car) # for VIF
#library(BMS)
library(BAS) # for Bayesian models
library(dplyr)


#######################
#### PRELIMINARIES ####
#######################

#### READ AND PREPARE DATA ####
setwd('C:/Users/pascalh/Documents/GitHub/Stickleback-parasites-2016')
data <- read.csv("data_2016.csv", sep=';')
environment <- read.csv("Environment_update.csv", sep=';')

#### ENVIRONMENTAL DATA ####
# look at correlation among water parameters
colnames(environment)
plot(environment$T_av); plot(density(environment$T_av))
plot(environment$pH_av) # correlated with Con -> leave out
#plot(environment$O2_av) # highly correlated with O2 and dependent on T -> leave out
plot(environment$O2_sat_av); plot(density((environment$O2_sat_av)^2))
plot(environment$Con_av); plot(density((environment$Con_av)^2))
plot(environment$Cl._av) # correlated with Con -> leave out
plot(environment$COD_av); plot(density(environment$COD_av))
plot(environment$KjN_av) # highly correlated with NH4 -> leave out
plot(environment$NH4._av); plot(density(log(environment$NH4._av)))
plot(environment$NO3._av); # correlated with Nt -> leave out
plot(environment$NO2._av) # correlated with NO3 -> leave out
plot(environment$Nt_av); plot(density(log(environment$Nt_av))) # highly correlated with NO3 -> leave out
plot(environment$Pt_av) # highly correlated with oPO4 -> leave out
plot(environment$oPO4_av) # highly correlated with NH4 -> leave out
plot(environment$SM_av); plot(density(log(environment$SM_av)))

#### FIELD DATA ####
field_data <- read.csv("field_data.csv", sep=',')
environment2 <- cbind(environment[,c(49,52:53,55,57,60,63)], field_data[-c(8,10,25,27,37),33:34])
environment2$pool_riffle <- as.factor(environment2$pool_riffle)
environment2$meander <- as.factor(environment2$meander)

#### SPATIAL DATA ####
spavar <- read.csv("space2.csv", sep=';') #spatial variables: network centrality and upstream distance
plot(spavar$netcen); plot(density(spavar$netcen))
plot(spavar$updist); plot(density(spavar$updist))

#### Matrix for PIP (Posterior Inclusion Probability) ####
PIP <- matrix(nrow=12, ncol=11)
rownames(PIP) <- c("Condition", "Length", "Temperature", "Oxygen saturation", "Conductivity", "COD", "Ammonium", "Total nitrogen", "Pool riffle pattern", "Meander", "Network centrality", "Upstream distance")
colnames(PIP) <- c("Condition", "Length", "Gyrodactylus abundance", "Gyrodactylus prevalence", "Gyrodactylus infection intensity", "Trichodina abundance", "Trichodina prevalence", "Trichodina infection intensity", "Glugea", "Contracaecum", "Aguillicola")  

#### Condition ####
# Effect of environment (average) on host condition
condition <- resid(lm(data$weight~data$length + data$Sex), na.action=na.exclude)
datao <- na.omit(data[,c(1,18,19,21)])
avcondition <- aggregate(condition, by = list(datao[,1]), function(x){mean(x, na.rm =T)})[,2]
summary(avcondition); plot(density(avcondition))

model <- lm(avcondition ~ T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(Nt_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2)
summary(model)
#plot(model)
vif(model)

step.model <- stepAIC(lm(avcondition ~ T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(Nt_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2),
                      direction = "both", 
                      trace = FALSE)
step.model

summary(model <- lm(avcondition ~ pool_riffle + meander, data=environment2))
#plot(model)
vif(model)
res <- resid(lm(avcondition ~ pool_riffle, data=environment2))
boxplot(res ~ environment2$meander)
res <- resid(lm(avcondition ~ meander, data=environment2))
boxplot(res ~ environment2$meander)

#Bayesian approach
bas.model <- bas.lm(avcondition ~ T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(Nt_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2, prior="JZS")
summary(bas.model)
image(bas.model, rotate=F)
coef.model <- coef(bas.model)
abs(coef.model$postmean)-2*coef.model$postsd > 0
plot(confint(coef.model, parm = 2:11))

pip <- summary(bas.model)
PIP[c(3:12),1] <- pip[2:11,1]*sign(coef.model$postmean[2:11])

#### Length ####
# Effect of environment (average) on host length
avlength <- aggregate(data$length, by = list(data[,1]), function(x){mean(x, na.rm =T)})[,2]
summary(avlength); plot(density(avlength))

model <- lm(avlength ~ T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(Nt_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2)
summary(model)
#plot(model)
vif(model)

step.model <- stepAIC(lm(avlength ~ T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(Nt_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2),
                      direction = "both", 
                      trace = FALSE)
step.model

summary(model <- lm(avlength ~ Con_av + meander + spavar$netcen, data=environment2))
#plot(model)
vif(model)


# Have a look at a Bayesian approach
bas.model <- bas.lm(avlength ~ T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(Nt_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2, prior="JZS")
summary(bas.model)
image(bas.model, rotate=F)
coef.model <- coef(bas.model)
abs(coef.model$postmean)-2*coef.model$postsd > 0
plot(confint(coef.model, parm = 2:11))

pip <- summary(bas.model)
PIP[c(3:12),2] <- pip[2:11,1]*sign(coef.model$postmean[2:11])

#Figure
BPM <- predict(bas.model, estimator = "BPM", se.fit = TRUE)
BPM$best.vars
fit <- lm(avlength ~ O2_sat_av + Con_av^2 + spavar$netcen + spavar$updist, data=environment2)
d <- data.frame(cbind(spavar$netcen, fit$fitted.values))
p <- ggplot(d, aes(x = spavar$netcen, y = fit$fitted.values)) +
  theme_bw() +
  geom_smooth(method=lm, color="red", fill="#69b3a2", se=TRUE) +
  geom_smooth(method=loess, se=FALSE, linetype="dashed") +
  geom_point(data = environment2, aes(x = spavar$netcen, y = fit$fitted.values)) +
  labs(x=expression("Network peripherality [m]"), y=expression("Average host length [mm SL]"))
p

fit <- lm(avlength ~ O2_sat_av + Con_av^2 + spavar$netcen + spavar$updist, data=environment2)
d <- data.frame(cbind(environment2$Con_av, fit$fitted.values))
p <- ggplot(d, aes(x = environment2$Con_av, y = fit$fitted.values)) +
  theme_bw() +
  geom_smooth(method=lm, color="red", fill="#69b3a2", se=TRUE) +
  geom_smooth(method=loess, se=FALSE, linetype="dashed") +
  geom_point(data = environment2, aes(x = environment2$Con_av, y = fit$fitted.values)) +
  labs(x=expression("Conductivity ["*mu*"S/cm]"), y=expression("Average host length [mm SL]"))
p



#### CALCULATE PARAMETERS ####
names(data)
#parasite data is overdispersed (mostly so for Trichodina), if using average abundance data, species matrix needs to be transformed
datao <- na.omit(data[,c(1,23:25,27:33)])
ddata <- dispweight(datao[,-1])
avab <- aggregate(ddata, by = list(datao[,1]), function(x){mean(x, na.rm =T)})
prev = aggregate(data[,c(23:25,27:33)], by = list(data[,1]), function(x){sum(x >0, na.rm = T)/length(x)})
medin = aggregate(data[,c(23:25,27:33)], by = list(data[,1]), function(x){median(x[x >0], na.rm = T)}) 
pa = aggregate(data[,c(23:25,27:33)], by = list(data[,1]), function(x){ifelse(mean(x, na.rm =T)>0,1,0)}) 

avab[is.na(avab)] <- 0
prev[is.na(prev)] <- 0
medin[is.na(medin)] <- 0

#### Gyrodactylus ####
#### Average abundance ####
plot(avab$gyro); plot(density(avab$gyro))
plot((avab$gyro)^(1/3)); plot(density((avab$gyro)^(1/3)))

model <- lm((avab$gyro)^(1/3) ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(Nt_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2)
summary(model)
#plot(model)                      

step.model <- stepAIC(lm((avab$gyro)^(1/3) ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(NO3._av) + log(SM_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2),
                      direction = "both", 
                      trace = FALSE)
step.model

model <- lm((avab$gyro)^(1/3) ~ avlength + T_av + COD_av + log(Nt_av) + pool_riffle, data=environment2)
summary(model)
#plot(model)                      
res <- resid(lm((avab$gyro)^(1/3) ~ T_av + COD_av + log(Nt_av) + pool_riffle, data=environment2))
plot(res ~ avlength)
lines(lowess(res ~ avlength), col=3)
cor.test(res, avlength, method="spearman")
res <- resid(lm((avab$gyro)^(1/3) ~ avlength + T_av + log(Nt_av) + pool_riffle, data=environment2))
plot(res ~ environment2$COD_av)
lines(lowess(res ~ environment2$COD_av), col=3)
cor.test(res, environment2$COD_av, method="spearman")


# Have a look at a Bayesian approach
bas.model <- bas.lm(avab$gyro ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(Nt_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2, prior="JZS")
summary(bas.model)
image(bas.model, rotate=F)
coef.model <- coef(bas.model)
abs(coef.model$postmean)-2*coef.model$postsd > 0
confint(coef.model)
plot(confint(coef.model, parm = 2:12))

pip <- summary(bas.model)
PIP[c(1:12),3] <- pip[2:13,1]*sign(coef.model$postmean[2:13])


#### Prevalence ####
plot(prev$gyro); plot(density(prev$gyro))

model <- lm(prev$gyro ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(Nt_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2)
summary(model)
#plot(model)                      

step.model <- stepAIC(lm(prev$gyro ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(NO3._av) + log(SM_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2),
                      direction = "both", 
                      trace = FALSE)
step.model

model <- lm(prev$gyro ~ Con_av + COD_av + pool_riffle + meander, data=environment2)
summary(model)
#plot(model)                      
res <- resid(lm(prev$gyro ~ COD_av + pool_riffle + meander, data=environment2))
plot(res ~ environment$Con_av)
lines(lowess(res ~ environment$Con_av), col=3)
cor.test(res, environment2$Con_av, method="spearman")
res <- resid(lm(prev$gyro ~ Con_av + pool_riffle + meander, data=environment2))
plot(res ~ environment$COD_av)
lines(lowess(res ~ environment$COD_av), col=3)
cor.test(res, environment2$COD_av, method="spearman")

# Have a look at a Bayesian approach
bas.model <- bas.lm(prev$gyro ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(Nt_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2, prior="JZS")
summary(bas.model)
image(bas.model, rotate=F)
coef.model <- coef(bas.model)
abs(coef.model$postmean)-2*coef.model$postsd > 0
confint(coef.model)
plot(confint(coef.model, parm = 2:12))

pip <- summary(bas.model)
PIP[c(1:12),4] <- pip[2:13,1]*sign(coef.model$postmean[2:13])

#### Median abundance ####
plot(medin$gyro); plot(density(medin$gyro))
plot((medin$gyro)^(1/3)); plot(density((medin$gyro)^(1/3)))

model <- lm((medin$gyro)^(1/3) ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(Nt_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2)
summary(model)
#plot(model)                      

step.model <- stepAIC(lm((medin$gyro)^(1/3) ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(NO3._av) + log(SM_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2),
                      direction = "both", 
                      trace = FALSE)
step.model

model <- lm((medin$gyro)^(1/3) ~ T_av + COD_av + spavar$netcen, data=environment2)
summary(model)
#plot(model)                      


# Have a look at a Bayesian approach
bas.model <- bas.lm((medin$gyro)^(1/3) ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(Nt_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2, prior="JZS")
summary(bas.model)
image(bas.model, rotate=F)
coef.model <- coef(bas.model)
abs(coef.model$postmean)-2*coef.model$postsd > 0
confint(coef.model)
plot(confint(coef.model, parm = 2:12))

pip <- summary(bas.model)
PIP[c(1:12),5] <- pip[2:13,1]*sign(coef.model$postmean[2:13])

#### Trichodina ####
#### Average abundance ####
plot(avab$tricho); plot(density(avab$tricho))
plot((avab$tricho)^(1/3)); plot(density((avab$tricho)^(1/3)))

model <- lm((avab$tricho)^(1/3) ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(NO3._av) + log(SM_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2)
summary(model)
#plot(model)                      

step.model <- stepAIC(lm((avab$tricho)^(1/3) ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(NO3._av) + log(SM_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2),
                      direction = "both", 
                      trace = FALSE)
step.model

model <- lm((avab$tricho)^(1/3) ~ Con_av + COD_av + pool_riffle + meander + spavar$netcen, data=environment2)
summary(model)
#plot(model)                      
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

# Bayesian
bas.model <- bas.lm((avab$tricho)^(1/3) ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(Nt_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2, prior="JZS")
summary(bas.model)
image(bas.model, rotate=F)
coef.model <- coef(bas.model)
abs(coef.model$postmean)-2*coef.model$postsd > 0
confint(coef.model)
plot(confint(coef.model, parm = 2:12))

pip <- summary(bas.model)
PIP[c(1:12),6] <- pip[2:13,1]*sign(coef.model$postmean[2:13])

#### Prevalence ####
plot(prev$tricho); plot(density(prev$tricho))

model <- lm(prev$tricho ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(NO3._av) + log(SM_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2)
summary(model)
#plot(model)                      

step.model <- stepAIC(lm(prev$tricho ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(NO3._av) + log(SM_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2),
                      direction = "both", 
                      trace = FALSE)
step.model

model <- lm(prev$tricho ~ Con_av + log(NH4._av), data=environment2)
summary(model)
#plot(model)                      
res <- resid(lm(prev$tricho ~ log(NH4._av), data=environment2))
plot(res ~ environment2$Con_av)
lines(lowess(res ~ environment2$Con_av), col=3)
cor.test(res, environment2$Con_av, method="spearman")
res <- resid(lm(prev$tricho ~ Con_av, data=environment2))
plot(res ~ log(environment2$NH4._av))
lines(lowess(res ~ log(environment2$NH4._av)), col=3)
cor.test(res, environment2$NH4._av, method="spearman")

# Bayesian
bas.model <- bas.lm(prev$tricho ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(Nt_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2, prior="JZS")
summary(bas.model)
image(bas.model, rotate=F)
coef.model <- coef(bas.model)
abs(coef.model$postmean)-2*coef.model$postsd > 0
confint(coef.model)
plot(confint(coef.model, parm = 2:12))

pip <- summary(bas.model)
PIP[c(1:12),7] <- pip[2:13,1]*sign(coef.model$postmean[2:13])


#### Median abundance ####
plot(medin$tricho); plot(density(medin$tricho))
plot((medin$tricho)^(1/3)); plot(density((medin$tricho)^(1/3)))

model <- lm((medin$tricho)^(1/3) ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(Nt_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2)
summary(model)
#plot(model)                      

step.model <- stepAIC(lm((medin$tricho)^(1/3) ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(Nt_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2),
                      direction = "both", 
                      trace = FALSE)
step.model

model <- lm((medin$tricho)^(1/3) ~ Con_av^2 + COD_av + log(NH4._av) + pool_riffle + meander + spavar$netcen, data=environment2)
summary(model)
#plot(model)                      
res <- resid(lm(prev$tricho ~ COD_av + log(NH4._av) + pool_riffle + meander + spavar$netcen, data=environment2))
plot(res ~ environment2$Con_av)
lines(lowess(res ~ environment2$Con_av), col=3)
cor.test(res, environment2$Con_av^2, method="spearman")

# Bayesian
bas.model <- bas.lm((medin$tricho)^(1/3) ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(Nt_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2, prior="JZS")
summary(bas.model)
image(bas.model, rotate=F)
coef.model <- coef(bas.model)
abs(coef.model$postmean)-2*coef.model$postsd > 0
confint(coef.model)
plot(confint(coef.model, parm = 2:12))

pip <- summary(bas.model)
PIP[c(1:12),8] <- pip[2:13,1]*sign(coef.model$postmean[2:13])

BPM <- predict(bas.model, estimator = "BPM")
variable.names(BPM)


BPM <- predict(bas.model, estimator = "BPM", se.fit = TRUE)
conf.fit <- confint(BPM, parm = "mean")
conf.pred <- confint(BPM, parm = "pred")
plot(conf.pred[,3]^3 ~ environment2$Con_av, xlab=expression(paste("Conductivity [", mu, "S/cm]")), ylab=expression(paste(italic("Trichodina"), " sp. [median infection intensity]")))
lines(lowess(conf.pred[,3]^3 ~ environment2$Con_av), lt=2)
abline(lm(conf.pred[,3]^3 ~ poly(environment2$Con_av, degree=-3)), lt=2)
abline(glm(conf.pred[,1]^3 ~ environment2$Con_av), lt=2)
abline(glm(conf.pred[,2]^3 ~ environment2$Con_av), lt=2)
# for comparison, plot the raw data
plot(medin$tricho ~ environment2$Con_av, xlab=expression(paste("Conductivity [", mu, "S/cm]")), ylab=expression(paste(italic("Trichodina"), " sp. [median infection intensity]")))
lines(lowess(medin$tricho ~ environment2$Con_av), lt=2)

fit <- lm(medin$tricho ~ T_av + Con_av + COD_av + NH4._av + spavar$updist, data=environment2)
d <- data.frame(cbind(environment2$Con_av, fit$fitted.values))
p <- ggplot(d, aes(x = environment2$Con_av, y = fit$fitted.values)) +
  theme_bw() +
  geom_smooth(method=lm, color="red", fill="#69b3a2", se=TRUE) +
  geom_point(data = environment2, aes(x = environment2$Con_av, y = fit$fitted.values))
p

library(ggplot2)
fit <- lm((medin$tricho)^(1/3) ~ T_av + Con_av^2 + COD_av + log(NH4._av) + spavar$updist, data=environment2)
fit <- lm((medin$tricho)^(1/3) ~ T_av + COD_av + log(NH4._av) + spavar$updist, data=environment2)
newdata <- as.data.frame(cbind(rep(mean(environment2$T_av),37), seq(from=range(environment2$Con_av)[1], to=range(environment2$Con_av)[2], length.out=37), rep(mean(environment2$COD_av),37), rep(mean(environment2$NH4._av),37), rep(mean(spavar$updist),37)))
colnames(newdata) <- c("T_av", "Con_av", "COD_av", "NH4._av", "spavar$updist")
summary(fit)
pred <- predict(fit)
pred <- predict(fit, newdata)
d <- data.frame(cbind(environment2$Con_av, (pred)^3))
p <- ggplot(d, aes(x = newdata$Con_av, y = (pred)^3)) +
  theme_bw() +
  geom_smooth(method=lm, color="red", fill="#69b3a2", se=TRUE) +
  geom_smooth(method=loess, se=FALSE, linetype="dashed") +
  geom_point(data = environment2, aes(x = environment2$Con_av, y = (pred)^3)) +
  labs(x=expression("Conductivity ["*mu*"S/cm]"), y=expression(italic(Trichodina)~"sp. [median infection intensity]"))
p

# This is the right graphic #

#install.packages("ggeffects")
#library(ggeffects)
updist <- spavar$updist
fit <- lm((medin$tricho)^(1/3) ~ T_av + Con_av^2 + COD_av + log(NH4._av) + updist, data=environment2)
predict <- ggpredict(fit, terms = "Con_av")
ggplot(predict, aes(x, predicted^3)) +
  theme_bw() +
  geom_line(color="red") +
  geom_ribbon(aes(ymin = conf.low^3, ymax = conf.high^3), alpha = .1) +
  geom_point(data = environment2, aes(x=environment2$Con_av, y=medin$tricho)) +
  labs(x=expression("Conductivity ["*mu*"S/cm]"), y=expression(italic(Trichodina)~"sp. [median infection intensity]"))

updist <- spavar$updist
fit <- lm((medin$tricho)^(1/3) ~ T_av + Con_av^2 + COD_av + log1p(NH4._av) + updist, data=environment2) # add a small number to the ammonium value to avoid infinite values in subsequent calculations
predict <- ggpredict(fit, terms = "NH4._av")
?log1p
#predicted[1] <- 0; predict$conf.high[1] <- 0; predict$conf.low[1] <- 0
ggplot(predict, aes(x, predicted^3)) +
  theme_bw() +
  geom_line(color="red") +
  geom_ribbon(aes(ymin = conf.low^3, ymax = conf.high^3), alpha = .1) +
  geom_point(data = environment2, aes(x=environment2$NH4._av, y=medin$tricho)) +
  labs(x=expression("NH4"), y=expression(italic(Trichodina)~"sp. [median infection intensity]"))

# This is the end of the right graphic #


fit1 = lm( log(mpg) ~ disp + hp + drat + wt, data = mtcars)
summary(fit1)

( mod_vars = all.vars( formula(fit1) )[-1] )

preddat_fun = function(data, allvars, var) {
  sums = summarise_at(data, 
                      vars( one_of(allvars), -one_of(var) ), 
                      median) 
  cbind( select_at(data, var), sums)
}

head( preddat_fun(mtcars, mod_vars, "disp") )

pred_dats = mod_vars %>%
  set_names() %>%
  map( ~preddat_fun(mtcars, mod_vars, .x) )
str(pred_dats)

preds = pred_dats %>%
  map(~augment(fit1, newdata = .x) ) %>%
  map(~mutate(.x, 
              lower = exp(.fitted - 2*.se.fit),
              upper = exp(.fitted + 2*.se.fit),
              pred = exp(.fitted) ) )

str(preds$disp)
ggplot(data = preds$disp, aes(x = disp, y = pred) ) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
  geom_rug(sides = "b") +
  theme_bw(base_size = 14) +
  labs(x = "Displacement (cu.in.)",
       y = "Miles/(US) gallon") +
  ylim(10, 32)

xlabs = c("Displacement (cu.in.)", 
          "Gross horsepower",
          "Rear axle ratio", 
          "Weight (1000 lbs)")

pred_plot = function(data, variable, xlab) {
  ggplot(data, aes(x = .data[[variable]], y = pred) ) +
    geom_line(size = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .25) +
    geom_rug(sides = "b") +
    theme_bw(base_size = 14) +
    labs(x = xlab,
         y = "Miles/(US) gallon") +
    ylim(10, 32)
}

pred_plot(preds[[1]], mod_vars[1], xlabs[1])

all_plots = pmap( list(preds, mod_vars, xlabs), pred_plot)
all_plots

cowplot::plot_grid(plotlist = all_plots,
                   labels = "AUTO",
                   align = "hv")



newdata <- as.data.frame(cbind(rep(mean(environment2$T_av),37), seq(from=range(environment2$Con_av)[1], to=range(environment2$Con_av)[2], length.out=37), rep(mean(environment2$COD_av),37), rep(mean(environment2$NH4._av),37), rep(mean(spavar$updist),37)))
colnames(newdata) <- c("T_av", "Con_av", "COD_av", "NH4._av", "spavar$updist")
summary(fit)
pred <- predict(fit)
pred <- predict(fit, newdata)
d <- data.frame(cbind(environment2$Con_av, (pred)^3))
p <- ggplot(d, aes(x = newdata$Con_av, y = (pred)^3)) +
  theme_bw() +
  geom_smooth(method=lm, color="red", fill="#69b3a2", se=TRUE) +
  geom_smooth(method=loess, se=FALSE, linetype="dashed") +
  geom_point(data = environment2, aes(x = environment2$Con_av, y = (pred)^3)) +
  labs(x=expression("Conductivity ["*mu*"S/cm]"), y=expression(italic(Trichodina)~"sp. [median infection intensity]"))
p


data1 <- as.data.frame(cbind((medin$tricho)^(1/3), environment2$T_av, environment2$Con_av^2, environment2$COD_av, log(environment2$NH4._av), spavar$updist))
colnames(data1) <- c("tricho", "T_av", "Con_av", "COD_av", "NH4._av", "updist")
fit <- lm(tricho ~ T_av + Con_av + COD_av + NH4._av + updist, data=data1)
effect_plot(model = fit, pred = Con_av, interval = TRUE, partial.residuals = TRUE)
??effect_plot
install.packages("effects")
library(effects)

e1.lm1 <- predictorEffect("Con_av", fit)
plot(e1.lm1)
install.packages("jtools")
library(jtools)

?effect_plot

BPM <- predict(bas.model, estimator = "BPM")
variable.names(BPM)
BPM <- predict(bas.model, estimator = "BPM", se.fit = TRUE)
conf.fit <- confint(BPM, parm = "mean")
conf.pred <- confint(BPM, parm = "pred")
plot(conf.pred[,3]^3 ~ environment2$NH4._av, xlab=expression(paste("NH4 [mg/l]")), ylab=expression(paste(italic("Trichodina"), " sp. [median infection intensity]")))
lines(lowess(conf.pred[,3]^3 ~ environment2$NH4._av), lt=2)
# for comparison, plot the raw data
plot(medin$tricho ~ environment2$NH4._av, xlab=expression(paste("NH4 [mg/l]")), ylab=expression(paste(italic("Trichodina"), " sp. [median infection intensity]")))
lines(lowess(medin$tricho ~ environment2$NH4._av), lt=2)

#### Glugea ####
model <- glm(pa$glugea ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(NO3._av) + log(SM_av) + meander, data=environment2, family=binomial(link="logit"))
summary(model)
#plot(model)   

model <- glm(pa$glugea ~ pool_riffle + spavar$netcen + spavar$updist, data=environment2, family=binomial(link="logit"))
summary(model)
#plot(model)   


step.model <- stepAIC(glm(pa$glugea ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(NO3._av) + log(SM_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2, family=binomial(link="logit")),
                      direction = "both", 
                      trace = FALSE)
step.model

model <- glm(pa$glugea ~ avlength + T_av + Con_av^2 + log(NH4._av) + pool_riffle + meander, data=environment2, family=binomial(link="logit"))
summary(model)
#plot(model) 

plot(as.factor(pa$glugea), avlength)
plot(as.factor(pa$glugea), environment2$Con_av)
table(as.factor(pa$glugea), environment2$meander)

bas.model <- bas.lm(pa$glugea ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(Nt_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2, prior="JZS")
summary(bas.model)
image(bas.model, rotate=F)
coef.model <- coef(bas.model)
abs(coef.model$postmean)-2*coef.model$postsd > 0
confint(coef.model)
plot(confint(coef.model, parm = 2:12))

pip <- summary(bas.model)
PIP[c(1:12),9] <- pip[2:13,1]*sign(coef.model$postmean[2:13])


#### Contracaecum ####
model <- glm(pa$contra ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(NO3._av) + log(SM_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2, family=binomial(link="logit"))
summary(model)
#plot(model)   

step.model <- stepAIC(glm(pa$contra ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(NO3._av) + log(SM_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2, family=binomial(link="logit")),
                      direction = "both", 
                      trace = FALSE)
step.model

model <- glm(pa$contra ~ avlength + O2_sat_av + Con_av^2 + log(NO3._av) + spavar$netcen, data=environment2, family=binomial(link="logit"))
summary(model)
#plot(model) 

plot(as.factor(pa$contra), avlength)
plot(as.factor(pa$contra), environment2$O2_sat_av)

bas.model <- bas.lm(pa$contra ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(Nt_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2, prior="JZS")
summary(bas.model)
image(bas.model, rotate=F)
coef.model <- coef(bas.model)
abs(coef.model$postmean)-2*coef.model$postsd > 0
confint(coef.model)
plot(confint(coef.model, parm = 2:12))

pip <- summary(bas.model)
PIP[c(1:12),10] <- pip[2:13,1]*sign(coef.model$postmean[2:13])

#### Anguillicola ####
model <- glm(pa$angui ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(NO3._av) + log(SM_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2, family=binomial(link="logit"))
summary(model)
#plot(model)   

# Model does not converge. Remove one predictor predictors.
model <- glm(pa$angui ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(NO3._av) + log(SM_av) + pool_riffle + meander, data=environment2, family=binomial(link="logit"))
summary(model)
#plot(model) 

step.model <- stepAIC(glm(pa$angui ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(NO3._av) + log(SM_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2, family=binomial(link="logit")),
                      direction = "both", 
                      trace = FALSE)
step.model

model <- glm(pa$angui ~ avlength + avcondition + T_av + Con_av^2 + COD_av + log(NH4._av) + log(SM_av) + pool_riffle + meander, data=environment2, family=binomial(link="logit"))
summary(model)
#plot(model) 

plot(as.factor(pa$angui), avlength)
plot(as.factor(pa$angui), environment2$T_av)
table(as.factor(pa$angui), environment2$meander)

bas.model <- bas.glm(pa$angui ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(Nt_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2, betaprior=g.prior(100), family=binomial)
summary(bas.model)
image(bas.model, rotate=F)
coef.model <- coef(bas.model)
abs(coef.model$postmean)-2*coef.model$postsd > 0
confint(coef.model)
plot(confint(coef.model, parm = 2:12))

pip <- summary(bas.model)
PIP[c(1:12),11] <- pip[2:13,1]*sign(coef.model$postmean[2:13])



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


avPI <- aggregate(PI, by = list(data[,1]), function(x){mean(x, na.rm =T)})[,2]
avPI_ecto <- aggregate(PI_ecto, by = list(data[,1]), function(x){mean(x, na.rm =T)})[,2]
avPI_endo <- aggregate(PI_endo, by = list(data[,1]), function(x){mean(x, na.rm =T)})[,2]


# Bayesian
bas.model <- bas.lm(avPI ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(Nt_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2, prior="JZS")
summary(bas.model)
image(bas.model, rotate=F)
coef.model <- coef(bas.model)
abs(coef.model$postmean)-2*coef.model$postsd > 0
confint(coef.model)
plot(confint(coef.model, parm = 2:12))

pip <- summary(bas.model)
PIP <- cbind(PIP,pip[2:13,1]*sign(coef.model$postmean[2:13]))


# Bayesian
bas.model <- bas.lm(avPI_ecto ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(Nt_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2, prior="JZS")
summary(bas.model)
image(bas.model, rotate=F)
coef.model <- coef(bas.model)
abs(coef.model$postmean)-2*coef.model$postsd > 0
confint(coef.model)
plot(confint(coef.model, parm = 2:12))

pip <- summary(bas.model)
PIP <- cbind(PIP,pip[2:13,1]*sign(coef.model$postmean[2:13]))


# Bayesian
bas.model <- bas.lm(avPI_endo ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(Nt_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2, prior="JZS")
summary(bas.model)
image(bas.model, rotate=F)
coef.model <- coef(bas.model)
abs(coef.model$postmean)-2*coef.model$postsd > 0
confint(coef.model)
plot(confint(coef.model, parm = 2:12))

pip <- summary(bas.model)
PIP <- cbind(PIP,pip[2:13,1]*sign(coef.model$postmean[2:13]))

write.matrix(PIP, file="PIP.txt", sep="\t")






# Bayesian Belief network
install.packages("bnlearn")
library(bnlearn)
bbn.data <- as.data.frame(cbind((medin$tricho)^1/3,avcondition,avlength,environment2,spavar$netcen,spavar$updist))
res = hc(bbn.data[,c("(medin$tricho)^1/3","avcondition","Con_av","COD_av")])
pvalues = arc.strength(res, data = bbn.data[,c("(medin$tricho)^1/3","avcondition","Con_av","COD_av")])
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rgraphviz")
BiocManager::install("RBGL")
library(gRain)
library(Rgraphviz)
strength.plot(res, strength = pvalues)

?hc


step.model <- stepAIC(glm(pa$angui ~ avlength + avcondition + T_av + O2_sat_av + Con_av^2 + COD_av + log(NH4._av) + log(NO3._av) + log(SM_av) + pool_riffle + meander + spavar$netcen + spavar$updist, data=environment2, family=binomial(link="logit")),
                      direction = "both", 
                      trace = FALSE)
step.model

model <- glm(pa$angui ~ avlength + O2_sat_av + Con_av^2 + log(NO3._av) + spavar$netcen, data=environment2, family=binomial(link="logit"))
summary(model)
#plot(model) 

plot(as.factor(pa$contra), avlength)
plot(as.factor(pa$contra), environment2$O2_sat_av)





res <- resid(glm(pa$contra ~ avlength + O2_sat_av + Con_av^2 + log(NO3._av) + spavar$netcen, data=environment2, family=binomial(link="logit")))
plot(res ~ environment$Con_av)
lines(lowess(res ~ environment$Con_av), col=3)





summary(model <- lm(avcondition[-2] ~ con_av + NO2_av + spavar$netcen[-2], data=env[-2,]))
#plot(model)
#Removing site 2 changes outcome of model -> rerun without location 2

step.model <- stepAIC(lm(avcondition[-2] ~ T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + spavar$netcen[-2] + spavar$updist[-2], data=env[-2,]),
                      direction = "both", 
                      trace = FALSE)
step.model

summary(model <- lm(avcondition[-2]^3 ~ O2_sat_av + NH4_av + NO2_av + spavar$updist[-2], data=env[-2,]))
#plot(model)

#model improved, reiterate selection with new approach

step.model <- stepAIC(lm(avcondition[-2]^3 ~ T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + spavar$netcen[-2] + spavar$updist[-2], data=env[-2,]),
                      direction = "both", 
                      trace = FALSE)
step.model

summary(model <- lm(avcondition[-c(2,26,34)]^3 ~ O2_sat_av + NH4_av + NO2_av + spavar$netcen[-c(2,26,34)] + spavar$updist[-c(2,26,34)], data=env[-c(2,26,34),]))
plot(model)
#Removing site 2 changes outcome of model -> rerun without location 2


res_O2_sat_av <- resid(lm(avcondition[-2] ~ NH4_av + NO2_av + spavar$netcen[-2] + spavar$updist[-2], data=env[-2,]))
plot(res_O2_sat_av ~ O2_sat_av, data=env[-2,])
abline(lm(res_O2_sat_av ~ O2_sat_av, data=env[-2,]))

res_NH4_av <- resid(lm(avcondition[-2] ~ O2_sat_av + NO2_av + spavar$netcen[-2] + spavar$updist[-2], data=env[-2,]))
plot(res_NH4_av ~ NH4_av, data=env[-2,])
abline(lm(res_NH4_av ~ NH4_av, data=env[-2,]))

res_NO2_av <- resid(lm(avcondition[-2] ~ O2_sat_av + NH4_av + spavar$netcen[-2] + spavar$updist[-2], data=env[-2,]))
plot(res_NO2_av ~ NO2_av, data=env[-2,])
abline(lm(res_NO2_av ~ NO2_av, data=env[-2,]))

res_updist <- resid(lm(avcondition[-2] ~ O2_sat_av + + NO2_av + NH4_av + spavar$netcen[-2], data=env[-2,]))
plot(res_updist ~  spavar$updist[-2])
abline(lm(res_updist ~  spavar$updist[-2]))


avlength <- aggregate(data$length, by = list(data[,1]), function(x){mean(x, na.rm =T)})[,2]



select <- buildglmmTMB(condition ~ T_av + con_av + O2_sat_av + Cl_av + COD_av + NH4_av + NO3_av + NO2_av + netcen + updist,
                       data=data,
                       include = ~ (1|site),
                       direction = c('order', 'forward'),
                       crit="LRT")
select

model <- lme(condition ~ con_av, random=~1|site, data=data, na.action=na.omit)
summary(model)



scatterplotMatrix(environment[,c(49,52:53,55,57:58,63)])
cor(environment[,c(49,52:53,55,57:58,63)])

data$site <- as.factor(data$site)
data$fish <- as.factor(data$fish)
data$length <- as.numeric(data$length)



#### CALCULATE PARAMETERS ####
avin = aggregate(data[,c(11:13,15:20)], by = list(data[,2]), function(x){mean(x[x >0], na.rm = T)}) 
avab = aggregate(data[,c(11:13,15:20)], by = list(data[,2]), function(x){mean(x, na.rm =T)})
prev = aggregate(data[,c(11:13,15:20)], by = list(data[,2]), function(x){sum(x >0, na.rm = T)/length(x)})
medin = aggregate(data[,c(11:13,15:20)], by = list(data[,2]), function(x){median(x[x >0], na.rm = T)}) 

avin
