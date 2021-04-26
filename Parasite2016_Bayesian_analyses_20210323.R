#ignore this first line


#The script runs on R v. 3.4.2

#This script codes the Bayesian Model Averging for the stickleback parasite study
#It also includes figures

library(ggplot2) # for making figures
library(ggeffects) # for making figures
library(vegan) # for function dispweight
library(MASS) # for step-wise model selection
library(car) # for VIF
library(BAS) # for Bayesian models
library(dplyr)


#######################
#### PRELIMINARIES ####
#######################

#### READ AND PREPARE DATA ####
#setwd('C:/Users/u0113095/Google Drive/PhD/2 Parasite/2016/Analysis_2020/Github/Data')
setwd('C:/Users/pascalh/Documents/GitHub/Stickleback-parasites-2016')
data <- read.csv("data_2016_2303.csv", sep=';')
data$site <- as.factor(data$site)
environment <- read.csv("Environment_update.csv", sep=';')
spavar <- read.csv("space2.csv", sep=';') #spatial variables: network centrality and upstream distance
plot(spavar$netcen); plot(density(spavar$netcen))
plot(spavar$updist); plot(density(spavar$updist))
field_data <- read.csv("field_data.csv", sep=',')
environment2 <- cbind(environment[,c(49,52:53,55,57,60,63)], field_data[-c(8,10,25,27,37),33:34], spavar[,c(2,3)])
environment2$pool_riffle <- as.factor(environment2$pool_riffle)
environment2$meander <- as.factor(environment2$meander)

#### Matrix for PIP (Posterior Inclusion Probability) ####
PIP <- matrix(nrow=12, ncol=14)
rownames(PIP) <- c("Condition", "Length", "Temperature", "Oxygen saturation", "Conductivity", "COD", "Ammonium", "Total nitrogen", "Pool riffle pattern", "Meander", "Network centrality", "Upstream distance")
colnames(PIP) <- c("Condition", "Length", "Gyrodactylus abundance", "Gyrodactylus prevalence", "Gyrodactylus infection intensity", "Trichodina abundance", "Trichodina prevalence", "Trichodina infection intensity", "Glugea", "Contracaecum", "Aguillicola",
                   "PI", "PI_ecto", "PI_endo")  

#### Condition ####

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


#### CALCULATE PARAMETERS ####
names(data)
#parasite data is overdispersed (mostly so for Trichodina), if using average abundance data, species matrix needs to be transformed
datao <- na.omit(data[,c(1,22:24,26:32)])
ddata <- dispweight(datao[,-1])
avab <- aggregate(ddata, by = list(datao[,1]), function(x){mean(x, na.rm =T)})
prev = aggregate(data[,c(22:24,26:32)], by = list(data[,1]), function(x){sum(x >0, na.rm = T)/length(x)})
medin = aggregate(data[,c(22:24,26:32)], by = list(data[,1]), function(x){median(x[x >0], na.rm = T)}) 
pa = aggregate(data[,c(22:24,26:32)], by = list(data[,1]), function(x){ifelse(mean(x, na.rm =T)>0,1,0)}) 

avab[is.na(avab)] <- 0
prev[is.na(prev)] <- 0
medin[is.na(medin)] <- 0

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
