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


#######################
#### PRELIMINARIES ####
#######################

#### READ AND PREPARE DATA ####
setwd('C:/Users/u0113095/Google Drive/PhD/2 Parasite/2016/Parasite2016_analysis')
data <- read.csv("data_2016.csv", sep=';')
data$site <- as.factor(data$site)
data$fish <- as.factor(data$fish)
data$length <- as.numeric(data$length)

#### CALCULATE PARAMETERS ####
avin = aggregate(data[,c(11:13,15:20)], by = list(data[,2]), function(x){mean(x[x >0], na.rm = T)}) 
avab = aggregate(data[,c(11:13,15:20)], by = list(data[,2]), function(x){mean(x, na.rm =T)})
prev = aggregate(data[,c(11:13,15:20)], by = list(data[,2]), function(x){sum(x >0, na.rm = T)/length(x)})
medin = aggregate(data[,c(11:13,15:20)], by = list(data[,2]), function(x){median(x[x >0], na.rm = T)}) 

#write.table(avin, "C:/Users/u0113095/Google Drive/PhD/2 Parasite/2016/_Analysis_2020/avin.txt", sep ='\t')
#write.table(avab, "C:/Users/u0113095/Google Drive/PhD/2 Parasite/2016/_Analysis_2020/avab.txt", sep ='\t')
#write.table(prev, "C:/Users/u0113095/Google Drive/PhD/2 Parasite/2016/_Analysis_2020/prev.txt", sep ='\t')
#write.table(medin, "C:/Users/u0113095/Google Drive/PhD/2 Parasite/2016/_Analysis_2020/medin.txt", sep ='\t')


##################################
#### 1. INDIVIDUAL BASED MODELS #
#################################

#### 1.1 INFECTION PRESENCE ####
data$tricho_pre=1 *(data$tricho >0) 
data$gyro_pre=1*(data$gyro >0)
data$glu_pre=1*(data$glugea>0)
data$con_pre=1*(data$contra>0)
data$angui_pre=1*(data$angui>0)

# glmer
library('lme4')

glmer_tri <- glmer(tricho_pre ~ site + Sex + length +(1|ecto_screener),family=binomial(link=logit),data=data, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)) )
summary(glmer_tri)
Anova(glmer_tri)

glmer_gy <- glmer(gyro_pre~ site + Sex + length +  (1|ecto_screener),family=binomial(link=logit),data=data)
summary(glmer_gy)
Anova(glmer_gy, type="III")

glmer_glu <- glmer(glu_pre~ site + Sex + length +  (1|ecto_screener),family=binomial(link=logit),data=data, control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
summary(glmer_glu)
Anova(glmer_glu, type="III")

glmer_con <- glmer(contra_pre~ site + Sex + length +  (1|endo_screener),family=binomial(link=cloglog),data=data)
summary(glmer_con)

Anova(glmer_con, type="III")
glmer_an <- glmer(angui_pre~ site + Sex + length+ (1|endo_screener),family=binomial(link=cloglog),data=data)
summary(glmer_an)
Anova(glmer_an, type="III")

# glmmPQL
library(MASS)

tr <- glmmPQL(tricho_pre ~ site + Sex + length + confactor,  ~1|ecto_screener, data=data, family = binomial)
summary(tr)
Anova(tr)

gy  <- glmmPQL(gyro_pre ~ site + Sex + length + confactor,  ~1|ecto_screener, data=data, family = binomial)
Anova(gy)

glu  <- glmmPQL(glu_pre ~ site + Sex + length + confactor,  ~1|ecto_screener, data=data, family = binomial)
Anova(glu)

con  <- glmmPQL(contra_pre ~ site + Sex + length + confactor,  ~1|endo_screener, data=data, family = binomial)
Anova(con)

Angui  <- glmmPQL(angui_pre ~ site + Sex + length + confactor,  ~1|endo_screener, data=data, family = binomial)
Anova(Angui)


#### 1.2 ABUNDANCE ####

# ZIGLMM (glmmTMB package)
library(glmmTMB)

fit_zipoisson <- glmmTMB(gyro ~ site + Sex + length + confactor + (1|ecto_screener),
                         data=data,
                         ziformula=~1,
                         family=poisson)
summary(fit_zipoisson)
Anova(fit_zipoisson)

fit_zinbinom <- glmmTMB(gyro ~ site + Sex + length + confactor + (1|ecto_screener),
                         data=data,
                         ziformula=~1,
                         family=nbinom1("log"))
summary(fit_zinbinom)
Anova(fit_zinbinom)

glu_zinbinom <- glmmTMB(glugea ~ site + Sex + length + confactor + (1|ecto_screener),
                        data=data,
                        ziformula=~1,
                        family=nbinom1("log"))
glu_zipoisson <- glmmTMB(glugea ~ site + Sex + length + confactor + (1|ecto_screener),
                        data=data,
                        ziformula=~1,
                        family=poisson)

contra_zinbinom <- glmmTMB(contra ~ site + Sex + length + confactor + (1|endo_screener),
                        data=data,
                        ziformula=~1,
                        family=nbinom1("log"))
contra_zipoisson <- glmmTMB(contra ~ site + Sex + length + confactor + (1|endo_screener),
                         data=data,
                         ziformula=~1,
                         family=poisson)

angui_zinbinom <- glmmTMB(angui ~ site + Sex + length + confactor + (1|endo_screener),
                           data=data,
                           ziformula=~1,
                           family=nbinom1("log"))
angui_zipoisson <- glmmTMB(angui ~ site + Sex + length + confactor + (1|endo_screener),
                            data=data,
                            ziformula=~1,
                            family=poisson)

# ZIGLM with pscl package
library('pscl')
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

ZI_tri <- zeroinfl(tricho ~ site + length + Sex + confactor|1, dist = "negbin", link = "logit", data=data)
Anova(ZI_tri, type='III')

ZI_tri_env <- zeroinfl(tricho ~ site + length + Sex + confactor|1, dist = "negbin", link = "logit", data=data)
Anova(ZI_tri_env, type='III')

ZI_gy <- zeroinfl(gyro ~ site + length + Sex + confactor|1, dist = "poisson", link = "logit", data=data)
Anova(ZI_gy, type='III')
overdisp_fun(ZI_gy)

ZI_glu <- zeroinfl(glugea ~ site + length + Sex + confactor|1, dist = "negbin", link = "logit", data=data)
Anova(ZI_glu, type='III')
overdisp_fun(ZI_glu)

ZI_contra <- zeroinfl(contra ~ site + length + Sex + confactor|1, dist = "poisson", link = "logit", data=data)
Anova(ZI_contra, type='III')
overdisp_fun(ZI_contra)

ZI_angui <- zeroinfl(angui ~ site + Sex + length + confactor|1, dist = "poisson", link = "logit", data=data)
Anova(ZI_angui, type='III')
boverdisp_fun(ZI_angui)


#### 1.3 Diversity ####
library(ggplot2)
g <- ggplot(data=data, aes(x=site, y=PI)) + geom_boxplot()
g + theme_minimal()

library(vegan)
library(tidyr)
data$PI_sqrt <- sqrt(data$PI)
hist(data$PI_sqrt)
data$PI2 = sign(data$PI) * abs(data$PI)^(1/3)
hist(data$PI2)

p <- ggplot(diversity, aes(x=site, y=species.richness)) + 
  geom_boxplot()
p +  theme_few() 

data_long <- gather(data, index, value, div_shannon:div_simpson, factor_key=TRUE)
ggplot(data_long, aes(x=site, y=value, fill=factor(index))) +
  geom_boxplot() + scale_fill_manual(values=c("#009E73", "#0072B2")) +  theme_few() 
summary(data)


sr <- lm(parspeciesrichness ~ length*site, data = data)

summary(sr)
Anova(sr)
plot(sr)

## diversity indices right skewed distribution


#### 1.4 DISSIMILARITIES IN PARASITE COMMUNITY COMPOSITION ####

library(vegan)
data_adonis <- read.csv("data_adonis.csv", sep=';')
data_adonis$site <- as.factor(data_adonis$site)
data_adonis$fish <- as.factor(data_adonis$site)
data_adonis$length <- as.numeric(data_adonis$length)
species.hel <- decostand(data_adonis[,9:18],method="hellinger")
Adonis1 <- adonis(species.hel ~ data_adonis$site*data_adonis$length + data_adonis$Sex, perm=999)  
Adonis1


