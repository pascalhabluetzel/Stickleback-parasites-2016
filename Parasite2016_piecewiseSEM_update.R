#Which version of R am I supposed to use here?

setwd('C:/Users/u0113095/Google Drive/PhD/2 Parasite/2016/Analysis_2020/Github')
setwd('C:/Users/pascalh/Documents/GitHub/Stickleback-parasites-2016')

#install.packages("magrittr")
library(magrittr)
install.packages("dplyr")
install.packages("ellipsis")
install.packages("vctrs")
library(dplyr)

env <- read.csv("Environment_R.csv", sep=';')
spavar <- read.csv("space2.csv", sep=';') 
spavar_exp <- spavar %>% slice(rep(1:n(), table(as.factor(data2$site)))) #Code is not reproducible. What is "data2"?
data2 <- read.csv('data_update.csv', sep=';')
data3 <- cbind(spavar_exp, data2)
data3$site <- as.factor(data3$site) #Code is not reproducible. What is "data3"?
data3 <- data3[,-1]

model <- psem(
  lme(PI_ecto ~ confactor + length + T_av +  con_av + COD_av + NH4_av + 
        Nt_av +  meander +  netcen, random=~1|site, data=data3),
  lme(confactor ~ T_av +  con_av + COD_av + NH4_av + 
        Nt_av +  meander + netcen, random=~1|site, data=data3),
  lme(length~ confactor + T_av +  con_av + COD_av + NH4_av + 
        Nt_av +  meander + netcen, random=~1|site, data=data3))
summary(model)

model <- lme(PI_ecto ~ confactor + length + T_av + O2_sat_av + con_av + COD_av + NH4_av + 
      Nt_av + poolriffle + meander + netcen + updist, random=~1|site, data=data3)


library(devtools)
install_github('sinhrks/ggfortify')
library(ggfortify); library(ggplot2)

pcaenv <- prcomp(data3[,c(12:17)])
summary(pca)
screeplot(pca)
png(file="Results/PCA.png", res=600, width=3000, height=3000)
autoplot(pcaenv, loadings.label = TRUE, loadings.label.size = 4, loadings.label.colour ='black',
         loadings.colour = 'black') + theme_bw()
dev.off()
pca <- pca$x[, 1]

data4 <- cbind(data3, pca)
data4$poolriffle <- as.factor(data4$poolriffle)
library(piecewiseSEM)
library(MASS)
library(car)
library(nlme)

model <- psem(
  lme(PI ~ confactor + length + pca + poolriffle + meander + updist + netcen, random=~1|site, data=data4),
  lme(confactor ~ length + pca  + poolriffle + meander + updist + netcen, random=~1|site, data=data4),
  lme(length~ pca  + poolriffle + meander + updist + netcen, random=~1|site, data=data4))
summary(model)
