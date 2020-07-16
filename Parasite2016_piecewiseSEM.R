library(devtools)
install_github("jslefche/piecewiseSEM@devel", build_vignette = TRUE)
library(piecewiseSEM)
library(nlme)
library(lme4)
library(nlme)

setwd('C:/Users/u0113095/Google Drive/PhD/2 Parasite/2016/Analysis_2020/data')
data <- read.csv("data_2016_SEM.csv", sep=';')
data$site <- as.factor(data$site)
data$fish <- as.factor(data$fish)
data$length <- as.numeric(data$length)
data$updist2 <- as.numeric(data$updist2)

data$PI_ecto2 = sign(data$PI_ecto) * abs(data$PI_ecto)^(1/3)
data$PI_endo2 = sign(data$PI_endo) * abs(data$PI_endo)^(1/3)
data$PI2 = sign(data$PI) * abs(data$PI)^(1/3)

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


