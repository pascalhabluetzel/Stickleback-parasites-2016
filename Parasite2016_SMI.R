setwd('C:/Users/u0113095/Google Drive/PhD/2 Parasite/2016/Analysis_2020/Github/data')
data <- read.csv("data_2016_1211.csv", sep=';')
library(ggplot2)
library(lmodel2)
model <- lm(weight ~ length, data = data)
model2 <- lm(weight ~ poly(length, 2, raw = TRUE), data = data)
model3 <- lm(weight ~ poly(length, 3, raw = TRUE), data = data) 
ggplot(data, aes(length, weight) ) +
  geom_point() +
  stat_smooth(method = lm, formula = y ~ poly(x, 2, raw = TRUE)) + theme_classic()
model <- lmodel2(log1p(weight)~log1p(length), data=data, nperm=999)
plot(model, "SMA")
model
