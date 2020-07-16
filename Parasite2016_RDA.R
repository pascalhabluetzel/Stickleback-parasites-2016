library(vegan)
library(ggvegan)
library(ggthemes)
library(raster)

setwd('C:/Users/u0113095/Google Drive/PhD/2 Parasite/2016/Parasite2016_analysis')

ab <- read.csv("abundance.csv", sep=";")
ab <- ab[,-1]
pre <- read.csv("prevalence.csv", sep=';')

spavar <- read.csv("space2.csv", sep=';') #spatial variables
spavar <- spavar[,-1]
spa <- read.csv("space.csv", sep=';')
plot(spa$long ~ spa$lat)
dist <- pointDistance(spa, allpairs = TRUE, lonlat = FALSE)

env <- read.csv("env_max.csv", sep=';')

# Hellinger transformation of the species dataset ####
spe.hel <- decostand(pre,"hellinger") #prevalence or abundance

#### Environmental variables ####
spe.rda <- rda(spe.hel, env)
plot(spe.rda, scaling = 1)
summary(spe.rda)
anova(spe.rda)
anova.cca(spe.rda, step=1000);
RsquareAdj(spe.rda)$adj.r.squared;
RsquareAdj(spe.rda)$r.squared

mod0 <- rda(spe.hel ~ 1, env)  # Model with intercept only
mod1 <- rda(spe.hel ~ ., env)  # Model with all explanatory variables
step.res <- ordiR2step(mod0, mod1, direction = "backward",perm.max = 200)
step.res$anova  # Summary table

g <- autoplot(spe.rda, arrows = FALSE, geom = c("point", "text")) + geom_text(size = 10) + theme_few() + xlim(c(-0.7,1.2)) + scale_color_manual(values = c("black","#009999","#006666"))
g + geom_hline(yintercept = 0, linetype="dotted") + geom_vline(xintercept = 0, linetype="dotted")
ggsave("environmentspecies.png", units="in", width=10, height=10, dpi=300)

# Spatial variables  ####
distance_matrix <- read.csv("distance_matrix.csv", sep=';')
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

spa.PCNM <- as.data.frame(xy.PCoA$points[1:nrow(spa),1:2])

# Distance matrix  ####
spa.rda <- rda(spe.hel,spa.PCNM)
summary(spa.rda)
plot(spa.rda)
anova(spa.rda)
anova.cca(spa.rda, step=1000)
R2 <- RsquareAdj(spa.rda)$r.squared;
R2 <- RsquareAdj(spa.rda)$adj.r.squared;
R2

mod0 <- rda(spe.hel ~ 1, spa2)  # Model with intercept only
mod1 <- rda(spe.hel ~ ., spa2)  # Model with all explanatory variables
step.res <- ordiR2step(mod0, mod1, direction = "backward",perm.max = 200)
step.res$anova  # Summary table


g <- autoplot(spa.rda, arrows = FALSE, geom = c("point", "text")) + geom_text(size = 16) + theme_few() + scale_color_manual(values = c("black","darkgrey","black"))
g + geom_hline(yintercept = 0, linetype="dotted") + geom_vline(xintercept = 0, linetype="dotted") 
ggsave("environmentspace2.png", units="in", width=10, height=10, dpi=300)

# Spatial variables  ####
spa.rda <- rda(env,space)
summary(spa.rda)
plot(spa.rda)
anova(spa.rda)
anova.cca(spa.rda, step=1000)
R2 <- RsquareAdj(spa.rda)$r.squared;
R2 <- RsquareAdj(spa.rda)$adj.r.squared;
R2

mod0 <- rda(spe.hel ~ 1, space)  # Model with intercept only
mod1 <- rda(spe.hel ~ ., space)  # Model with all explanatory variables
step.res <- ordiR2step(mod0, mod1, direction = "backward",perm.max = 200)
step.res$anova  # Summary table


g <- autoplot(spa.rda, arrows = FALSE, geom = c("point", "text")) + geom_text(size = 16) + theme_few() + scale_color_manual(values = c("black","darkgrey","black"))
g + geom_hline(yintercept = 0, linetype="dotted") + geom_vline(xintercept = 0, linetype="dotted") 
ggsave("environmentspace2.png", units="in", width=10, height=10, dpi=300)

#variance partitioning
spe.varpart1 <- varpart(spe.hel,space, env)
par(mfrow=c(1,2))
showvarparts(2)
plot(spe.varpart1,digits=2)
spe.varpart1



