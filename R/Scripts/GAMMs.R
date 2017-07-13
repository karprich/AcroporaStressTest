# Load libraries
library(gamm4)
library(dplyr)
library(MASS)
library(scales)
library(RColorBrewer)

# Define function to use for plotting confidence intervals
addpoly <- function(x, y1, y2, col=alpha("lightgrey", 0.8), ...){
  ii <- order(x)
  y1 <- y1[ii]
  y2 <- y2[ii]
  x <- x[ii]
  polygon(c(x, rev(x)), c(y1, rev(y2)), col=col, border=NA, ...)
}

# Load data
load("Output/all.f.RData")
all.f <- na.omit(all.f)
all.f <- droplevels(all.f)

# Create numeric vector for number of days since start date
all.f$Days <- with(all.f, as.numeric(Date - first(Date[order(Date)])))

# Fit generalized additive model
all.f$Genotype <- factor(all.f$Genotype)

mod <- gamm4(Y ~ Genotype + s(Days, k=4, by=Genotype), random=~(1|FragID), data=all.f)

#summary(mod$mer)

# Get fitted values
newdata <- expand.grid(Days=seq(0, max(all.f$Days)),
                       Genotype=levels(all.f$Genotype))
newdata$fit <- predict(mod$gam, newdata)

# Simulate from multivariate normal distribution of fitted model
set.seed(789)
Rbeta <- mvrnorm(n = 1000, coef(mod$gam), vcov(mod$gam))
Xp <- predict(mod$gam, newdata = newdata, type = "lpmatrix")
sim <- Xp %*% t(Rbeta)

# Extract 84% confidence intervals from simulation
newdata$lci <- apply(sim, 1, quantile, 0.08)
newdata$uci <- apply(sim, 1, quantile, 0.92)

# Get rid of fitted values beyond day range for which we have data
maxDays <- aggregate(all.f$Days, by=list(all.f$Genotype), FUN=max)
for(i in 1:nrow(maxDays)) {
  newdata <- newdata[!(newdata$Genotype==maxDays[i, "Group.1"] & newdata$Days>maxDays[i, "x"]),] 
}

# Plot results
library(lattice)
xyplot(fit ~ Days, groups= Genotype, data=newdata, type="l")
xyplot(fit ~ Days | Genotype, data=newdata, type="l")

# Plot single genotype with confidence interval
levels(newdata$Genotype)
G1721 <- subset(newdata, Genotype=="1721")

# Define function to plot fit and CI for specified genotypes
plotfits <- function(data, geno) {
  df <- droplevels(subset(data, Genotype %in% geno))
  getColor <- colorRampPalette(brewer.pal(length(geno), "Dark2"))
  colors <- getColor(length(geno))
  dfsplit <- split(df, f=df$Genotype)
  plot(NA, xlim=range(df$Days), ylim=range(df$fit), xaxs="i", xlab = "Days", ylab = "Y")
  for (i in 1:length(dfsplit)) {
    with(dfsplit[[i]], {
      lines(Days, fit)
      addpoly(Days, lci, uci, col=alpha(colors[i], 0.3))
    })
  }
  legend("bottomleft", legend=levels(df$Genotype), pch=15, bty="n",
         col=alpha(colors, 0.4), pt.cex=1.5,
         y.intersp=0.3)
}

# Plot...
par(mfrow = c(1,1))
plotfits(newdata, "1721")

abline(v=17, col = "red")
plotfits(newdata, "1722")
abline(v=17, col = "red")
plotfits(newdata, "1734")
abline(v=17, col = "red")

plotfits(newdata, c("1721", "1755", "1738"))
plotfits(newdata, c("1755", "1721", "1738"))

plotfits(newdata, c("1721", "1755", "1738", "1731", "1750"))
#Kelsey Genotypes
plotfits(newdata, c("1725", "1739", "1736", "1737", "1745"))
#Cooper
plotfits(newdata, c("1730"))
#Elkhorn
plotfits(newdata, c("1745"))
#Best
plotfits(newdata, c("1735", "1729", "1751", "1759"))
#Worst
plotfits(newdata, c("1732", "1734", "1722", "1738","1730"))
#Plot w/ Cooper, Kelsey, Elkhorn, Best, Worst
png(filename = "Output/Figures/diseasegenotypes.png", width=9, height = 9, units = "in", res = 300)
plotfits(newdata, c("1725", "1739", "1736", "1737", "1730", "1745", "1738", "1722","1735", "1729"))
legend("bottomleft", inset=c(.1, 0), legend=c("Worst-1", "Kelsey-1", "Best-1", "Cooper", "Best-2", "Kelsey-2", "Kelsey-3", "Worst-2", "Kelsey-4", "Elkhorn"), bty="n",
       y.intersp=0.3)
dev.off()

plotfits(newdata, c("1735", "1758", "1755", "1736", "1738", "1731", "1733", "1721", "1729"))
plotfits(newdata, c("1721", "1722", "1723", "1725", "1726", "1727", "1728", "1729", "1730", "1731", "1732", "1733", "1734", "1735", "1736", "1737", "1738", "1739", "1741", "1743", "1745", "1747", "1748", "1750","1751", "1752", "1753", "1755", "1757", "1758", "1759"))
abline(v=17, col = "red")
abline(h= 0)
plotfits(newdata, levels(all.f$Genotype))
# Panel xyplot()
png(filename = "Output/Figures/first14ipamscurve.png", width=7, height = 7, units = "in", res = 300)
xyplot(fit ~ Days | Genotype, data = newdata, type = "l")
dev.off()

? xyplot
# Plot Genotypes seperately with confident intervals and points

#Normalize data
for(geno in newdata$Genotype) {
  sdata <- newdata[which(newdata$Genotype==geno),]
  for(day in sdata$Days) {
    sdata[which(sdata$Days ==day), "rfvfm"] <- sdata[which(sdata$Days==day), "fit"] / sdata[which(sdata$Days==0), "fit"]
    sdata[which(sdata$Days ==day), "rlci"] <- sdata[which(sdata$Days==day), "lci"] / sdata[which(sdata$Days==0), "fit"]
    sdata[which(sdata$Days ==day), "ruci"] <- sdata[which(sdata$Days==day), "uci"] / sdata[which(sdata$Days==0), "fit"]
  }
  newdata[which(newdata$Genotype==geno), "rfvfm"] <- sdata[, "rfvfm"]
  newdata[which(newdata$Genotype==geno), "rlci"] <- sdata[, "rlci"]
  newdata[which(newdata$Genotype==geno), "ruci"] <- sdata[, "ruci"]
}

plotfitsnormal <- function(data, geno) {
  df <- droplevels(subset(data, Genotype %in% geno))
  getColor <- colorRampPalette(brewer.pal(length(geno), "Dark2"))
  colors <- getColor(length(geno))
  dfsplit <- split(df, f=df$Genotype)
  plot(NA, xlim=range(df$Days), ylim=range(df$rfvfm), xaxs="i", xlab = "Days", ylab = "Y")
  for (i in 1:length(dfsplit)) {
    with(dfsplit[[i]], {
      lines(Days, rfvfm)
      addpoly(Days, rlci, ruci, col=alpha(colors[i], 0.3))
    })
  }
  legend("bottomleft", legend=levels(df$Genotype), pch=15, bty="n",
         col=alpha(colors, 0.4), pt.cex=1.5,
         y.intersp=0.3)
}

plotfitsnormal(newdata, c("1725", "1739", "1736", "1737", "1730", "1745", "1738", "1722","1735", "1729"))

library(lattice)
xyplot(rfvfm ~ Days, groups= Genotype, data=newdata, type="l")
xyplot(rfvfm ~ Days | Genotype, data=newdata, type="l")
plotfitsnormal(newdata, c("1722", "1732", "1738", "1735"))
#Residuals OUtlier
residualsModel <- residuals(mod$mer)
sdmod <-sd(residuals(mod$mer))
outliersMod <- residualsModel[residualsModel<= 2*sdmod]
all.ff <- all.f[which(residualsModel < 2 * sdmod), ]
mod.ff <- gamm4(Y ~ Genotype + s(Days, k=4, by=Genotype), random=~(1|FragID), data=all.ff)
xyplot(Y~Date, groups=Genotype, data=all.ff)
xyplot(Y ~ Days | Genotype, data=all.ff)
xyplot(Y ~ Days | Genotype, data=all.f)
