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
all.f <- droplevels(all.f)

# Create numeric vector for number of days since start date
all.f$Days <- with(all.f, as.numeric(Date - first(Date[order(Date)])))

# Fit generalized additive model
all.f$Genotype <- factor(all.f$Genotype)

mod <- gamm4(Y ~ Genotype + s(Days, k=4, by=Genotype), random=~(1|FragID), data=all.f)

summary(mod$mer)

# Get fitted values
newdata <- expand.grid(Days=seq_len(max(all.f$Days)),
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
  plot(NA, xlim=range(df$Days), ylim=c(.2, .7), xaxs="i")
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
plotfits(newdata, c("1721", "1755", "1738", "1731", "1750"))
plotfits(newdata, c("1735", "1758", "1755", "1736", "1738", "1731", "1733", "1721"))
abline(v=17, col = "red")
plotfits(newdata, levels(all.f$Genotype))
# Panel xyplot()
xyplot(fit ~ Days | Genotype, data = newdata, type = "l")

? xyplot
