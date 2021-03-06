# Load libraries
library(gamm4)
library(dplyr)
library(MASS)
library(scales)
library(RColorBrewer)
library(mgcv)
library(zoo)
options(stringsAsFactors = FALSE)

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
all.f[which(all.f$Date=="2017-05-22"),"Tank"]<- NA

# Create numeric vector for number of days since start date
all.f$Days <- with(all.f, as.numeric(Date - first(Date[order(Date)])))

# Fit generalized additive model
all.f$Genotype <- factor(all.f$Genotype)

mod <- gamm4(Y ~ Genotype + s(Days, k=4, by=Genotype), random=~(1|FragID)+(1|Tank), data=all.f)
modMean <- gamm4(Y ~ s(Days, k=4), random=~(1|FragID)+(1|Tank), data=all.f)

#summary(mod$mer)

# Get fitted values
newdata <- expand.grid(Days=seq(0, max(all.f$Days)),
                       Genotype=levels(all.f$Genotype))
newdatamean <- expand.grid(Days=seq(0, max(all.f$Days)))
newdata$fit <- predict(mod$gam, newdata)
newdatamean$fit <-predict(modMean$gam, newdatamean)
#Find Derivative for mean
rMDeriv<- NA
lMDeriv<-NA
rMDeriv$Days<- newdatamean$Days+1e-6 # Right x value of point
lMDeriv$Days<- newdatamean$Days-1e-6 # Left x value of point
rMDeriv$val <- predict(modMean$gam, rMDeriv)
lMDeriv$val <- predict(modMean$gam, lMDeriv)
newdatamean$deriv <- (rMDeriv$val-lMDeriv$val)/(2e-6)

#Derivatives for all values and genotypes
rDeriv<- NA
lDeriv<-NA
rDeriv<-  newdata  
lDeriv<-  newdata
rDeriv$Days<- rDeriv$Days+1e-6 # Right x value of point
lDeriv$Days<- lDeriv$Days-1e-6 # Left x value of point
rDeriv$val <- predict(mod$gam, rDeriv)
lDeriv$val <- predict(mod$gam, lDeriv)
newdata$deriv <- (rDeriv$val-lDeriv$val)/(2e-6)



#Functions  that may work
#aTry$deriv <- predict.gam(modMean$gam, aTry, type="lpmatrix")
#aTry$ci <- predict.gam(modMean$gam, aTry, type="terms", se.fit=T)

# Simulate from multivariate normal distribution of fitted model
#For fitted model of all genos
set.seed(681) 
Rbeta <- mvrnorm(n = 1000, coef(mod$gam), vcov(mod$gam))
Xp <- predict(mod$gam, newdata = newdata, type = "lpmatrix")
rDeriv$val <- predict(mod$gam, rDeriv, type="lpmatrix")
lDeriv$val <- predict(mod$gam, lDeriv, type="lpmatrix")
Xpd <- (rDeriv$val-lDeriv$val)/(2e-6)
sim <- Xp %*% t(Rbeta)
simDeriv <- Xpd  %*% t(Rbeta)

#For mean not confinced this is necassary... for the derivs
set.seed(900)
RbetaMean <- mvrnorm(n = 1000, coef(modMean$gam), vcov(modMean$gam))
XpMean <- predict(modMean$gam, newdata = newdatamean, type = "lpmatrix")
rMDeriv$val <- predict(modMean$gam, rMDeriv, type="lpmatrix")
lMDeriv$val <- predict(modMean$gam, lMDeriv, type="lpmatrix")
Xpmd <- (rMDeriv$val-lMDeriv$val)/(2e-6)
simMean <- XpMean %*% t(RbetaMean)
simMderiv <- Xpmd  %*% t(RbetaMean)

# Extract 84% confidence intervals from simulation
newdata$lci <- apply(sim, 1, quantile, 0.08)
newdata$uci <- apply(sim, 1, quantile, 0.92)
newdata$lcideriv <- apply(simDeriv, 1, quantile, 0.08)
newdata$ucideriv <- apply(simDeriv, 1, quantile, 0.92)
newdatamean$lci <- apply(simMean, 1, quantile, 0.08)
newdatamean$uci <- apply(simMean, 1, quantile, 0.92)
newdatamean$lciDeriv <- apply(simMderiv, 1, quantile, 0.08)
newdatamean$uciDeriv <- apply(simMderiv, 1, quantile, 0.92)

# Get rid of fitted values beyond day range for which we have data

maxDays <- aggregate(all.f$Days, by=list(all.f$Genotype), FUN=max)
for(i in 1:nrow(maxDays)) {
  graphData <- newdata[!(newdata$Genotype==maxDays[i, "Group.1"] & newdata$Days>maxDays[i, "x"]),] 
}

# Plot results
library(lattice)
xyplot(fit ~ Days, groups= Genotype, data=newdata, type="l")
xyplot(deriv ~ Days, groups= Genotype, data=newdata, type="l")
xyplot(deriv~ Days | Genotype, data=newdata, type="l")
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


  #Normalize Mean
  for(day in newdatamean$Days) {
    newdatamean[which(newdatamean$Days ==day), "rfvfm"] <- newdatamean[which(newdatamean$Days==day), "fit"] / newdatamean[which(newdatamean$Days==0), "fit"]
    newdatamean[which(newdatamean$Days ==day), "rlci"] <- newdatamean[which(newdatamean$Days==day), "lci"] / newdatamean[which(newdatamean$Days==0), "fit"]
    newdatamean[which(newdatamean$Days ==day), "ruci"] <- newdatamean[which(newdatamean$Days==day), "uci"] / newdatamean[which(newdatamean$Days==0), "fit"]
  }
 
plotfitsnormal <- function(data, geno, title, key) {
  df <- droplevels(subset(data, Genotype %in% geno))
  getColor <- colorRampPalette(brewer.pal(length(geno), "Dark2"))
  colors <- getColor(length(geno))
  dfsplit <- split(df, f=df$Genotype)
  plot(NA, xlim=range(df$Days), ylim=range(df$rfvfm), xaxs="i", xlab = "Days", ylab = "Standardized Y", main=title)
  if(key==T) {
    par(new=T)
    plot(temp$Projected ~ temp$Days, type = "l", ylim=c(25, 35), xaxs="i", axes=F, lwd = 2, xlim=c(0,39), xlab="", ylab="")
    lines(temp$Thres~temp$Days, col = 'orange', lty = 2)
    lines(averageTemp$mean ~ averageTemp$Days, lwd = (.5), xlim=c(0, 39))
    lines(averageTemp$smooth ~ averageTemp$Days, col = 'red')
    axis(side=4)
    mtext("Tempearture (oC)", side=4)
    par(new=T)
    plot(NA, xlim=range(df$Days), ylim=range(df$rfvfm), xaxs="i", xlab = "Days", ylab = "Standardized Y", main=title)
  }
  for (i in 1:length(dfsplit)) {
    with(dfsplit[[i]], {
      lines(Days, rfvfm)
      addpoly(Days, rlci, ruci, col=alpha(colors[i], 0.3))
    })
  }
 
  legend("bottom", legend=levels(df$Genotype), pch=15, bty="n",
         col=alpha(colors, 0.4), pt.cex=1.5,
         y.intersp=0.3)
}

plotfitsnormal(newdata, c("1725", "1739", "1736", "1737", "1730", "1745", "1738", "1722","1735", "1729"), "", T)
lines(rfvfm~Days, data=newdatamean, lwd=2, col="red")
addpoly(newdatamean$Days, newdatamean$rlci, newdatamean$ruci, col=alpha("red", 0.3))
library(lattice)
xyplot(rfvfm ~ Days, groups= Genotype, data=newdata, type="l")
png(filename = "Output/Figures/normalsplit.png", width=9, height = 9, units = "in", res = 300)
xyplot(rfvfm ~ Days | Genotype, data=newdata, type="l")
dev.off()
plotfitsnormal(graphData, c("1722", "1732", "1738", "1735", "1729"), "Hello", T)
png(filename = "Output/Figures/normalall.png", width=9, height = 9, units = "in", res = 300)
plotfitsnormal(graphData, c("1721", "1722", "1723", "1725", "1726", "1727", "1728", "1729", "1730", "1731", "1732", "1733", "1734", "1735", "1736", "1737", "1738", "1739", "1741", "1743", "1745", "1747", "1748", "1750","1751", "1752", "1753", "1755", "1757", "1758", "1759"), "All Genotypes", T)
#lines(rfvfm~Days, data=newdatamean, lwd=2, col="red")
#addpoly(newdatamean$Days, newdatamean$rlci, newdatamean$ruci, col=alpha("red", 0.3))
abline(v=17, col= "red")
text(x= 21.4, y = 1.0, adj = c(.5), "32oC Reached")
dev.off()
abline(v=17, col= "red")
text(x= 21.4, y = 1.0, adj = c(.5), "32oC Reached")
#Residuals OUtlier
#residualsModel <- residuals(mod$mer)
#sdmod <-sd(residuals(mod$mer))
#outliersMod <- residualsModel[residualsModel<= 2*sdmod]
#all.ff <- all.f[which(residualsModel < 2 * sdmod), ]
#mod.ff <- gamm4(Y ~ Genotype + s(Days, k=4, by=Genotype), random=~(1|FragID), data=all.ff)
xyplot(Y~Date, groups=Genotype, data=all.ff)
xyplot(Y ~ Days | Genotype, data=all.ff)
xyplot(Y ~ Days | Genotype, data=all.f)

#Stats on NormalPlot
library(multcomp)
library(multcompView)
library(lsmeans)

#Day 0 initial Values
dayZ <- newdata[which(newdata$Days==0),]
dayZ <- dayZ[order(dayZ$Genotype),]

#Two Days Post Heat Stress Day 22

dayTT <- newdata[which(newdata$Days=="22"),]
dayTT <-dayTT[order(dayTT$rfvfm),]
dayTT$score
for(i in 1:nrow(dayTT)) {
  dayTT[i, "score"]<-i
}

#Determine Overlap between confidence intervals to genotypes at Day 22 using normal
compareGeno<-data.frame(NA)
dayTT <-dayTT[order(dayTT$Genotype),]
for(i in 1:nrow(dayTT)){
  compareGeno[i, "Genotype"]<-dayTT[i, "Genotype"]
  compareGeno[i, "Initial"]<-dayZ[i, "fit"]
    if(dayTT[i, "rlci"]>newdatamean[which(newdatamean$Days==22), "ruci"]) {
      compareGeno[i, "Normalized Result Day 22"]<-"Above Average"
    }
    else if(dayTT[i, "ruci"]<newdatamean[which(newdatamean$Days==22), "rlci"]){
      compareGeno[i, "Normalized Result Day 22"]<-"Below Average"
    }
    else{
      compareGeno[i, "Normalized Result Day 22"]<-"Not Signifcantly Different"
    }
}
compareGeno <- compareGeno[, 2:4]

#Determine Overlap of non normalized

for(i in 1:nrow(dayTT)){
  
  if(dayTT[i, "lci"]>newdatamean[which(newdatamean$Days==22), "uci"]) {
    compareGeno[i, "Result Day 22"]<-"Above Average"
  }
  else if(dayTT[i, "uci"]<newdatamean[which(newdatamean$Days==22), "lci"]){
    compareGeno[i, "Result Day 22"]<-"Below Average"
  }
  else{
    compareGeno[i, "Result Day 22"]<-"Not Signifcantly Different"
  }
}

for(i in 1:nrow(dayTT)){
  
  if(dayTT[i, "lcideriv"]>newdatamean[which(newdatamean$Days==22), "uciDeriv"]) {
    compareGeno[i, "Slope Day 22"]<-"Above Average"
  }
  else if(dayTT[i, "ucideriv"]<newdatamean[which(newdatamean$Days==22), "lciDeriv"]){
    compareGeno[i, "Slope Day 22"]<-"Below Average"
  }
  else{
    compareGeno[i, "Slope Day 22"]<-"Not Signifcantly Different"
  }
}

#five <-five[order(five$Genotype),]
#scorefive <-aggregate(five$score, by = list(five$Genotype), FUN = paste)

#Day 27
dayTS <- newdata[which(newdata$Days=="27"),]
dayTS <-dayTS[order(dayTS$rfvfm),]
dayTS$score
for(i in 1:nrow(dayTS)) {
  dayTS[i, "score"]<-i
}

#scoreten <-aggregate(ten$score, by = list(ten$Genotype), FUN = paste)
#dayTS <-dayTS[order(dayTS$Genotype),]
#for(i in 1:nrow(dayTS)){
#    if(dayTS[i, "rlci"]>newdatamean[which(newdatamean$Days==27), "ruci"]) {
#    compareGeno[i, "Normalized Result Day 27"]<-"Above Average"
#  }
#  else if(dayTS[i, "ruci"]<newdatamean[which(newdatamean$Days==27), "rlci"]){
#    compareGeno[i, "Normalized Result Day 27"]<-"Below Average"
#  }
#  else{
#    compareGeno[i, "Normalized Result Day 27"]<-"Not Signifcantly Different"
#  }
#}

#for(i in 1:nrow(dayTS)){
  
#  if(dayTS[i, "lci"]>newdatamean[which(newdatamean$Days==27), "uci"]) {
#    compareGeno[i, "Result Day 27"]<-"Above Average"
#  }
#  else if(dayTS[i, "uci"]<newdatamean[which(newdatamean$Days==27), "lci"]){
#   compareGeno[i, "Result Day 27"]<-"Below Average"
#  }
#  else{
#   compareGeno[i, "Result Day 27"]<-"Not Signifcantly Different"
#  }
#}

# Six Days post heat Day 26
dayTwS <- newdata[which(newdata$Days=="26"),]
dayTwS <-dayTwS[order(dayTwS$rfvfm),]
for(i in 1:nrow(dayTwS)) {
  dayTwS[i, "score"]<-i
}
dayTwS <-dayTwS[order(dayTwS$Genotype),]

dayTwS <-dayTwS[order(dayTwS$Genotype),]
for(i in 1:nrow(dayTwS)){
  if(dayTwS[i, "rlci"]>newdatamean[which(newdatamean$Days==26), "ruci"]) {
    compareGeno[i, "Normalized Result Day 26"]<-"Above Average"
  }
  else if(dayTwS[i, "ruci"]<newdatamean[which(newdatamean$Days==26), "rlci"]){
    compareGeno[i, "Normalized Result Day 26"]<-"Below Average"
  }
  else{
    compareGeno[i, "Normalized Result Day 26"]<-"Not Signifcantly Different"
  }
}

for(i in 1:nrow(dayTwS)){
  
  if(dayTwS[i, "lci"]>newdatamean[which(newdatamean$Days==26), "uci"]) {
    compareGeno[i, "Result Day 26"]<-"Above Average"
  }
  else if(dayTwS[i, "uci"]<newdatamean[which(newdatamean$Days==26), "lci"]){
    compareGeno[i, "Result Day 26"]<-"Below Average"
  }
  else{
    compareGeno[i, "Result Day 26"]<-"Not Signifcantly Different"
  }
}

for(i in 1:nrow(dayTwS)){
  
  if(dayTwS[i, "lcideriv"]>newdatamean[which(newdatamean$Days==26), "uciDeriv"]) {
    compareGeno[i, "Slope Day 26"]<-"Above Average"
  }
  else if(dayTwS[i, "ucideriv"]<newdatamean[which(newdatamean$Days==26), "lciDeriv"]){
    compareGeno[i, "Slope Day 26"]<-"Below Average"
  }
  else{
    compareGeno[i, "Slope Day 26"]<-"Not Signifcantly Different"
  }
}

plotfitsnormal(newdata, c("1729", "1735", "1739", "1741", "1751", "1759", "1722", "1730","1731", "1721", "1738", "1743", "1747", "1750", "1755"), "", T)
lines(rfvfm~Days, data=newdatamean, lwd=2, col="red")
addpoly(newdatamean$Days, newdatamean$rlci, newdatamean$ruci, col=alpha("red", 0.3))
abline(v=17, col= "red")
text(x= 21.4, y = 1.0, adj = c(.5), "32oC Reached")
abline(v=32, col= "red")
abline(v=26, col= "red")
text(x= 32.4, y = 1.0, adj = c(.5), "Statistically Different")

#Nine Days Post Heat Stress 29
dayTN <- newdata[which(newdata$Days=="29"),]
dayTN <-dayTN[order(dayTN$rfvfm),]
for(i in 1:nrow(dayTN)) {
  dayTN[i, "score"]<-i
}
dayTN <-dayTN[order(dayTN$Genotype),]

for(i in 1:nrow(dayTN)){
  if(dayTN[i, "rlci"]>newdatamean[which(newdatamean$Days==29), "ruci"]) {
    compareGeno[i, "Normalized Result Day 29"]<-"Above Average"
  }
  else if(dayTN[i, "ruci"]<newdatamean[which(newdatamean$Days==29), "rlci"]){
    compareGeno[i, "Normalized Result Day 29"]<-"Below Average"
  }
  else{
    compareGeno[i, "Normalized Result Day 29"]<-"Not Signifcantly Different"
  }
}

for(i in 1:nrow(dayTN)){
  
  if(dayTN[i, "lci"]>newdatamean[which(newdatamean$Days==29), "uci"]) {
    compareGeno[i, "Result Day 29"]<-"Above Average"
  }
  else if(dayTN[i, "uci"]<newdatamean[which(newdatamean$Days==29), "lci"]){
    compareGeno[i, "Result Day 29"]<-"Below Average"
  }
  else{
    compareGeno[i, "Result Day 29"]<-"Not Signifcantly Different"
  }
}

for(i in 1:nrow(dayTN)){
  
  if(dayTN[i, "lcideriv"]>newdatamean[which(newdatamean$Days==29), "uciDeriv"]) {
    compareGeno[i, "Slope Day 29"]<-"Above Average"
  }
  else if(dayTN[i, "ucideriv"]<newdatamean[which(newdatamean$Days==29), "lciDeriv"]){
    compareGeno[i, "Slope Day 29"]<-"Below Average"
  }
  else{
    compareGeno[i, "Slope Day 29"]<-"Not Signifcantly Different"
  }
}
plotfitsnormal(newdata, c("1759", "1735", "1738"))
plotfitsnormal(newdata, c("1729", "1735", "1739", "1741", "1751", "1759", "1722", "1730","1731", "1721", "1738", "1743", "1747", "1750", "1755"), "", T)
lines(rfvfm~Days, data=newdatamean, lwd=2, col="red")
addpoly(newdatamean$Days, newdatamean$rlci, newdatamean$ruci, col=alpha("red", 0.3))
abline(v=17, col= "red")
text(x= 21.4, y = 1.0, adj = c(.5), "32oC Reached")
abline(v=32, col= "red")

#twelve Post Heat Stress Day 32
dayThT <- newdata[which(newdata$Days==32),]
dayThT <-dayThT[order(dayThT$rfvfm),]
for(i in 1:nrow(dayThT)) {
  dayThT[i, "score"]<-i
}
dayThT <-dayThT[order(dayThT$Genotype),]


for(i in 1:nrow(dayThT)){
  if(dayThT[i, "rlci"]>newdatamean[which(newdatamean$Days==32), "ruci"]) {
    compareGeno[i, "Normalized Result Day 32"]<-"Above Average"
  }
  else if(dayThT[i, "ruci"]<newdatamean[which(newdatamean$Days==32), "rlci"]){
    compareGeno[i, "Normalized Result Day 32"]<-"Below Average"
  }
  else{
    compareGeno[i, "Normalized Result Day 32"]<-"Not Signifcantly Different"
  }
}



for(i in 1:nrow(dayThT)){
  
  if(dayThT[i, "lci"]>newdatamean[which(newdatamean$Days==32), "uci"]) {
    compareGeno[i, "Result Day 32"]<-"Above Average"
  }
  else if(dayThT[i, "uci"]<newdatamean[which(newdatamean$Days==32), "lci"]){
    compareGeno[i, "Result Day 32"]<-"Below Average"
  }
  else{
    compareGeno[i, "Result Day 32"]<-"Not Signifcantly Different"
  }
}

for(i in 1:nrow(dayThT)){
  
  if(dayThT[i, "lcideriv"]>newdatamean[which(newdatamean$Days==32), "uciDeriv"]) {
    compareGeno[i, "Slope Day 32"]<-"Above Average"
  }
  else if(dayThT[i, "ucideriv"]<newdatamean[which(newdatamean$Days==32), "lciDeriv"]){
    compareGeno[i, "Slope Day 32"]<-"Below Average"
  }
  else{
    compareGeno[i, "Slope Day 32"]<-"Not Signifcantly Different"
  }
}
#png(filename = "Output/Figures/bestworstday15.png", width=7, height = 7, units = "in", res = 300)
#dev.off()
plotfitsnormal(newdata, c("1729", "1735", "1739", "1741", "1751", "1753", "1759", "1736", "1721", "1722", "1731", "1738", "1747", "1750", "1755", "1743"), "", F)
lines(rfvfm~Days, data=newdatamean, lwd=2, col="red")
addpoly(newdatamean$Days, newdatamean$rlci, newdatamean$ruci, col=alpha("red", 0.3))
abline(v=20, col= "red")
text(x= 21.4, y = 1.0, pos =2, "32oC Reached")
abline(v=32, col= "red")
text(x= 36.4, y = 1.0, pos=2, "Statistically Different @ Day 12")

dev.off()
#Fifteen Days Post Heat Stress Day 35
#Ls Trends using mod mer modgam
#Overall FvFm
#Initial FvFm vs. Performance

dayThF <- newdata[which(newdata$Days==35),]
dayThF <-dayThF[order(dayThF$rfvfm),]
for(i in 1:nrow(dayThF)) {
  dayThF[i, "score"]<-i
}
dayThF <-dayThF[order(dayThF$Genotype),]

for(i in 1:nrow(dayThF)){
  if(dayThF[i, "rlci"]>newdatamean[which(newdatamean$Days==35), "ruci"]) {
    compareGeno[i, "Normalized Result Day 35"]<-"Above Average"
  }
  else if(dayThF[i, "ruci"]<newdatamean[which(newdatamean$Days==35), "rlci"]){
    compareGeno[i, "Normalized Result Day 35"]<-"Below Average"
  }
  else{
    compareGeno[i, "Normalized Result Day 35"]<-"Not Signifcantly Different"
  }
}

for(i in 1:nrow(dayThF)){
  
  if(dayThF[i, "lci"]>newdatamean[which(newdatamean$Days==35), "uci"]) {
    compareGeno[i, "Result Day 35"]<-"Above Average"
  }
  else if(dayThF[i, "uci"]<newdatamean[which(newdatamean$Days==35), "lci"]){
    compareGeno[i, "Result Day 35"]<-"Below Average"
  }
  else{
    compareGeno[i, "Result Day 35"]<-"Not Signifcantly Different"
  }
}

for(i in 1:nrow(dayThF)){
  
  if(dayThF[i, "lcideriv"]>newdatamean[which(newdatamean$Days==35), "uciDeriv"]) {
    compareGeno[i, "Slope Day 35"]<-"Above Average"
  }
  else if(dayThF[i, "ucideriv"]<newdatamean[which(newdatamean$Days==35), "lciDeriv"]){
    compareGeno[i, "Slope Day 35"]<-"Below Average"
  }
  else{
    compareGeno[i, "Slope Day 35"]<-"Not Signifcantly Different"
  }
}

plotfitsnormal(graphData, c("1729", "1735", "1739", "1741", "1751", "1753", "1759", "1723", "1736", "1745", "1728", "1732", "1734","1757","1758", "1721", "1722", "1731", "1738", "1747", "1750", "1755"), "", T)
lines(rfvfm~Days, data=newdatamean, lwd=2, col="red")
addpoly(newdatamean$Days, newdatamean$rlci, newdatamean$ruci, col=alpha("red", 0.3))
abline(v=20, col= "red")
text(x= 21.4, y = 1.0, pos =2, "32oC Reached")
abline(v=35, col= "red")
text(x= 36.4, y = 1.0, pos=2, "Statistically Different @ Day 18")

#eighteen days post heat stress Day 38
dayThE <- newdata[which(newdata$Days==38),]
dayThE <-dayThE[order(dayThE$rfvfm),]
for(i in 1:nrow(dayThE)) {
  dayThE[i, "score"]<-i
}
dayThE <-dayThE[order(dayThE$Genotype),]

for(i in 1:nrow(dayThE)){
  if(dayThE[i, "rlci"]>newdatamean[which(newdatamean$Days==38), "ruci"]) {
    compareGeno[i, "Normalized Result Day 38"]<-"Above Average"
  }
  else if(dayThE[i, "ruci"]<newdatamean[which(newdatamean$Days==38), "rlci"]){
    compareGeno[i, "Normalized Result Day 38"]<-"Below Average"
  }
  else{
    compareGeno[i, "Normalized Result Day 38"]<-"Not Signifcantly Different"
  }
}

for(i in 1:nrow(dayThE)){
  
  if(dayThE[i, "lci"]>newdatamean[which(newdatamean$Days==38), "uci"]) {
    compareGeno[i, "Result Day 38"]<-"Above Average"
  }
  else if(dayThE[i, "uci"]<newdatamean[which(newdatamean$Days==38), "lci"]){
    compareGeno[i, "Result Day 38"]<-"Below Average"
  }
  else{
    compareGeno[i, "Result Day 38"]<-"Not Signifcantly Different"
  }
}

for(i in 1:nrow(dayThE)){
  
  if(dayThE[i, "lcideriv"]>newdatamean[which(newdatamean$Days==38), "uciDeriv"]) {
    compareGeno[i, "Slope Day 38"]<-"Above Average"
  }
  else if(dayThE[i, "ucideriv"]<newdatamean[which(newdatamean$Days==38), "lciDeriv"]){
    compareGeno[i, "Slope Day 38"]<-"Below Average"
  }
  else{
    compareGeno[i, "Slope Day 38"]<-"Not Signifcantly Different"
  }
}

total<-NA
total$Genotype<-nine$Genotype
total$Score<-nine$score+twelve$score+fifteen$score+teen$score+tweone$score
total<-total[,2:3]
View(total)

View(aggregate(total, by = list(twenty$Genotype), FUN = paste))
totalTable <-cbind(geno3table$Genotype, total)
View(totalTable)
boxplot(total)
histogram(total)

