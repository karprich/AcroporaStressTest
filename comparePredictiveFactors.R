#Check for correlations among genotypes

# Load libraries
library(gamm4)
library(plyr)
library(dplyr)
library(MASS)
library(scales)
library(RColorBrewer)
library(mgcv)
library(zoo)
library(multcomp)
library(multcompView)
library(lsmeans)
library(PMCMR)
options(stringsAsFactors = FALSE)

# Define function to use for plotting confidence intervals
addpoly <- function(x, y1, y2, col=alpha("lightgrey", 0.8), ...){
  ii <- order(x)
  y1 <- y1[ii]
  y2 <- y2[ii]
  x <- x[ii]
  polygon(c(x, rev(x)), c(y1, rev(y2)), col=col, border=NA, ...)
}
load("Output/all.f.Rdata")
all.f <- na.omit(all.f)
all.f <- droplevels(all.f)
all.f[which(all.f$Date=="2017-05-22"),"Tank"]<- NA

# Create numeric vector for number of days since start date
all.f$Days <- with(all.f, as.numeric(Date - first(Date[order(Date)])))

# Fit generalized additive model
all.f$Genotype <- factor(all.f$Genotype)

#Load Gamm model
load("Output/kmodel_6.RData")
load("Output/kmodelmean_6.RData")




# Get fitted values
newdata <- expand.grid(Days=seq(0, max(all.f$Days)),
                       Genotype=levels(all.f$Genotype))
newdatamean <- expand.grid(Days=seq(0, max(all.f$Days)))
newdata$fit <- predict(mod$gam, newdata) 
newdatamean$fit <-predict(modMean$gam, newdatamean) 

#For fitted model of all genos
set.seed(681) 
Rbeta <- mvrnorm(n = 10000, coef(mod$gam), vcov(mod$gam))
Xp <- predict(mod$gam, newdata = newdata, type = "lpmatrix")
sim <- Xp  %*% t(Rbeta)

#For mean not confinced this is necassary... for the derivs
set.seed(900)
RbetaMean <- mvrnorm(n = 10000, coef(modMean$gam), vcov(modMean$gam))
XpMean <- predict(modMean$gam, newdata = newdatamean, type = "lpmatrix")
simMean <- XpMean %*% t(RbetaMean)

# Extract 99% confidence intervals from simulation
newdata$lci <- apply(sim, 1, quantile, 0.005)
newdata$uci <- apply(sim, 1, quantile, 0.995)
newdatamean$lci <- apply(simMean, 1, quantile, 0.005)
newdatamean$uci <- apply(simMean, 1, quantile, 0.995)

#NormalizeData
for(geno in newdata$Genotype) {
  sdata <- NA
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

# Get rid of fitted values beyond day range for which we have data

maxDays <- aggregate(all.f$Days, by=list(all.f$Genotype), FUN=max)
maxDays$Group.1<-factor(maxDays$Group.1)
graphData <-newdata
for(i in 1:nrow(maxDays)) {
  graphData[which(graphData$Genotype==maxDays[i, "Group.1"] & newdata$Days>maxDays[i, "x"]), "Days"]<-NA
}
graphData <- na.omit(graphData)

#Get day zero fit and raw values
dayZero <- newdata[which(newdata$Days==0),]
rawDayZero <- all.f[which(all.f$Days==0),]
rawDZ <- aggregate(rawDayZero$Y, by=list(rawDayZero$Genotype), FUN=mean)
colnames(rawDZ) <- c("Genotype", "AverageRawY")
combineDayZero <- merge(dayZero, rawDZ, by= c("Genotype"))
combinerawDZ <- merge(rawDZ, compareGeno, by=c("Genotype"))
mod_rawdz <- lm(combinerawDZ$AverageRawY ~ combinerawDZ$Result)
anova(mod_rawdz)
lsmeans(mod_rawdz, specs= "Result")

#twelve Post Heat Stress Day 32
compareGeno <- data.frame(NA)
dayThT <- newdata[which(newdata$Days==32),]

dayThT <-dayThT[order(dayThT$Genotype),]


for(i in 1:nrow(dayThT)){
  compareGeno[i, "Genotype"]<-dayThT[i, "Genotype"]
  if(dayThT[i, "rlci"]>newdatamean[which(newdatamean$Days==32), "rfvfm"]) {
    compareGeno[i, "Result"]<-"above"
  }
  else if(dayThT[i, "ruci"]<newdatamean[which(newdatamean$Days==32), "rfvfm"]){
    compareGeno[i, "Result"]<-"below"
  }
  else{
    compareGeno[i, "Result"]<-"not"
  }
}

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
      #points(Y~Days, data=all.f[which(all.f$Genotype==names(dfsplit)[i]),], col=colors[i])
    })
  }
  legend("bottomleft", legend=levels(df$Genotype), pch=15, bty="n",
         col=alpha(colors, 0.4), pt.cex=1.5,
         y.intersp=0.3)
}

#Define normal plotting function
plotfitsnormal <- function(data, geno, title, key) {
  df <- droplevels(subset(data, Genotype %in% geno))
  getColor <- colorRampPalette(brewer.pal(length(geno), "Dark2"))
  colors <- getColor(length(geno))
  dfsplit <- split(df, f=df$Genotype)
  plot(NA, xlim=range(df$Days), ylim=c(0, 1), xaxs="i", yaxs="i", xlab = "Days", ylab = "Standardized Y", main=title)
  if(key==T) {
    par(new=T)
    plot(temp$Projected ~ temp$Days, type = "l", ylim=c(25, 35), xaxs="i", axes=F, lwd = 2, xlim=c(0,39), xlab="", ylab="")
    lines(temp$Thres~temp$Days, col = 'orange', lty = 2)
    lines(averageTemp$mean ~ averageTemp$Days, lwd = (.5), xlim=c(0, 39))
    lines(averageTemp$smooth ~ averageTemp$Days, col = 'orange')
    axis(side=4)
    mtext("Tempearture (oC)", side=4)
    par(new=T)
    plot(NA, xlim=range(df$Days), ylim=c(0, 1), xaxs="i", yaxs="i", xlab = "Days", ylab = "Standardized Y", main=title)
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
compareGeno <- compareGeno[, 2:3]
plotfits(graphData, c("1721", "1735"))
results_DayThT<-compareGeno[, "Genotype"]
plotfitsnormal(graphData, results_DayThT, NA, F)
lines(rfvfm~Days, data=newdatamean, lwd=2, col="red")
addpoly(newdatamean$Days, newdatamean$rlci, newdatamean$ruci, col=alpha("red", 0.3))


genoInfo <- read.csv("Data/CollectionR/genoinfowMMM.csv")

combineInfo <- merge(compareGeno, genoInfo, by = c("Genotype"))
combineInfo$Result <- factor(combineInfo$Result, levels=c("above", "not", "below"))
combineInfo$Nursery <- as.factor(combineInfo$Nursery)

#Take out struggle bus 1732
combineInfo[which(combineInfo$Genotype=="1732"),] <- NA

#Nursery
plot(Result~ Nursery, data=combineInfo)
plot(Nursery~ Result, data=combineInfo)
combineInfo$Years.in.nursery <- as.integer(combineInfo$Years.in.nursery)
boxplot(combineInfo$Years.in.nursery ~ combineInfo$Result)
mod_years <- lm(Years.in.nursery~Result, data=combineInfo)
anova(mod_years)
TukeyHSD(aov(mod_years))
cld(lsmeans(mod_years, "Result"))
#Nursery MMM
plot(nurseMMM~ Result, data=combineInfo)
mod_mmmnurse <- lm(nurseMMM~Result, data=combineInfo)
anova(mod_mmmnurse)
TukeyHSD(aov(mod_mmmnurse))
cld(lsmeans(mod_mmmnurse, "Result"))
# Depth
boxplot(Depth ~ Result, data = combineInfo)

mod_depth <- lm(Depth ~ Result, data = combineInfo)
anova(mod_depth)
TukeyHSD(aov(mod_depth))
cld(lsmeans(mod_depth, "Result"))

#Lat
#Look at non parametric tests for groups
boxplot(combineInfo$Lat ~ combineInfo$Result, ylab="Latitude Collected", xlab= "Result")
plot(combineInfo$Lat ~ combineInfo$Result, type="p", ylab="Latitude Collected", xlab= "Result")
mod_lat <- lm(Lat ~ Result, data = combineInfo)
anova(mod_lat)
TukeyHSD(aov(mod_lat))
cld(lsmeans(mod_lat, "Result"))
check <- lsmeans(mod_lat, "Result")
pairs(check) #Look at pairs help

kruskal.test(Lat ~ Result, data=combineInfo)
posthoc.kruskal.nemenyi.test(combineInfo$Lat, combineInfo$Result, dist = "Chisquare")
pairwise.wilcox.test(combineInfo$Lat, combineInfo$Result)

mod_latdep <- lm(Lat * Depth ~ Result, data = combineInfo)
anova(mod_latdep)

lat_summary <- ddply(combineInfo, c("Result"), summarise, N    = length(Lat),
                          mean = mean(Lat),
                          sd   = sd(Lat),
                          se   = sd / sqrt(N)
)
lat_summary$number <- seq.int(nrow(lat_summary))
plot(mean ~ number, data=lat_summary, xlim=c(1,3), xaxt="n", xlab="Result", ylab="Latitude", ylim=c(25.35, 25.7), pch=19)
arrows(lat_summary$number, lat_summary$mean-lat_summary$se, lat_summary$number, lat_summary$mean+lat_summary$se, length=0.05, angle=90, code=3)
axis(side=1,at=lat_summary$number,labels=lat_summary$Result)
shapiro.test(combineInfo$Lat)
bartlett.test(Lat~Result, data=combineInfo)
leveneTest(Lat~Result, data=combineInfo)
fligner.test(Lat~Result, data=combineInfo)

#MMM
plot(combineInfo$MMM ~ combineInfo$Result, ylab="MMM", xlab= "Result")
combineInfo$number <- seq.int(nrow(combineInfo))
plot(MMM~number, data=combineInfo[which(combineInfo$Result=="above"),], col="blue", pch=19, ylim=c(29.5, 30.5) )
points(MMM~number, data=combineInfo[which(combineInfo$Result=="not"),], col="green", pch=19)
points(MMM~number, data=combineInfo[which(combineInfo$Result=="below"),], col="red", pch=19)
mod_MMM <- lm(MMM ~ Result, data = combineInfo)
anova(mod_MMM)
TukeyHSD(aov(mod_MMM))
cld(lsmeans(mod_MMM, "Result"))
kruskal.test(MMM ~ Result, data=combineInfo)
mmm_summary <- ddply(combineInfo, c("Result"), summarise, N    = length(MMM),
                     mean = mean(MMM),
                     sd   = sd(MMM),
                     se   = sd / sqrt(N)
)
mmm_summary$number <- seq.int(nrow(mmm_summary))
plot(mean ~ number, data=mmm_summary, xlim=c(1,3), xaxt="n", xlab="Result", ylab="MMM", ylim=c(29.75, 30.25), pch=19)
arrows(mmm_summary$number, mmm_summary$mean-mmm_summary$se, mmm_summary$number, mmm_summary$mean+mmm_summary$se, length=0.05, angle=90, code=3)
axis(side=1,at=lat_summary$number,labels=lat_summary$Result)
#Long
boxplot(combineInfo$Long ~ combineInfo$Result, ylab="Latitude Collected", xlab= "Result")
plot(combineInfo$Long ~ combineInfo$Result, type="p", ylab="Latitude Collected", xlab= "Result")
mod_long <- lm(Long ~ Result, data = combineInfo)
anova(mod_long)
TukeyHSD(aov(mod_long))
cld(lsmeans(mod_long, "Result"))

#Test andrew's email hyptohesis
a <- combineInfo[which(combineInfo$Result=="above"),]
hist(a$Years.in.nursery)
hist(a$Lat)
plot(a$Nursery)
#z test 
z.test = function(x, mu, var){
  zeta = (mean(x) - mu) / (sqrt(var / length(x)))
  score = round((mean(x)-mu)/(var/sqrt(length(x))),3)
  return(zeta)
}
# Test 25.42 +- .03 as most southerly
z.test(a$Lat, 25.42, .03)
rnorm(11, 25.42, .03)
t.test(a$Lat, alternative =  "two.sided", mu = 25.42, var.equal = T, conf.level = .05)
# Test 25.55 +- .06 as central
t.test(a$Lat, alternative =  "two.sided", mu = 25.55, var.equal = T, conf.level = .05)
z.test(a$Lat, 25.55, .06)
#Test 25.63 +- .08
t.test(a$Lat, alternative =  "two.sided", mu = 25.63, var.equal = T, conf.level = .05)
z.test(a$Lat, 25.63, .08)

#below
b <- combineInfo[which(combineInfo$Result=="below"),]
hist(b$Years.in.nursery)
plot(b$Nursery)
# Test 25.42 +- .03 as most southerly
t.test(b$Lat, alternative =  "two.sided", mu = 25.42, var.equal = T, conf.level = .05)
z.test(b$Lat, 25.42, .03)
# Test 25.55 +- .06 as central
z.test(b$Lat, 25.55, .06)
t.test(b$Lat, alternative =  "two.sided", mu = 25.55, var.equal = T, conf.level = .05)
#Test 25.63 +- .08
t.test(b$Lat, alternative =  "two.sided", mu = 25.63, var.equal = T, conf.level = .05)
z.test(b$Lat, 25.63, .08)

#not
n <- combineInfo[which(combineInfo$Result=="not"),]
hist(n$Years.in.nursery)
# Test 25.42 +- .03 as most southerly
t.test(n$Lat, alternative =  "two.sided", mu = 25.42, var.equal = T, conf.level = .05)
z.test(b$Lat, 25.42, .03)
# Test 25.55 +- .06 as central
t.test(n$Lat, alternative =  "two.sided", mu = 25.55, var.equal = T, conf.level = .05)
z.test(n$Lat, 25.55, .06)
#Test 25.63 +- .08
t.test(n$Lat, alternative =  "two.sided", mu = 25.63, var.equal = T, conf.level = .05)
z.test(n$Lat, 25.63, .08)


#combine rfvfm to genoinfo to do linear regression on Lat day 32 values
combinePredicValueThT <- merge(dayThT, genoInfo, by = c("Genotype"))
plot(Lat ~ rfvfm, data=combinePredicValueThT)

fitLatThT<- lm(Lat ~ rfvfm, data=combinePredicValueThT)
abline(fitLatThT, col="red")
anova(fitLatThT)
library(effects)
plot(Effect("rfvfm", fitLatThT))
abline(fitLatThT, col="red")

#Do on MMM
plot(MMM ~ rfvfm, data=combinePredicValueThT, pch=19)

fitmmmThT<- lm(MMM ~ rfvfm, data=combinePredicValueThT)
abline(fitmmmThT, col="red")
anova(fitmmmThT)
library(effects)
plot(Effect("rfvfm", fitmmmThT))
abline(fitLatThT, col="red")

#combine rfvfm to genoinfo to do linear regression on Lat day 0 values
combinePredicValueZ <- merge(combineDayZero, genoInfo, by = c("Genotype"))
plot(AverageRawY ~ Lat, data=combinePredicValueZ)

fitLatZ<- lm(AverageRawY ~ Lat, data=combinePredicValueZ)
anova(fitLatZ)
abline(fitLatZ, col="red")


#Redo everything with struggle bus removed


#Initial FvFm 

intialsFvFm <- graphData[which(graphData$Days==0),]
merge.initals <- merge(compareGeno, intialsFvFm, by=c("Genotype"))
boxplot(fit ~ Result, data=merge.initals)

#Lat vs MMM
plot(MMM ~ Lat, data=combineInfo[which(combineInfo$Result=="above"),], col=alpha(colors[1], .5), pch=19, xlab="Latitude", ylab="MMM", ylim=c(29.8, 30.25), xlim=c(25.2, 25.8))
points(MMM ~ Lat, data=combineInfo[which(combineInfo$Result=="not"),], col=alpha(colors[2], .5), pch=19)
points(MMM ~ Lat, data=combineInfo[which(combineInfo$Result=="below"),], col=alpha(colors[3], .5), pch=19)
legend("topright", legend = c("above", "not", "below"), pch = 19, col=c(alpha(colors[1], .5),alpha(colors[2], .5), alpha(colors[3], .5)))
plot(MMM ~ Lat, data=combineInfo, type="p", xlab="Latitude Collected", ylab= "MMM")
mod_latmmm <- lm(MMM ~ Lat, data=combineInfo)
anova(mod_latmmm)
abline(mod_latmmm, col="red")

mod_latMMM <- lm(MMM ~ Result * Lat, data = combineInfo)
anova(mod_latMMM)
TukeyHSD(aov(mod_latMMM))
cld(lsmeans(mod_latMMM, "Result"))

##Square MMM
squareMMM<-read.csv("Data/CollectionR/squareMMM.csv")
plot(MMM~Lat, data = squareMMM, pch=19)
mod_latmmmsquare <- lm(MMM ~ Lat, data=squareMMM)
anova(mod_latmmmsquare)
abline(mod_latmmmsquare, col="red")
modMMMLat <- gamm4(MMM ~ s(Lat, k=5), random=~(1|MMM), data=squareMMM)
plot(modMMMLat$gam)
summary(modMMMLat$gam)
anova(modMMMLat$gam)
