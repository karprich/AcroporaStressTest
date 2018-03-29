#Variables and inclusion
library(lattice)
library(reshape2)
library(lsmeans)
library(gamm4)
library(dplyr)
library(MASS)
library(scales)
library(RColorBrewer)
library(mgcv)
library(zoo)
library(multcomp)
library(multcompView)
library(effects)
library(lme4)
library(lsmeans)

options(stringsAsFactors = FALSE)

#Import CSV for blast data, chl, and symbiont?
options(stringsAsFactors = FALSE)
blast <- read.csv("Data/SampleExtraction/rData/Blastedcoralsinfo.csv", header = T)
chlmeas <- read.csv("Data/SampleExtraction/rData/chlmeas.csv", header = T)
cellcount <- read.csv("Data/SampleExtraction/rData/Cellcountdata.csv", header=T)


#Date sampled as posix to days
blast$DateSampled <- as.POSIXct(blast$DateSampled, format = "%d/%m/%y")
blast$Date.Blasted <- as.POSIXct(blast$Date.Blasted, format = "%d/%m/%y")

#Calculate SA of coral fragment
blast$averageD <- (blast$WidthO+blast$WidthT)/2 #average D
#Widest D
for(i in 1:nrow(blast)) {
if(blast[i, "WidthO"]>=blast[i, "WidthT"]){
  blast[i, "maxD"]<-blast[i, "WidthO"]
}else {
  blast[i, "maxD"] <- blast[i, "WidthT"]
}
}
blast$SA <- blast$Length * blast$averageD * pi #* .44 # average diameter * pi * h *.44+- .05 conversion factor from Nauman 2009 paper take out factor
#convert SA from mm^2 to cm^2
blast$SA <- blast$SA/100

#make total volume 4.5 for all NAs
blast[which(is.na(blast$TotalVolume)), "TotalVolume"] <- 4.5 

#apply MeOH pheophytin correction of .01
chlmeas$Pheophytin.ug.mL. <- as.numeric(as.character(chlmeas$Pheophytin.ug.mL.)) - .01

#Merge spread sheets by genotype, frag ID and time point
combineExtraction <- merge(blast, chlmeas, by = c("TimePoint", "Genotype", "FragID")) #With Frag ID as factor. Misses 19 obs
combineExtraction <- merge(combineExtraction, compareGeno, by=c("Genotype"))
combineExtraction$Days <- as.integer(combineExtraction$DateSampled-as.POSIXct("2017-05-22"))

#Conservative Model and removal
#Remove time point 2
combineExtraction[which(combineExtraction$TimePoint==2), ] <- NA

#Remove time point 3
combineExtraction[which(combineExtraction$TimePoint==3), ] <- NA
#Remove time point 8
combineExtraction[which(combineExtraction$TimePoint==8), ] <- NA
#Remove time point 5 blasted on date October 11 2017 due to reverse filter
combineExtraction[which(combineExtraction$Date.Blasted==as.POSIXct("2017-10-11")), ] <- NA
#Remove time point 7 from 1530-1630 due to reverse filter
combineExtraction[which(combineExtraction$Date.Blasted==as.POSIXct("2017-10-30") & combineExtraction$Genotype=="1725"), ] <- NA
combineExtraction[which(combineExtraction$Date.Blasted==as.POSIXct("2017-10-30") & combineExtraction$Genotype=="1726"), ] <- NA
combineExtraction[which(combineExtraction$Date.Blasted==as.POSIXct("2017-10-30") & combineExtraction$Genotype=="1727"), ] <- NA
combineExtraction[which(combineExtraction$Date.Blasted==as.POSIXct("2017-10-30") & combineExtraction$Genotype=="1729"), ] <- NA
combineExtraction[which(combineExtraction$Date.Blasted==as.POSIXct("2017-10-30") & combineExtraction$Genotype=="1731"), ] <- NA
combineExtraction[which(combineExtraction$Date.Blasted==as.POSIXct("2017-10-30") & combineExtraction$Genotype=="1733"), ] <- NA
combineExtraction[which(combineExtraction$Date.Blasted==as.POSIXct("2017-10-30") & combineExtraction$Genotype=="1751"), ] <- NA
#Remove outlier 1757 time point 1
combineExtraction[which(combineExtraction$TimePoint==1 & combineExtraction$Genotype=="1757"), ] <- NA
#Remove Time point 5 time point 5 as a high value
combineExtraction[which(combineExtraction$TimePoint==5 & combineExtraction$Genotype=="1735"), ] <- NA


#Rearrange Table 
combineExtraction$Result <- as.factor(combineExtraction$Result)
combineExtraction$Genotype <- as.factor(combineExtraction$Genotype)
combineExtraction$Days <- as.integer(combineExtraction$Days)


#Determine ug in sample for chl and phe
combineExtraction$totalchl <- combineExtraction$Chlorophyl.ug.mL. * combineExtraction$mLMeOH * (1/combineExtraction$Dilution)
combineExtraction$totalchl <- combineExtraction$totalchl * combineExtraction$TotalVolume / combineExtraction$SampleFiltered
#determine chl and phe per cm^2
combineExtraction$totalchlcm <- combineExtraction$totalchl / combineExtraction$SA


#Plot points
plot(totalchlcm ~ Days, data=combineExtraction[which(combineExtraction$Result=="not"),])
points(totalchlcm ~ Days, data=combineExtraction[which(combineExtraction$Result=="above"),], col="blue")
points(totalchlcm ~ Days, data=combineExtraction[which(combineExtraction$Result=="below"),], col="red")

boxplot(totalchlcm ~ Days, data=combineExtraction[which(combineExtraction$Result=="not"),])
boxplot(totalchlcm ~ Days, data=combineExtraction[which(combineExtraction$Result=="above"),])
boxplot(totalchlcm ~ Days, data=combineExtraction[which(combineExtraction$Result=="below"),])
boxplot(totalchlcm~Days, data = combineExtraction)
#Plots

plot(totalchlcm ~Days, data = combineExtraction)
boxplot(totalchlcm ~ Days, data=combineExtraction)
xyplot(totalchlcm ~ Days | Genotype, data=combineExtraction, type=c("p"))
xyplot(totalchlcm ~ Days | Result, data=combineExtraction, type=c("p"))
xyplot(totalchlcm ~ TimePoint | Genotype, data=combineExtraction, type=c("p"))


#Combine with results

combineExtraction$TimePoint <- as.factor(combineExtraction$TimePoint)
combineExtraction$Days <- as.factor(combineExtraction$Days)

xyplot(totalchlcm ~ Days | Result, data=combineExtraction, type=c("p"))


xyplot(totalchlcm ~ Result | Days, data=combineExtraction, type=c("p"))
xyplot(totalchlcm ~ Result | TimePoint, data=combineExtraction, type=c("p"))

mod_chlbasic <- lm(totalchlcm ~ Result, data=combineExtraction)
anova(mod_chlbasic)
TukeyHSD(aov(mod_chlbasic))
cld(lsmeans(mod_chlbasic, "Result"))

mod_chldet <- lm(totalchlcm ~ factor(Days) * Result, data=combineExtraction)
anova(mod_chldet)
TukeyHSD(aov(mod_chldet))
cld(lsmeans(mod_chldet, specs="Days"))
cld(lsmeans(mod_chldet, "Result"))

View(cld(lsmeans(mod_chldet, specs = c("Days", "Result"))))

#Model by genotype
mod_geno <- lm(totalchlcm ~ Genotype, data=combineExtraction)
mod_genoDays <- lm(totalchlcm ~ Days * Genotype, data=combineExtraction)
mod_days <- lm(totalchlcm ~ Days, data= combineExtraction)
mod_chlresult <- lm(totalchlcm ~ Days * Result, data= combineExtraction)
anova(mod_geno)
TukeyHSD(aov(mod_geno))
cld(lsmeans(mod_geno, "Genotype"))




anova(mod_days)
TukeyHSD(aov(mod_days))
cld(lsmeans(mod_days, "Days"))

anova(mod_genoDays)
TukeyHSD(aov(mod_genoDays))
cld(lsmeans(mod_genoDays, spec= c("Days", "Genotype")))

combineExtraction$Days <- as.integer(combineExtraction$Days)
plot(NA, xlim=c(0,40), ylim = c(0, 10), xlab="Days", ylab="Chl per cm^2", xaxs="i", yaxs="i")
points(totalchlcm ~ Days, data = combineExtraction[which(combineExtraction$Result=="above"),], type="p", col="blue", xlim=c(0,45))
points(totalchlcm ~ Days, data = combineExtraction[which(combineExtraction$Result=="below"),], type="p", col="red")
points(totalchlcm ~ Days, data = combineExtraction[which(combineExtraction$Result=="not"),], type="p", col="green")


#Mean SE plots

above_chl <- combineExtraction[which(combineExtraction$Result=="above"),]
above_chl$Relative <- above_chl$totalchlcm/chl_summary[which(chl_summary$Result=="above" & chl_summary$Days==2), "mean"]
above_chl$Days <- as.integer(above_chl$DateSampled-as.POSIXct("2017-05-22"))
abovechl_summary <- ddply(above_chl, c("Days"), summarise, N    = length(totalchlcm),
                           mean = mean(totalchlcm),
                           sd   = sd(totalchlcm),
                           se   = sd / sqrt(N),
                           rmean = mean(Relative),
                           rsd= sd(Relative),
                           rse=rsd/sqrt(N)
)
abovechl_summary$Result <- "above"
mod_aboveChlsum <- lm(mean ~ Days, data=abovechl_summary)
plot(mean ~ Days, abovechl_summary, xaxs="i", yaxs="i", ylim=c(0, 5), xlim=c(0,40))
abline(mod_aboveChlsum)
arrows(abovechl_summary$Days, abovechl_summary$mean-abovechl_summary$se, abovechl_summary$Days, abovechl_summary$mean+abovechl_summary$se, length=0.05, angle=90, code=3)

modrel_aboveChlsum <- lm(rmean ~ Days, data=abovechl_summary)
plot(rmean ~ Days, abovechl_summary)
abline(modrel_aboveChlsum)

#Below
below_chl <- combineExtraction[which(combineExtraction$Result=="below"),]
below_chl$Relative <- below_chl$totalchlcm/chl_summary[which(chl_summary$Result=="below" & chl_summary$Days==2), "mean"]
below_chl$Days <- as.integer(below_chl$DateSampled-as.POSIXct("2017-05-22"))
belowchl_summary <- ddply(below_chl, c("Days"), summarise, N    = length(totalchlcm),
                          mean = mean(totalchlcm),
                          sd   = sd(totalchlcm),
                          se   = sd / sqrt(N),
                          rmean = mean(Relative),
                          rsd= sd(Relative),
                          rse=rsd/sqrt(N)
)
belowchl_summary$Result <- "below"
mod_belowChlsum <- lm(mean ~ Days, data=belowchl_summary)
plot(mean ~ Days, belowchl_summary)
abline(mod_belowChlsum)
arrows(belowchl_summary$Days, belowchl_summary$mean-belowchl_summary$se, belowchl_summary$Days, belowchl_summary$mean+belowchl_summary$se, length=0.05, angle=90, code=3)

modrel_belowChlsum <- lm(rmean ~ Days, data=belowchl_summary)
plot(rmean ~ Days, belowchl_summary)
abline(modrel_belowChlsum)

#Not
not_chl <- combineExtraction[which(combineExtraction$Result=="not"),]
not_chl$Relative <- not_chl$totalchlcm/chl_summary[which(chl_summary$Result=="not" & chl_summary$Days==2), "mean"]
not_chl$Days <- as.integer(not_chl$DateSampled-as.POSIXct("2017-05-22"))
notchl_summary <- ddply(not_chl, c("Days"), summarise, N    = length(totalchlcm),
                          mean = mean(totalchlcm),
                          sd   = sd(totalchlcm),
                          se   = sd / sqrt(N),
                        rmean = mean(Relative),
                        rsd= sd(Relative),
                        rse=rsd/sqrt(N)
)
notchl_summary$Result <- "not"
mod_notChlsum <- lm(totalchlcm ~ Days, data=not_chl)
plot(mean ~ Days, notchl_summary)
abline(mod_notChlsum)
arrows(notchl_summary$Days, notchl_summary$mean-notchl_summary$se, notchl_summary$Days, notchl_summary$mean+notchl_summary$se, length=0.05, angle=90, code=3)

modrel_notChlsum <- lm(rmean ~ Days, data=notchl_summary)
plot(rmean ~ Days, notchl_summary)
abline(modrel_notChlsum)
arrows(notchl_summary$Days, notchl_summary$mean-notchl_summary$se, notchl_summary$Days, notchl_summary$mean+notchl_summary$se, length=0.05, angle=90, code=3)


# Test if slopes are differnt
combineExtraction$Days <- as.integer(combineExtraction$DateSampled-as.POSIXct("2017-05-22"))
chl_summary <- ddply(combineExtraction, c("Days", "Result"), summarise, N    = length(totalchlcm),
                          mean = mean(totalchlcm),
                          sd   = sd(totalchlcm),
                          se   = sd / sqrt(N)
)
mod_chl_summary <- lm(totalchlcm ~ Days *Result, data=combineExtraction)
anova(mod_chl_summary)


chl_lst <- lstrends(mod_chl_summary, "Result", var="Days")
pairs(chl_lst)
cld(lstrends(mod_chl_summary, "Result", var="Days"))
TukeyHSD(aov(mod_chl_summary))
cld(lsmeans(mod_chl_summary, "Result"))

#check relative
relcombinechl<- rbind(abovechl_summary, notchl_summary, belowchl_summary)
mod_chl_rel <- lm(rmean ~ Days *Result, data=relcombinechl)
anova(mod_chl_rel)
chl_rellst <- lstrends(mod_chl_rel, "Result", var="Days")
pairs(chl_rellst)
cld(lstrends(mod_chl_rel, "Result", var="Days"))

#plot on log scale
par(mfrow = c(1, 3))
points(mean ~ Days, data=abovechl_summary, ylab="Chl A (ug per cm^2)", xlab="Days", main="Above average genotypes", yaxs="i", col=alpha(colors[1], .9), pch=19, log="y")
abline(mod_aboveChlsum)
arrows(abovechl_summary$Days, abovechl_summary$mean-abovechl_summary$se, abovechl_summary$Days, abovechl_summary$mean+abovechl_summary$se, length=0.05, angle=90, code=3)
text(x = 32, y = 3, labels = expression("Adj. R"^2*"=0.7095"))
points(mean ~ Days, notchl_summary, ylim=c(0, 4), ylab="Chl A (ug per cm^2)", xlab="Days", main="Not different genotypes", yaxs="i", col=alpha(colors[2], .9), pch=19)
abline(mod_notChlsum)
arrows(notchl_summary$Days, notchl_summary$mean-notchl_summary$se, notchl_summary$Days, notchl_summary$mean+notchl_summary$se, length=0.05, angle=90, code=3)
text(x = 32, y = 3, labels = expression("Adj. R"^2*"=0.4903"))

plot(mean ~ Days, belowchl_summary, ylab="Chl A (ug per cm^2)", xlab="Days", main="Below average genotypes", yaxs="i", col=alpha(colors[3], .9), pch=19, log="y")
abline(mod_belowChlsum)
arrows(belowchl_summary$Days, belowchl_summary$mean-belowchl_summary$se, belowchl_summary$Days, belowchl_summary$mean+belowchl_summary$se, length=0.05, angle=90, code=3)
text(x = 32, y = 3, labels = expression("Adj. R"^2*"=0.8462"))
par(mfrow = c(1, 1))

#Relative chl plot

plot(rmean ~ Days, data=abovechl_summary, ylim=c(0, 1.5), ylab="Relative Chl A (ug per cm^2)", xlab="Days", main="Relative Chl A", yaxs="i", col=alpha(colors[1], .5), pch=19)
arrows(abovechl_summary$Days, abovechl_summary$rmean-abovechl_summary$rse, abovechl_summary$Days, abovechl_summary$rmean+abovechl_summary$rse, length=0.05, angle=90, code=3)
points(rmean ~ Days, notchl_summary, ylim=c(0, 4), ylab="Chl A (ug per cm^2)", xlab="Days", main="Not different genotypes", yaxs="i", col=alpha(colors[2], .5), pch=19)
arrows(notchl_summary$Days, notchl_summary$rmean-notchl_summary$rse, notchl_summary$Days, notchl_summary$rmean+notchl_summary$rse, length=0.05, angle=90, code=3)

points(rmean ~ Days, belowchl_summary, ylim=c(0, 4), ylab="Chl A (ug per cm^2)", xlab="Days", main="Below average genotypes", yaxs="i", col=alpha(colors[3], .5), pch=19)
arrows(belowchl_summary$Days, belowchl_summary$rmean-belowchl_summary$rse, belowchl_summary$Days, belowchl_summary$rmean+belowchl_summary$rse, length=0.05, angle=90, code=3)

legend("topright", legend = c("above", "not", "below"), pch = 19, col=c(alpha(colors[1], .5),alpha(colors[2], .5), alpha(colors[3], .5)))

chlrel <- rbind(notchl_summary, abovechl_summary, belowchl_summary)
xyplot(rmean ~ Days, data=chlrel, groups = Result)

#Gamm model chl
addpoly <- function(x, y1, y2, col=alpha("lightgrey", 0.8), ...){
  ii <- order(x)
  y1 <- y1[ii]
  y2 <- y2[ii]
  x <- x[ii]
  polygon(c(x, rev(x)), c(y1, rev(y2)), col=col, border=NA, ...)
}

#Omit NA values
combineExtraction$Days <- as.integer(combineExtraction$DateSampled-as.POSIXct("2017-05-22"))
combineExtraction <- combineExtraction[which(!is.na(combineExtraction$totalchlcm)),]
combineExtraction[which(combineExtraction$totalchlcm==0), "totalchlcm"] <- 0.001 # ChANGE ALL 0s to .001
combineExtraction$logchl <- log10(combineExtraction$totalchlcm)

mod_gamChl <- gamm4(logchl ~ Genotype + s(Days, k=4, by=Genotype), random=~(1|FragID), data=combineExtraction)
modMean_gamChl <- gamm4(logchl ~ s(Days, k=4), random=~(1|Genotype), data=combineExtraction)
modResult_gam <- gamm4(logchl ~ Result + s(Days, k=4, by=Result), random=~(1|Genotype)+(1|FragID), data=combineExtraction)

newdata_chl <- expand.grid(Days=seq(0, max(combineExtraction$Days)),
                       Genotype=levels(combineExtraction$Genotype))
newdatamean_chl <- expand.grid(Days=seq(0, max(combineExtraction$Days)))
data_chlResult <- expand.grid(Days=seq(0, max(combineExtraction$Days)), Result=levels(combineExtraction$Result))
newdata_chl$fit <- predict(mod_gamChl$gam, newdata_chl) 
newdatamean_chl$fit <-predict(modMean_gamChl$gam, newdatamean_chl)
data_chlResult$fit <-predict(modResult_gam$gam, data_chlResult) 
#For fitted model of all genos
set.seed(231) 
Rbeta <- mvrnorm(n = 10000, coef(mod_gamChl$gam), vcov(mod_gamChl$gam))
Xp <- predict(mod_gamChl$gam, newdata = newdata_chl, type = "lpmatrix")
sim_chl <- Xp %*% t(Rbeta)


#For mean not confinced this is necassary... for the derivs
set.seed(543)
RbetaMean <- mvrnorm(n = 10000, coef(modMean_gamChl$gam), vcov(modMean_gamChl$gam))
XpMean <- predict(modMean_gamChl$gam, newdata = newdatamean_chl, type = "lpmatrix")
simMean_chl <- XpMean %*% t(RbetaMean)

#For mean not confinced this is necassary... for the derivs
set.seed(729)
RbetaMean <- mvrnorm(n = 10000, coef(modResult_gam$gam), vcov(modResult_gam$gam))
XpMean <- predict(modResult_gam$gam, newdata = data_chlResult, type = "lpmatrix")
simchl_Result <- XpMean %*% t(RbetaMean)


# Extract 90% confidence intervals from simulation
newdata_chl$lci <- apply(sim_chl, 1, quantile, 0.05)
newdata_chl$uci <- apply(sim_chl, 1, quantile, 0.95)
newdatamean_chl$lci <- apply(simMean_chl, 1, quantile, 0.05)
newdatamean_chl$uci <- apply(simMean_chl, 1, quantile, 0.95)
data_chlResult$lci <- apply(simchl_Result, 1, quantile, 0.05)
data_chlResult$uci <- apply(simchl_Result, 1, quantile, 0.95)

#Plot fits chl function
plotfitschlres <- function(data, res, title, key) {
  df <- droplevels(subset(data, Result %in% res))
  getColor <- colorRampPalette(brewer.pal(length(res), "Dark2"))
  colors <- getColor(length(res))
  dfsplit <- split(df, f=df$Result)
  plot(NA, xlim=c(0, 40), ylim=range(-3, 1), xaxs="i", yaxs="i", xlab = "Days", ylab = "log(chlorophyll A) (ug per cm^2)", main=title)
  if(key==T) {
    par(new=T)
    plot(temp$Projected ~ temp$Days, type = "l", ylim=c(25, 35), xaxs="i", axes=F, lwd = 2, xlim=c(0,39), xlab="", ylab="")
    lines(temp$Thres~temp$Days, col = 'orange', lty = 2)
    lines(averageTemp$mean ~ averageTemp$Days, lwd = (.5), xlim=c(0, 39))
    lines(averageTemp$smooth ~ averageTemp$Days, col = 'red')
    axis(side=4)
    mtext("Tempearture (oC)", side=4)
    par(new=T)
    plot(NA, xlim=c(0, 40), ylim=c(0, 10), xaxs="i", yaxs="i", xlab = "Days", ylab = "Chl mg per cm^2", main=title)
  }
  for (i in 1:length(dfsplit)) {
    with(dfsplit[[i]], {
      lines(Days, fit)
      addpoly(Days, lci, uci, col=alpha(colors[i], 0.3))
    })
  }
  
  legend("bottom", legend=levels(df$Result), pch=15, bty="n",
         col=alpha(colors, 0.4), pt.cex=1.5,
         y.intersp=0.3)
}

#Normalize plotting
plotfitschlresnormal <- function(data, res, title, key) {
  df <- droplevels(subset(data, Result %in% res))
  getColor <- colorRampPalette(brewer.pal(length(res), "Dark2"))
  colors <- getColor(length(res))
  dfsplit <- split(df, f=df$Result)
  plot(NA, xlim=c(0, 40), ylim=c(-3,1), xaxs="i", yaxs="i", xlab = "Days", ylab = "Relative Chl ug per cm^2", main=title)
  if(key==T) {
    par(new=T)
    plot(temp$Projected ~ temp$Days, type = "l", ylim=c(25, 35), xaxs="i", axes=F, lwd = 2, xlim=c(0,39), xlab="", ylab="")
    lines(temp$Thres~temp$Days, col = 'orange', lty = 2)
    lines(averageTemp$mean ~ averageTemp$Days, lwd = (.5), xlim=c(0, 39))
    lines(averageTemp$smooth ~ averageTemp$Days, col = 'red')
    axis(side=4)
    mtext("Tempearture (oC)", side=4)
    par(new=T)
    plot(NA, xlim=c(0, 40), ylim=c(0, 10), xaxs="i", yaxs="i", xlab = "Days", ylab = "relative Chl ug per cm^2", main=title)
  }
  for (i in 1:length(dfsplit)) {
    with(dfsplit[[i]], {
      lines(Days, rchl)
      addpoly(Days, rlci, ruci, col=alpha(colors[i], 0.3))
    })
  }
  
  legend("bottom", legend=levels(df$Result), pch=15, bty="n",
         col=alpha(colors, 0.4), pt.cex=1.5,
         y.intersp=0.3)
}

#Plot fits chl function
plotfitschlgeno <- function(data, geno, title, key) {
  df <- droplevels(subset(data, Genotype %in% geno))
  getColor <- colorRampPalette(brewer.pal(length(geno), "Dark2"))
  colors <- getColor(length(geno))
  dfsplit <- split(df, f=df$Genotype)
  plot(NA, xlim=c(0, 40), ylim=range(df$fit), xaxs="i", xlab = "Days", ylab = "Chl mg per cm^2", main=title)
  if(key==T) {
    par(new=T)
    plot(temp$Projected ~ temp$Days, type = "l", ylim=c(25, 35), xaxs="i", axes=F, lwd = 2, xlim=c(0,39), xlab="", ylab="")
    lines(temp$Thres~temp$Days, col = 'orange', lty = 2)
    lines(averageTemp$mean ~ averageTemp$Days, lwd = (.5), xlim=c(0, 39))
    lines(averageTemp$smooth ~ averageTemp$Days, col = 'red')
    axis(side=4)
    mtext("Tempearture (oC)", side=4)
    par(new=T)
    plot(NA, xlim=c(0, 40), ylim=range(df$fit), xaxs="i", xlab = "Days", ylab = "Chl mg per cm^2", main=title)
  }
  for (i in 1:length(dfsplit)) {
    with(dfsplit[[i]], {
      lines(Days, fit)
      addpoly(Days, lci, uci, col=alpha(colors[i], 0.3))
    })
  }
  
  legend("bottom", legend=levels(df$Genotype), pch=15, bty="n",
         col=alpha(colors, 0.4), pt.cex=1.5,
         y.intersp=0.3)
}

#Normalize chl value
for(res in data_chlResult$Result) {
  sdata <- NA
  sdata <- data_chlResult[which(data_chlResult$Result==res),]
  for(day in sdata$Days) {
    sdata[which(sdata$Days ==day), "rchl"] <- sdata[which(sdata$Days==day), "fit"] / sdata[which(sdata$Days==0), "fit"]
    sdata[which(sdata$Days ==day), "rlci"] <- sdata[which(sdata$Days==day), "lci"] / sdata[which(sdata$Days==0), "fit"]
    sdata[which(sdata$Days ==day), "ruci"] <- sdata[which(sdata$Days==day), "uci"] / sdata[which(sdata$Days==0), "fit"]
  }
  data_chlResult[which(data_chlResult$Result==res),"rchl"] <- sdata[, "rchl"]
  data_chlResult[which(data_chlResult$Result==res), "rlci"] <- sdata[, "rlci"]
  data_chlResult[which(data_chlResult$Result==res), "ruci"] <- sdata[, "ruci"]
}

#Normalize Mean
for(day in newdatamean_chl$Days) {
  newdatamean_chl[which(newdatamean_chl$Days ==day), "rchl"] <- newdatamean_chl[which(newdatamean_chl$Days==day), "fit"] / newdatamean_chl[which(newdatamean_chl$Days==0), "fit"]
  newdatamean_chl[which(newdatamean_chl$Days ==day), "rlci"] <- newdatamean_chl[which(newdatamean_chl$Days==day), "lci"] / newdatamean_chl[which(newdatamean_chl$Days==0), "fit"]
  newdatamean_chl[which(newdatamean_chl$Days ==day), "ruci"] <- newdatamean_chl[which(newdatamean_chl$Days==day), "uci"] / newdatamean_chl[which(newdatamean_chl$Days==0), "fit"]
}
#Plot
plotfitschlgeno(newdata_chl, c("1738", "1721", "1759", "1750", "1741"), NA, F)
plotfitschlgeno(newdata_chl, c("1721", "1735", "1748"), NA, F)
plotfitschlgeno()
#Plotting
xyplot(fit ~ Days | Genotype, data=newdata_chl, type="l")

xyplot(totalchlcm ~ Days | Genotype, data = combineExtraction)

plot(fit~Days, data = newdatamean_chl, type="l", xaxs="i", yaxs="i", ylim=c(0,10))
plot(fit~ Days, data=data_chlResult[which(data_chlResult$Result=="above"),], col="blue", type="l", ylim=c(0,10))

lines(fit~ Days, data=data_chlResult[which(data_chlResult$Result=="below"),], col="red")
lines(fit~ Days, data=data_chlResult[which(data_chlResult$Result=="not"),], col= "green")
legend("topright", legend = c("Mean", "Above", "Below", "Not"), lwd  = 2, col = c("black", "blue", "red", "green"))
plot(totalchlcm~ Days, data=combineExtraction[which(combineExtraction$Result=="not"),], col= "green")
#Rchl plots
plot(rchl~Days, data = newdatamean_chl, type="l", xaxs="i", yaxs="i", ylim=c(0,2))
lines(rchl~ Days, data=data_chlResult[which(data_chlResult$Result=="above"),], col="blue")

lines(rchl~ Days, data=data_chlResult[which(data_chlResult$Result=="below"),], col="red")
lines(rchl~ Days, data=data_chlResult[which(data_chlResult$Result=="not"),], col= "green")

plotfitschlresnormal(data_chlResult, c("above", "below", "not"), "Relative Chlorophyll A for the groups", F)

xyplot(fit ~ Days, data=data_chlResult, groups=Result, type="l")
xyplot(rchl ~ Days, data=data_chlResult, groups=Result, type="l")
# Result with mean in legend
plotfitschlres(data_chlResult, c("above", "below", "not"), "Categories of chlorophyll throughout heat stress", F)
lines(fit~Days, data = newdatamean_chl, ylim=c(0, 4), type="l", xaxs="i", yaxs="i", col="red")
legend(x = c(15.2,20), y=c(1.35, 1.45), legend = c("Mean"), pch=15, bty="n",
       col="red", pt.cex=1.5,
       y.intersp=0.3)
addpoly(newdatamean_chl$Days, newdatamean_chl$lci, newdatamean_chl$uci, alpha("red", .4))
lines(fit~Days, newdata[which(newdata$Genotype=="1735"),])

plotfitschlres(data_chlResult, c("above", "below", "not"), "Categories of chlorophyll A on log scale", F)

plotfits(newdata, c("1721, 1735", "1727"))
#Run analysis on Day 33 differences Day 32 iPam



chlDayThT <- combineExtraction[which(combineExtraction$Days==33),]
dayThT_chllm <- lm(totalchlcm ~ Result, data = chlDayThT)
dayThT_chllmGeno <- lm(totalchlcm ~ Genotype, data = chlDayThT)
anova(dayThT_chllm)
TukeyHSD(aov(dayThT_chllm))
cld(lsmeans(dayThT_chllm, "Result"))

anova(dayThT_chllmGeno)
TukeyHSD(aov(dayThT_chllmGeno))
cld(lsmeans(dayThT_chllmGeno, "Genotype"))

par(mfrow = c(1, 3))
plot(totalchlcm ~ Days, data=combineExtraction[which(combineExtraction$Result=="above"),], col="green", yaxs="i", xlab="Days", ylab="Ug Chl A per cm^2", main="Above average genotypes", ylim=c(0, 14))
boxplot(totalchlcm ~ Days, data=combineExtraction[which(combineExtraction$Result=="not"),], col="yellow", yaxs="i", xlab="Days", ylab="Ug Chl A per cm^2", main="Average genotypes", ylim=c(0, 14))
boxplot(totalchlcm ~ Days, data=combineExtraction[which(combineExtraction$Result=="below"),], col="red", yaxs="i", xlab="Days", ylab="Ug Chl A per cm^2", main="Below average genotypes", ylim=c(0, 14))
par(mfrow = c(1, 1))

#Cell count processing

cellcount$Genotype <- as.factor(cellcount$Genotype)
cellcount$totalCellspsq <- cellcount$Average.with.Dilution / .0001 # each box counted is .0001 mL
cellcount$totalcell <- cellcount$totalCellspsq * (cellcount$Volume.Lugols +.0375) 

combineCells <- merge(blast, cellcount, by=c("TimePoint", "Genotype"))
combineCells <- merge(combineCells, compareGeno, by=c("Genotype"))
combineCells$Result <- factor(combineCells$Result, levels = c("above", "not", "below"))
combineCells$Genotype <- as.factor(combineCells$Genotype)
combineCells$Days <- as.integer(combineCells$DateSampled-as.POSIXct("2017-05-22"))
combineCells$totalcellvol <- combineCells$totalcell * combineCells$TotalVolume / combineCells$Volume.Lugols
combineCells$totalcellcm <- combineCells$totalcellvol / combineCells$SA
plot(totalcellcm ~ Days, data=combineCells)
boxplot(totalcellcm ~ TimePoint, data=combineCells)
boxplot(totalcellcm ~ Result, data=combineCells[which(combineCells$TimePoint==7),])
xyplot(totalcellcm ~ Days | Result, data=combineCells)



#Symbionts per cm sq
par(mfrow = c(1, 3))
boxplot(totalcellcm ~ Days, data=combineCells[which(combineCells$Result=="not"),], col="green")
boxplot(totalcellcm ~ Days, data=combineCells[which(combineCells$Result=="above"),], col="yellow")
boxplot(totalcellcm ~ Days, data=combineCells[which(combineCells$Result=="below"),], col="red")
par(mfrow = c(1, 1))

#Stats
mod_symResult <- lm(totalcellcm ~ Result, data=combineCells )
anova(mod_symResult)
TukeyHSD(aov(mod_symResult))

mod_symResulttime <- lm(totalcellcm ~factor(Days) * Result, data=combineCells )
anova(mod_symResulttime)
TukeyHSD(aov(mod_symResulttime))
View(cld(lsmeans(mod_symResulttime, c("Result", "Days"))))


#Mean SE plots for sym

above_sym <- combineCells[which(combineCells$Result=="above"),]
above_sym$Relative <- above_sym$totalcellcm/sym_summary[which(sym_summary$Days==2 & sym_summary$Result=="above"), "mean"] 
above_sym$Days <- as.integer(above_sym$DateSampled-as.POSIXct("2017-05-22"))
abovesym_summary <- ddply(above_sym, c("Days"), summarise, N    = length(totalcellcm),
                          mean = mean(totalcellcm),
                          sd   = sd(totalcellcm),
                          se   = sd / sqrt(N),
                          rmean = mean(Relative),
                          rsd= sd(Relative),
                          rse = rsd/ sqrt(N)
)
mod_abovesymsum <- lm(mean ~ Days, data=abovesym_summary)
plot(mean ~ Days, abovesym_summary)
abline(mod_abovesymsum)
arrows(abovesym_summary$Days, abovesym_summary$mean-abovesym_summary$se, abovesym_summary$Days, abovesym_summary$mean+abovesym_summary$se, length=0.05, angle=90, code=3)

#Below
below_sym <- combineCells[which(combineCells$Result=="below"),]
below_sym$Relative <- below_sym$totalcellcm/sym_summary[which(sym_summary$Days==2 & sym_summary$Result=="below"), "mean"] 
below_sym$Days <- as.integer(below_sym$DateSampled-as.POSIXct("2017-05-22"))
belowsym_summary <- ddply(below_sym, c("Days"), summarise, N    = length(totalcellcm),
                          mean = mean(totalcellcm),
                          sd   = sd(totalcellcm),
                          se   = sd / sqrt(N),
                          rmean = mean(Relative),
                          rsd= sd(Relative),
                          rse = rsd/ sqrt(N)
)
mod_belowsymsum <- lm(mean ~ Days, data=belowsym_summary)
plot(mean ~ Days, belowsym_summary)
abline(mod_belowsymsum)
arrows(belowsym_summary$Days, belowsym_summary$mean-belowsym_summary$se, belowsym_summary$Days, belowsym_summary$mean+belowsym_summary$se, length=0.05, angle=90, code=3)
#Not
not_sym <- combineCells[which(combineCells$Result=="not"),]
not_sym$Relative <- not_sym$totalcellcm/sym_summary[which(sym_summary$Days==2 & sym_summary$Result=="not"), "mean"] 
not_sym$Days <- as.integer(not_sym$DateSampled-as.POSIXct("2017-05-22"))
notsym_summary <- ddply(not_sym, c("Days"), summarise, N    = length(totalcellcm),
                          mean = mean(totalcellcm),
                          sd   = sd(totalcellcm),
                          se   = sd / sqrt(N),
                          rmean = mean(Relative),
                          rsd= sd(Relative),
                          rse = rsd/ sqrt(N)

)
mod_notsymsum <- lm(mean ~ Days, data=notsym_summary)
plot(mean ~ Days, notsym_summary)
abline(mod_notsymsum)
arrows(notsym_summary$Days, notsym_summary$mean-notsym_summary$se, notsym_summary$Days, notsym_summary$mean+notsym_summary$se, length=0.05, angle=90, code=3)

#Plot
#par(mfrow = c(1, 3))
plot(mean ~ Days, data=abovesym_summary, ylim=c(0, 6E6), ylab="Symbionts per cm^2", xlab="Days", main="Above average genotypes", yaxs="i", col=alpha(colors[1], .9), pch=19)
#abline(mod_abovesymsum)
arrows(abovesym_summary$Days, abovesym_summary$mean-abovesym_summary$se, abovesym_summary$Days, abovesym_summary$mean+abovesym_summary$se, length=0.05, angle=90, code=3)
#text(x = 32, y = 5E6, labels = expression("Adj. R"^2*"=0.5393"))
points(mean ~ Days, notsym_summary, ylim=c(0, 6E6), ylab="Symbionts per cm^2", xlab="Days", main="Not different genotypes", yaxs="i", col=alpha(colors[2], .9), pch=19)
#abline(mod_notsymsum)
arrows(notsym_summary$Days, notsym_summary$mean-notsym_summary$se, notsym_summary$Days, notsym_summary$mean+notsym_summary$se, length=0.05, angle=90, code=3)
#text(x = 32, y = 5E6, labels = expression("Adj. R"^2*"=0.7462"))

points(mean ~ Days, belowsym_summary, ylim=c(0, 6E6), ylab="Symbionts per cm^2", xlab="Days", main="Below average genotypes", yaxs="i", col=alpha(colors[3], .9), pch=19)
#abline(mod_belowsymsum)
arrows(belowsym_summary$Days, belowsym_summary$mean-belowsym_summary$se, belowsym_summary$Days, belowsym_summary$mean+belowsym_summary$se, length=0.05, angle=90, code=3)
#text(x = 32, y = 5E6, labels = expression("Adj. R"^2*"=0.6943"))
#par(mfrow = c(1, 1))

#Differences at t2 for syms
symDayTon <- combineCells[which(combineCells$Days>17),]
symDayTon <- symDayTon[which(symDayTon$Days!=34),]
mod_symDayTon <- lm(totalcellcm ~ Days *Result, data=symDayTon)
anova(mod_symDayTon)
sym_lst <- lstrends(mod_symDayTon, "Result", var="Days")
pairs(sym_lst)
cld(lstrends(mod_symDayTon, "Result", var="Days"))

xyplot(totalcellcm ~ Days, data=symDayTon, groups = Result)

#Relative sym plots
plot(rmean ~ Days, data=abovesym_summary, ylim=c(0, 2.5), ylab="Relative Symbionts per cm^2", xlab="Days", main="Relative Symbionts for each group", yaxs="i", col=alpha(colors[1], .5), pch=19)
arrows(abovesym_summary$Days, abovesym_summary$rmean-abovesym_summary$rse, abovesym_summary$Days, abovesym_summary$rmean+abovesym_summary$rse, length=0.05, angle=90, code=3)
points(rmean ~ Days, notsym_summary, ylim=c(0, 6E6), ylab="Symbionts per cm^2", xlab="Days", main="Not different genotypes", yaxs="i", col=alpha(colors[2], .5), pch=19)
#abline(mod_notsymsum)
arrows(notsym_summary$Days, notsym_summary$rmean-notsym_summary$rse, notsym_summary$Days, notsym_summary$rmean+notsym_summary$rse, length=0.05, angle=90, code=3)
#text(x = 32, y = 5E6, labels = expression("Adj. R"^2*"=0.7462"))

points(rmean ~ Days, belowsym_summary, ylim=c(0, 6E6), ylab="Symbionts per cm^2", xlab="Days", main="Below average genotypes", yaxs="i", col=alpha(colors[3], .5), pch=19)
arrows(belowsym_summary$Days, belowsym_summary$rmean-belowsym_summary$rse, belowsym_summary$Days, belowsym_summary$rmean+belowsym_summary$rse, length=0.05, angle=90, code=3)
legend("topright", legend = c("above", "not", "below"), pch = 19, col=c(alpha(colors[1], .5),alpha(colors[2], .5), alpha(colors[3], .5)))
# Test if slopes are differnt
sym_summary <- ddply(combineCells, c("Days", "Result"), summarise, N    = length(totalcellcm),
                     mean = mean(totalcellcm),
                     sd   = sd(totalcellcm),
                     se   = sd / sqrt(N)
)
xyplot(mean ~ Days | Result, sym_summary, type=c("p", "r"))
mod_sym_summary <- lm(totalcellcm ~ Days *Result, data=combineCells)
anova(mod_sym_summary)
sym_lst <- lstrends(mod_sym_summary, "Result", var="Days")
pairs(sym_lst)
cld(lstrends(mod_sym_summary, "Result", var="Days"))
TukeyHSD(aov(mod_sym_summary))
cld(lsmeans(mod_sym_summary, "Result"))

combineAll <- merge(combineExtraction, combineCells, by=c("TimePoint", "Genotype"))
combineAll$chlpsym <-  combineAll$totalchlcm/combineAll$totalcellcm
boxplot(chlpsym ~ Days.x, data=combineAll)
boxplot(chlpsym ~ TimePoint, data=combineAll)
boxplot(totalcellcm ~ Result.x, data=combineAll[which(combineAll$TimePoint==6),])

xyplot(chlpsym ~ Days.x | Result.x, data=combineAll)
xyplot(totalchlcm ~ Days | Result, data=combineExtraction)

xyplot(chlpsym ~ Days.x | Genotype, data=combineAll)
xyplot(totalcellcm ~ Days | Genotype, data=combineCells)


#chl per sym 
plot(chlpsym ~ Days.x, data=combineAll[which(combineAll$Result.x=="not"),], col="green")
points(chlpsym ~ Days.x, data=combineAll[which(combineAll$Result.x=="above"),], col="blue")
points(chlpsym ~ Days.x, data=combineAll[which(combineAll$Result.x=="below"),], col="red")

mod_cpsym <- lm(chlpsym ~ Result.x * Days.x, data=combineAll )
anova(mod_cpsym)
TukeyHSD(aov(mod_cpsym))
cld(lsmeans(mod_cpsym, specs = "Result.x"))

#Boxplot
par(mfrow=c(1,3))


boxplot(chlpsym ~ Days.x, data=combineAll[which(combineAll$Result.x=="above"),], col="green", ylim=c(0, 3e-6), main="Above average genotypes", xlab="Days", ylab="Ug Chl A per symbiont")
boxplot(chlpsym ~ Days.x, data=combineAll[which(combineAll$Result.x=="not"),], col="yellow", ylim=c(0, 3e-6), main="Not different genotypes", xlab="Days", ylab="Ug Chl A per symbiont")
boxplot(chlpsym ~ Days.x, data=combineAll[which(combineAll$Result.x=="below"),], col="red", ylim=c(0, 3e-6), main="Below average genotypes", xlab="Days", ylab="Ug Chl A per symbiont")
par(mfrow=c(1,1))

#Run analysis on Day 33 differences Day 32 iPam
symDayThT <- combineCells[which(combineCells$Days==33),]
dayThT_symlm <- lm(totalcellcm ~ Result, data = symDayThT)
dayThT_symlmGeno <- lm(totalcellcm ~ Genotype, data = symDayThT)
anova(dayThT_symlm)
TukeyHSD(aov(dayThT_symlm))
cld(lsmeans(dayThT_symlm, "Result"))

anova(dayThT_symlmGeno)
TukeyHSD(aov(dayThT_symlmGeno))
cld(lsmeans(dayThT_symlmGeno, "Genotype"))

#Run analysis on Day 33 differences Day 32 iPam

chlsymDayThT <- combineAll[which(combineAll$Days.x==33),]
dayThT_chlsymlm <- lm(chlpsym ~ Result.x, data = chlsymDayThT)
dayThT_chlsymlmGeno <- lm(chlpsym ~ Genotype, data = chlsymDayThT)
anova(dayThT_chlsymlm)
TukeyHSD(aov(dayThT_chlsymlm))
cld(lsmeans(dayThT_chlsymlm, "Result.x"))

anova(dayThT_chlsymlmGeno)
TukeyHSD(aov(dayThT_chlsymlmGeno))
cld(lsmeans(dayThT_chlsymlmGeno, "Genotype"))

#gamm model by result NOT DONE
combineCells$logcell <- log10(combineCells$totalcellcm)

mod_chlsymgamResult <- gamm4(chlpsym ~ Result.x + s(Days.x, k=5, by=Result.x), random=~(1|FragID.x)+(1|Genotype), data=combineAll)
modMean_gamChl <- gamm4(totalchlcm ~ s(Days, k=5), random=~(1|Genotype), data=combineExtraction)
mod_symgam <- gamm4(totalcellcm ~ Result + s(Days, k=5, by=Result), random=~(1|Genotype)+(1|FragID.x), data=combineCells)
mod_symGenogam <- gamm4(totalcellcm ~ Genotype + s(Days, k=5, by=Genotype), random=~(1|FragID.x), data=combineCells)

data_symResult <- expand.grid(Days=seq(0, 40), Result=levels(combineCells$Result))
data_symGeno <- expand.grid(Days=seq(0, 40), Genotype=levels(combineCells$Genotype))
data_chlpsymResult <- expand.grid(Days=seq(0, 40), Result=levels(combineAll$Result.x))

data_chlpsymResult$fit <-predict(mod_chlsymgamResult$gam, data_chlpsymResult)
data_symResult$fit <-predict(mod_symgam$gam, data_symResult) 
data_symGeno$fit <-predict(mod_symGenogam$gam, data_symGeno) 
xyplot(fit ~ Days |Genotype, data=data_symGeno, type="l")
#For fitted model of all genos
set.seed(563) 
Rbeta <- mvrnorm(n = 10000, coef(mod_symgam$gam), vcov(mod_symgam$gam))
Xp <- predict(mod_symgam$gam, newdata = data_symResult, type = "lpmatrix")
sim_symres <- Xp %*% t(Rbeta)


#For mean not confinced this is necassary... for the derivs
set.seed(543)
Rbetageno <- mvrnorm(n = 10000, coef(mod_symGenogam$gam), vcov(mod_symGenogam$gam))
Xpgeno <- predict(mod_symGenogam$gam, newdata = data_symGeno, type = "lpmatrix")
simsym_geno <- Xpgeno %*% t(Rbetageno)

#For mean not confinced this is necassary... for the derivs
set.seed(729)
RbetaMean <- mvrnorm(n = 10000, coef(modResult_gam$gam), vcov(modResult_gam$gam))
XpMean <- predict(modResult_gam$gam, newdata = data_symGeno, type = "lpmatrix")
simchl_Result <- XpMean %*% t(RbetaMean)


# Extract 90% confidence intervals from simulation
data_symResult$lci <- apply(sim_symres, 1, quantile, 0.05)
data_symResult$uci <- apply(sim_symres, 1, quantile, 0.95)

data_symGeno$lci<- apply(simsym_geno, 1, quantile, 0.05)
data_symGeno$uci <- apply(simsym_geno, 1, quantile, 0.95)

plotfitssymres <- function(data, res, title, key) {
  df <- droplevels(subset(data, Result %in% res))
  getColor <- colorRampPalette(brewer.pal(length(res), "Dark2"))
  colors <- getColor(length(res))
  dfsplit <- split(df, f=df$Result)
  plot(NA, xlim=c(0, 40), ylim=c(0, 5.2e6), xaxs="i", yaxs="i", xlab = "Days", ylab = "Symbionts per cm^2", main=title)
  if(key==T) {
    par(new=T)
    plot(temp$Projected ~ temp$Days, type = "l", ylim=c(25, 35), xaxs="i", axes=F, lwd = 2, xlim=c(0,39), xlab="", ylab="")
    lines(temp$Thres~temp$Days, col = 'orange', lty = 2)
    lines(averageTemp$mean ~ averageTemp$Days, lwd = (.5), xlim=c(0, 39))
    lines(averageTemp$smooth ~ averageTemp$Days, col = 'red')
    axis(side=4)
    mtext("Tempearture (oC)", side=4)
    par(new=T)
    plot(NA, xlim=c(0, 40), ylim=c(0, 10), xaxs="i", yaxs="i", xlab = "Days", ylab = "Chl mg per cm^2", main=title)
  }
  for (i in 1:length(dfsplit)) {
    with(dfsplit[[i]], {
      lines(Days, fit)
      addpoly(Days, lci, uci, col=alpha(colors[i], 0.3))
    })
  }
  
  legend("bottom", legend=levels(df$Result), pch=15, bty="n",
         col=alpha(colors, 0.4), pt.cex=1.5,
         y.intersp=0.3)
}


plotfitssymres(data_symResult, c("above", "not", "below"), "Categories of symbiont density", F)

for(res in data_symResult$Result) {
  sdata <- NA
  sdata <- data_symResult[which(data_symResult$Result==res),]
  for(day in sdata$Days) {
    sdata[which(sdata$Days ==day), "rsym"] <- sdata[which(sdata$Days==day), "fit"] / sdata[which(sdata$Days==0), "fit"]
    sdata[which(sdata$Days ==day), "rlci"] <- sdata[which(sdata$Days==day), "lci"] / sdata[which(sdata$Days==0), "fit"]
    sdata[which(sdata$Days ==day), "ruci"] <- sdata[which(sdata$Days==day), "uci"] / sdata[which(sdata$Days==0), "fit"]
  }
  data_symResult[which(data_symResult$Result==res),"rsym"] <- sdata[, "rsym"]
  data_symResult[which(data_symResult$Result==res), "rlci"] <- sdata[, "rlci"]
  data_symResult[which(data_symResult$Result==res), "ruci"] <- sdata[, "ruci"]
}
#Normal sym fitted model
plotfitssymresnormal <- function(data, res, title, key) {
  df <- droplevels(subset(data, Result %in% res))
  getColor <- colorRampPalette(brewer.pal(length(res), "Dark2"))
  colors <- getColor(length(res))
  dfsplit <- split(df, f=df$Result)
  plot(NA, xlim=c(0, 40), ylim=c(0, 2), xaxs="i", yaxs="i", xlab = "Days", ylab = "Symbionts per cm^2", main=title)
  if(key==T) {
    par(new=T)
    plot(temp$Projected ~ temp$Days, type = "l", ylim=c(25, 35), xaxs="i", axes=F, lwd = 2, xlim=c(0,39), xlab="", ylab="")
    lines(temp$Thres~temp$Days, col = 'orange', lty = 2)
    lines(averageTemp$mean ~ averageTemp$Days, lwd = (.5), xlim=c(0, 39))
    lines(averageTemp$smooth ~ averageTemp$Days, col = 'red')
    axis(side=4)
    mtext("Tempearture (oC)", side=4)
    par(new=T)
    plot(NA, xlim=c(0, 40), ylim=c(0, 10), xaxs="i", yaxs="i", xlab = "Days", ylab = "Chl mg per cm^2", main=title)
  }
  for (i in 1:length(dfsplit)) {
    with(dfsplit[[i]], {
      lines(Days, rsym)
      addpoly(Days, rlci, ruci, col=alpha(colors[i], 0.3))
    })
  }
  
  legend("bottom", legend=levels(df$Result), pch=15, bty="n",
         col=alpha(colors, 0.4), pt.cex=1.5,
         y.intersp=0.3)
}

plotfitssymresnormal(data_symResult, c("above", "not", "below"), "Categories of symbiont density", F)

xyplot(fit ~ Days | Result, data=data_symResult,type="l")
xyplot(fit ~ Days | Genotype, data=data_symGeno,type="l")

#Check for cell counts done
View(table(combineCells$TimePoint))
View(table(combineCells$Genotype))

#Mean se plots for chl per sym

#Mean SE plots for sym

above_chlsym <- combineAll[which(combineAll$Result.x=="above"),]
above_chlsym$Days <- as.integer(above_chlsym$DateSampled.x-as.POSIXct("2017-05-22"))
abovechlsym_summary <- ddply(above_chlsym, c("Days"), summarise, N    = length(chlpsym),
                          mean = mean(chlpsym),
                          sd   = sd(chlpsym),
                          se   = sd / sqrt(N)
)
mod_abovechlpsymsum <- lm(mean ~ Days, data=abovechlsym_summary)
plot(mean ~ Days, abovechlsym_summary)
abline(mod_abovesymsum)
arrows(abovesym_summary$Days, abovesym_summary$mean-abovesym_summary$se, abovesym_summary$Days, abovesym_summary$mean+abovesym_summary$se, length=0.05, angle=90, code=3)

#Below
below_sym <- combineCells[which(combineCells$Result=="below"),]
below_sym$Days <- as.integer(below_sym$DateSampled-as.POSIXct("2017-05-22"))
belowsym_summary <- ddply(below_sym, c("Days"), summarise, N    = length(totalcellcm),
                          mean = mean(totalcellcm),
                          sd   = sd(totalcellcm),
                          se   = sd / sqrt(N)
)
mod_belowsymsum <- lm(mean ~ Days, data=belowsym_summary)
plot(mean ~ Days, belowsym_summary)
abline(mod_belowsymsum)
arrows(belowsym_summary$Days, belowsym_summary$mean-belowsym_summary$se, belowsym_summary$Days, belowsym_summary$mean+belowsym_summary$se, length=0.05, angle=90, code=3)
#Not
not_sym <- combineCells[which(combineCells$Result=="not"),]
not_sym$Days <- as.integer(not_sym$DateSampled-as.POSIXct("2017-05-22"))
notsym_summary <- ddply(not_sym, c("Days"), summarise, N    = length(totalcellcm),
                        mean = mean(totalcellcm),
                        sd   = sd(totalcellcm),
                        se   = sd / sqrt(N)
)
mod_notsymsum <- lm(mean ~ Days, data=notsym_summary)
plot(mean ~ Days, notsym_summary)
abline(mod_notsymsum)
arrows(notsym_summary$Days, notsym_summary$mean-notsym_summary$se, notsym_summary$Days, notsym_summary$mean+notsym_summary$se, length=0.05, angle=90, code=3)

#Plot
par(mfrow = c(1, 3))
plot(mean ~ Days, data=abovesym_summary, ylim=c(0, 1.3E7), ylab="Symbionts per cm^2", xlab="Days", main="Above average genotypes", yaxs="i", col=alpha(colors[1], .9), pch=19)
#abline(mod_abovesymsum)
arrows(abovesym_summary$Days, abovesym_summary$mean-abovesym_summary$se, abovesym_summary$Days, abovesym_summary$mean+abovesym_summary$se, length=0.05, angle=90, code=3)
#text(x = 32, y = 8, labels = expression("Adj. R"^2*"=0.7095"))
plot(mean ~ Days, notsym_summary, ylim=c(0, 1.3E7), ylab="Symbionts per cm^2", xlab="Days", main="Not different genotypes", yaxs="i", col=alpha(colors[2], .9), pch=19)
abline(mod_notsymsum)
arrows(notsym_summary$Days, notsym_summary$mean-notsym_summary$se, notsym_summary$Days, notsym_summary$mean+notsym_summary$se, length=0.05, angle=90, code=3)
#text(x = 32, y = 8, labels = expression("Adj. R"^2*"=0.4903"))

plot(mean ~ Days, belowsym_summary, ylim=c(0, 1.3E7), ylab="Symbionts per cm^2", xlab="Days", main="Below average genotypes", yaxs="i", col=alpha(colors[3], .9), pch=19)
abline(mod_belowsymsum)
arrows(belowsym_summary$Days, belowsym_summary$mean-belowsym_summary$se, belowsym_summary$Days, belowsym_summary$mean+belowsym_summary$se, length=0.05, angle=90, code=3)
#text(x = 32, y = 8, labels = expression("Adj. R"^2*"=0.8462"))
par(mfrow = c(1, 1))


# Test if slopes are differnt
chlsym_summary <- ddply(combineAll, c("Days.x", "Result.x"), summarise, N    = length(chlpsym),
                     mean = mean(chlpsym),
                     sd   = sd(chlpsym),
                     se   = sd / sqrt(N)
)
xyplot(mean ~ Days.x | Result.x, chlsym_summary, type=c("p", "r"))
mod_chlsym_summary <- lm(chlpsym ~ Days.x *Result.x, data=combineAll)
anova(mod_chlsym_summary)
chlsym_lst <- lstrends(mod_chlsym_summary, "Result.x", var="Days.x")
pairs(chlsym_lst)

#Cells to chl plot
plot(logcell ~ totalchlcm, data=combineAll[which(combineAll$Result.x=="above"),], col=alpha("blue", .5), pch=19)
points(logcell ~ totalchlcm, data=combineAll[which(combineAll$Result.x=="not"),], col=alpha("green", .5), pch=19)
points(logcell ~ totalchlcm, data=combineAll[which(combineAll$Result.x=="below"),], col=alpha("red", .5), pch=19)
combineAll$logchl <- log10(combineAll$totalchlcm)
combineAll$logcell <- log10(combineAll$totalcellcm)

plot(totalcellcm ~ totalchlcm, data=combineAll[which(combineAll$Result.x=="above"),], col=alpha("blue", .5), pch=19)
points(totalcellcm ~ totalchlcm, data=combineAll[which(combineAll$Result.x=="not"),], col=alpha("green", .5), pch=19)
points(totalcellcm ~ totalchlcm, data=combineAll[which(combineAll$Result.x=="below"),], col=alpha("red", .5), pch=19)


mod_cellsym <- lm(totalchlcm ~ totalcellcm, data=combineAll)
abline(mod_cellsym)
anova(mod_cellsym)
summary(mod_cellsym)
xyplot(totalcellcm ~ totalchlcm |Genotype, data=combineAll, type=c("p", "r"))
xyplot(totalcellcm ~ Days.x |Genotype, data=combineAll, type=c("p", "r"))
xyplot(totalchlcm ~ Days.x |Genotype, data=combineAll, type=c("p", "r"))
