# Load libraries
library(lattice)
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
library(effects)
library(car)
#library(PMCMR)
#Test for differences in variance
#http://www.sthda.com/english/wiki/compare-multiple-sample-variances-in-r
options(stringsAsFactors = FALSE)

genopastT <- read.csv("Data/TemperatureProjections/GenoCollectwTproj.csv")
genopastT <- merge(compareGeno, genopastT, by=c("Genotype"))
squarepastT <- read.csv("Data/TemperatureProjections/SquareProjections.csv")

#MMM Square
plot(MMM~Lat, data=squarepastT, pch=19)
modSquareMMM <- lm(MMM~Lat, data=squarepastT)
abline(modSquareMMM)
summary(modSquareMMM)

#Trend Century
plot(TrendCentury~Lat, data=squarepastT, main="oC per Century Trend", pch=19)
modTrendCent <- lm(TrendCentury~Lat, data=squarepastT)
abline(modTrendCent)
summary(modTrendCent)

#Amp Square
plot(AmpaCycle~Lat, data=squarepastT, pch=19, main="Amplitude of Annual Cycle")
modAmp <- lm(AmpaCycle~Lat, data=squarepastT)
abline(modAmp)
summary(modAmp)

#Warm Season Trend
plot(WarmSeasonTrend~Lat, data=squarepastT, pch=19, main="Warmest Season Trend in oC per Decade ")
modWST <- lm(WarmSeasonTrend~Lat, data=squarepastT)
abline(modWST)
summary(modWST)

#Warm Month Trend
plot(WarmMonthTrend~Lat, data=squarepastT, pch=19, main="Warmest Month Trend oC per Decade")
modWMT <- lm(WarmMonthTrend~Lat, data=squarepastT)
abline(modWMT)
summary(modWMT)

#DHW Zero
plot(DHWZero~Lat, data=squarepastT, pch=19)
modDHWZ <- lm(DHWZero~Lat, data=squarepastT)
abline(modDHWZ)
summary(modDHWZ)

#DHW Four
plot(DHWFour~Lat, data=squarepastT, pch=19, main="DHW > 4")
modDHWF <- lm(DHWFour~Lat, data=squarepastT)
abline(modDHWF)
summary(modDHWF)

#DHW Eight
plot(DHWEight~Lat, data=squarepastT, pch=19)
modDHWE <- lm(DHWEight~Lat, data=squarepastT)
abline(modDHWE)
summary(modDHWE)

#Factor Genos and result
#str(genopastT)
genopastT$Result<-as.factor(genopastT$Result)
genopastT$Result <- factor(genopastT$Result, levels=c("above", "not", "below"))

#Genos by MMM
plot(MMM~Result, data=genopastT)
modMMM_Result <- lm(MMM~Result, data=genopastT)
anova(modMMM_Result)

plot(MMM~Lat, data=genopastT[which(genopastT$Result=="above"),], col=alpha(colors[1], 0.5), ylim=c(29.65,30.2), pch=19)
points(MMM~Lat, data=genopastT[which(genopastT$Result=="not"),], col=alpha(colors[2], 0.5), pch=19)
points(MMM~Lat, data=genopastT[which(genopastT$Result=="below"),], col=alpha(colors[3], 0.5), pch=19)
points(mean~lat, data=MMM_summary, pch=19)
legend("topright", legend = c("above", "not", "below", "mean"), pch = 19, col=c(alpha(colors[1], .5),alpha(colors[2], .5), alpha(colors[3], .5), "black"))
MMM_summary <- ddply(genopastT, c("Result"), summarise, N    = length(MMM),
                       mean = mean(MMM),
                       sd   = sd(MMM),
                       se   = sd / sqrt(N),
                       lat  = mean(Lat)
)
MMM_summary$number <- seq.int(nrow(MMM_summary))

plot(mean ~ number, data=MMM_summary, xlim=c(1,3), xaxt="n", xlab="Result", ylab="MMM", pch=19, ylim=c(29.2, 30.3))
arrows(MMM_summary$number, MMM_summary$mean-MMM_summary$se, MMM_summary$number, MMM_summary$mean+MMM_summary$se, length=0.05, angle=90, code=3)
axis(side=1,at=lat_summary$number,labels=lat_summary$Result)



shapiro.test(genopastT$MMM)
bartlett.test(MMM~Result, data=genopastT)
leveneTest(MMM~Result, data=genopastT)
fligner.test(MMM~Result, data=genopastT)
#DHW Four
plot(DHWFour~Result, data=genopastT)
mod_DHWF <- lm(DHWFour~Result, data=genopastT)
anova(mod_DHWF)

DHWF_summary <- ddply(genopastT, c("Result"), summarise, N    = length(DHWFour),
                     mean = mean(DHWFour),
                     sd   = sd(DHWFour),
                     se   = sd / sqrt(N),
                     lat  = mean(Lat)
)
DHWF_summary$number <- seq.int(nrow(DHWF_summary))
plot(mean ~ number, data=DHWF_summary, xlim=c(1,3), xaxt="n", xlab="Result", ylab="DHW > 4", pch=19, ylim=c(1.5, 2.2))
arrows(DHWF_summary$number, DHWF_summary$mean-DHWF_summary$se, DHWF_summary$number, DHWF_summary$mean+DHWF_summary$se, length=0.05, angle=90, code=3)
axis(side=1,at=lat_summary$number,labels=lat_summary$Result)

shapiro.test(genopastT$DHWFour)

bartlett.test(DHWFour~Result, data=genopastT)
leveneTest(DHWFour~Result, data=genopastT)
fligner.test(DHWFour~Result, data=genopastT)


#Warm Season Trend Geno
plot(WarmSeasonTrend~Result, data=genopastT)
modWST_geno <- lm(WarmSeasonTrend~Lat, data=squarepastT)
anova(modWST_geno)

wst_summary <- ddply(genopastT, c("Result"), summarise, N    = length(WarmSeasonTrend),
                       mean = mean(WarmSeasonTrend),
                       sd   = sd(WarmSeasonTrend),
                       se   = sd / sqrt(N),
                       lat  = mean(Lat)
)
wst_summary$number <- seq.int(nrow(wst_summary))
plot(mean ~ number, data=wst_summary, xlim=c(1,3), xaxt="n", xlab="Result", ylab="WST oC", pch=19, ylim=c(.17, .19))
arrows(wst_summary$number, wst_summary$mean-wst_summary$se, wst_summary$number, wst_summary$mean+wst_summary$se, length=0.05, angle=90, code=3)
axis(side=1,at=lat_summary$number,labels=lat_summary$Result)

shapiro.test(genopastT$WarmSeasonTrend)

bartlett.test(WarmSeasonTrend~Result, data=genopastT)
leveneTest(WarmSeasonTrend~Result, data=genopastT)
fligner.test(WarmSeasonTrend~Result, data=genopastT)

#TrendCentury Geno
mean(genopastT$TrendCentury)
plot(TrendCentury~Result, data=genopastT)
modTrend_geno <- lm(TrendCentury~Lat, data=squarepastT)
anova(modTrend_geno)

trend_summary <- ddply(genopastT, c("Result"), summarise, N    = length(TrendCentury),
                     mean = mean(TrendCentury),
                     sd   = sd(TrendCentury),
                     se   = sd / sqrt(N),
                     lat  = mean(Lat)
)
trend_summary$number <- seq.int(nrow(trend_summary))
plot(mean ~ number, data=trend_summary, xlim=c(1,3), xaxt="n", xlab="Result", ylab="Trend oC Century", pch=19, ylim=c(1.6, 2.4))
arrows(trend_summary$number, trend_summary$mean-trend_summary$se, trend_summary$number, trend_summary$mean+trend_summary$se, length=0.05, angle=90, code=3)
axis(side=1,at=lat_summary$number,labels=lat_summary$Result)

plot(TrendCentury~Lat, data=genopastT[which(genopastT$Result=="above"),], col=alpha(colors[1], 0.5), pch=19, ylim=c(1.4, 2.3))
points(TrendCentury~Lat, data=genopastT[which(genopastT$Result=="not"),], col=alpha(colors[2], 0.5), pch=19)
points(TrendCentury~Lat, data=genopastT[which(genopastT$Result=="below"),], col=alpha(colors[3], 0.5), pch=19)


shapiro.test(genopastT$TrendCentury)

bartlett.test(TrendCentury~Result, data=genopastT)
leveneTest(TrendCentury~Result, data=genopastT)
fligner.test(TrendCentury~Result, data=genopastT)
#Merge with rFvFv
temprfvfm <- merge(dayThT, genopastT, by=c("Genotype"))

#MMM w SD!! NOT SE
plot(MMM~rfvfm, data=temprfvfm, main="MMM vs rFvFm")
points(MMM~rfvfm, data=temprfvfm[which(temprfvfm$Result=="above"),], col=alpha(colors[1], 0.5), pch=19)
points(MMM~rfvfm, data=temprfvfm[which(temprfvfm$Result=="not"),], col=alpha(colors[2], 0.5), pch=19)
points(MMM~rfvfm, data=temprfvfm[which(temprfvfm$Result=="below"),], col=alpha(colors[3], 0.5), pch=19)
points(mean~rfvfm, data=MMM_summary_fvfm, pch=19)
arrows(MMM_summary_fvfm$rfvfm, MMM_summary_fvfm$mean-MMM_summary_fvfm$sd, MMM_summary_fvfm$rfvfm, MMM_summary_fvfm$mean+MMM_summary_fvfm$sd, length=0.05, angle=90, code=3)
legend("bottomright", legend = c("above", "not", "below", "mean"), pch = 19, col=c(alpha(colors[1], .5),alpha(colors[2], .5), alpha(colors[3], .5), "black"))

MMM_summary_fvfm <- ddply(temprfvfm, c("Result"), summarise, N    = length(MMM),
                       mean = mean(MMM),
                       sd   = sd(MMM),
                       se   = sd / sqrt(N),
                       rfvfm  = mean(rfvfm)
)

mod_rfvfmMMM<- lm(MMM~rfvfm, data=temprfvfm)
anova(mod_rfvfmMMM)
abline(mod_rfvfmMMM)
xyplot(MMM~Lat|Result, data=temprfvfm)
#Trend plot with SD!!
plot(TrendCentury~rfvfm, data=temprfvfm, ylim=c(1.0, 3.5), ylab="Trend (oC/ Century)", main="Century Trend vs. rfvfm")
points(TrendCentury~rfvfm, data=temprfvfm[which(temprfvfm$Result=="above"),], col=alpha(colors[1], 0.5), pch=19)
points(TrendCentury~rfvfm, data=temprfvfm[which(temprfvfm$Result=="not"),], col=alpha(colors[2], 0.5), pch=19)
points(TrendCentury~rfvfm, data=temprfvfm[which(temprfvfm$Result=="below"),], col=alpha(colors[3], 0.5), pch=19)
points(mean~rfvfm, data=trend_summary_fvfm, pch=19)
arrows(trend_summary_fvfm$rfvfm, trend_summary_fvfm$mean-trend_summary_fvfm$sd, trend_summary_fvfm$rfvfm, trend_summary_fvfm$mean+trend_summary_fvfm$sd, length=0.05, angle=90, code=3)
legend("topright", legend = c("above", "not", "below", "mean"), pch = 19, col=c(alpha(colors[1], .5),alpha(colors[2], .5), alpha(colors[3], .5), "black"))
trend_summary_fvfm <- ddply(temprfvfm, c("Result"), summarise, N    = length(TrendCentury),
                          mean = mean(TrendCentury),
                          sd   = sd(TrendCentury),
                          se   = sd / sqrt(N),
                          rfvfm  = mean(rfvfm)
)


mod_rfvfmTrend<- lm(TrendCentury~rfvfm, data=temprfvfm)
anova(mod_rfvfmTrend)
abline(mod_rfvfmTrend)

#MMM vs Trend
xyplot(TrendCentury~MMM | Result, data=genopastT, type="p", ylab="Trend oC/Century", main="Trend vs MMM")
xyplot(TrendCentury~MMM, groups= Result, data=genopastT, type="p", pch=19)
plot(TrendCentury~MMM, data=genopastT, ylab="Trend oC/Century", main="Trend vs MMM", ylim=c(1, 4), xlim=c(29, 31))
points(TrendCentury~MMM, data=temprfvfm[which(temprfvfm$Result=="above"),], col=alpha(colors[1], 0.5), pch=19)
points(TrendCentury~MMM, data=temprfvfm[which(temprfvfm$Result=="not"),], col=alpha(colors[2], 0.5), pch=19)
points(TrendCentury~MMM, data=temprfvfm[which(temprfvfm$Result=="below"),], col=alpha(colors[3], 0.5), pch=19)
points(meanTrend~meanMMM, data=trendmmm_summary[which(trendmmm_summary=="above"),], pch=19)
arrows(trendmmm_summary[which(trendmmm_summary=="above"), "meanMMM"], trendmmm_summary[which(trendmmm_summary=="above"), "meanTrend"]-trendmmm_summary[which(trendmmm_summary=="above"), "sdTrend"], trendmmm_summary[which(trendmmm_summary=="above"), "meanMMM"], trendmmm_summary[which(trendmmm_summary=="above"), "meanTrend"]+trendmmm_summary[which(trendmmm_summary=="above"), "sdTrend"], length=0.05, angle=90, code=3)
arrows(trendmmm_summary[which(trendmmm_summary=="above"), "meanMMM"]-trendmmm_summary[which(trendmmm_summary=="above"), "sdMMM"], trendmmm_summary[which(trendmmm_summary=="above"), "meanTrend"], trendmmm_summary[which(trendmmm_summary=="above"), "meanMMM"]+trendmmm_summary[which(trendmmm_summary=="above"), "sdMMM"], trendmmm_summary[which(trendmmm_summary=="above"), "meanTrend"], length=0.05, angle=90, code=3)
legend("topright", legend = c("above", "not", "below", "mean Above"), pch = 19, col=c(alpha(colors[1], .5),alpha(colors[2], .5), alpha(colors[3], .5), "black"))

trendmmm_summary <- ddply(temprfvfm, c("Result"), summarise, N    = length(TrendCentury),
                            meanTrend = mean(TrendCentury),
                            sdTrend   = sd(TrendCentury),
                            seTrend   = sdTrend / sqrt(N),
                            meanMMM = mean(MMM),
                            sdMMM   = sd(MMM),
                            seMMM   = sdMMM / sqrt(N)
)


mod_mmmtrend<- lm(TrendCentury~MMM, data=genopastT)
anova(mod_mmmtrend)
abline(mod_mmmtrend)
summary(mod_mmmtrend)
mmmTrend_lst <- lstrends(mod_mmmtrend, "Result", var="MMM")
pairs(mmmTrend_lst)
cld(lstrends(mod_mmmtrend, "Result", var="MMM"))
TukeyHSD(aov(mod_sym_summary))
cld(lsmeans(mod_sym_summary, "Result"))

#DHW Four
plot(DHWFour~rfvfm, data=temprfvfm)
mod_rfvfmDHWF<- lm(DHWFour~rfvfm, data=temprfvfm)
anova(mod_rfvfmDHWF)
abline(mod_rfvfmDHWF)

#Warm Season
plot(WarmSeasonTrend~rfvfm, data=temprfvfm)
mod_rfvfmWST<- lm(WarmSeasonTrend~rfvfm, data=temprfvfm)
anova(mod_rfvfmWST)
abline(mod_rfvfmWST)

#Warm Month
plot(WarmMonthTrend~rfvfm, data=temprfvfm)
mod_rfvfmWMT<- lm(WarmMonthTrend~rfvfm, data=temprfvfm)
anova(mod_rfvfmWMT)
abline(mod_rfvfmWMT)


#Lat vs rfvfm
plot(Lat~rfvfm, data=temprfvfm)
points(Lat~rfvfm, data=temprfvfm[which(temprfvfm$Result=="above"),], col=alpha(colors[1], 0.5), pch=19)
points(Lat~rfvfm, data=temprfvfm[which(temprfvfm$Result=="not"),], col=alpha(colors[2], 0.5), pch=19)
points(Lat~rfvfm, data=temprfvfm[which(temprfvfm$Result=="below"),], col=alpha(colors[3], 0.5), pch=19)

mod_Latrfvfm<-lm(Lat~rfvfm, data=temprfvfm)
abline(mod_Latrfvfm)
anova(mod_Latrfvfm)
plot(Effect("rfvfm", mod_Latrfvfm))



