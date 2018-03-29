#Growth data
library(zoo)
library(multcomp)
library(multcompView)
library(lattice)
library(lsmeans)
library(plyr)
options(stringsAsFactors = FALSE)



growthData <- read.csv("Data/CollectionR/growth_Ford.csv")
fordNames <- read.csv("Data/CollectionR/ford_genonames.csv")
fordNamesnocoop <- subset(fordNames, subset = fordNames$Geno!="1730") 


ford <- unique(growthData[, "Genotype"])
cooper <- subset(all,subset = all$Genotype=="1730")
cooper$Days <- with(cooper, as.numeric(Date - first(Date[order(Date)])))
cooper$normal <- cooper$Y/cooper[which(cooper$Days==0), "Y"]

combineGrowth <- merge(compareGeno, growthData, by = c("Genotype"))
combineGrowth$Result <- factor(combineGrowth$Result, levels = c("above", "not", "below"))
combineGrowth$Site <- as.factor(combineGrowth$Site)
combineGrowth$Result <- as.factor(combineGrowth$Result)
fordOverview <- unique(combineGrowth$Genotype) 
for(i in 1:length(fordNamesnocoop$Geno)) {
  fordNamesnocoop[i, "Result"] <- unique(combineGrowth[which(combineGrowth$Genotype==fordNamesnocoop[i, "Geno"]), "Result"])
}
write.csv(x = fordNamesnocoop, file = "fordResults.csv")


plot(GrowthRate ~ Genotype, data = combineGrowth)
plot(GrowthRate ~ Site, data = combineGrowth)
xyplot(GrowthRate ~ Genotype | Site, data=combineGrowth, type=c("p"), ylab="Growth Rate (cm/day)")

mod_growth <- lm(GrowthRate~Genotype, data=combineGrowth)
anova(mod_growth)
TukeyHSD(aov(mod_growth))
cld(lsmeans(mod_growth, "Genotype"))

mod_site <- lm(GrowthRate~Site, data=combineGrowth)
anova(mod_site)
TukeyHSD(aov(mod_site))
cld(lsmeans(mod_site, "Site"))

plot(GrowthRate ~ Result, data = combineGrowth)
xyplot(GrowthRate ~ Result, data=combineGrowth, type=c("p"))
xyplot(GrowthRate ~ Result | Site, data=combineGrowth, type=c("p", "r"), ylab="Growth Rate (cm/day)")
xyplot(GrowthRate ~ Site | Result, data=combineGrowth, type=c("p"))

mod_res <- lm(GrowthRate~Result, data=combineGrowth)
anova(mod_res)
TukeyHSD(aov(mod_res))
cld(lsmeans(mod_res, "Result"))

mod_resite <- lm(GrowthRate~Result + Site, data=combineGrowth)
anova(mod_resite)
TukeyHSD(aov(mod_resite))
cld(lsmeans(mod_resite, "Result"))

growth_summary <- ddply(combineGrowth, c("Result"), summarise, N    = length(GrowthRate),
                     mean = mean(GrowthRate),
                     sd   = sd(GrowthRate),
                     se   = sd / sqrt(N)
)
growth_summary$number <- seq.int(nrow(growth_summary))
plot(mean ~ number, data=growth_summary, xlim=c(1,3), xaxt="n", xlab="Result", ylab="Growth Rate (cm/day)", ylim=c(0, .035), pch=19)
arrows(growth_summary$number, growth_summary$mean-growth_summary$se, growth_summary$number, growth_summary$mean+growth_summary$se, length=0.05, angle=90, code=3)
axis(side=1,at=growth_summary$number,labels=growth_summary$Result)

#Ford Plots

plotfitsnormal(newdata, ford, "Ford Genotypes", F)
points(normal~Days, data= cooper)
lines(rfvfm~Days, data=newdatamean, lwd=2, col="red")
addpoly(newdatamean$Days, newdatamean$rlci, newdatamean$ruci, col=alpha("red", 0.3))
abline(v=32, col= "red")
text(x= 32, y = 1.0, adj = c(.5), "Statistically Different")

legend("bottomleft", inset=c(.1, 0), legend=fordNamesnocoop$Name, bty="n",
       y.intersp=0.3)

#Correlattion of growth to rfvfm

combineGrowthPredict <- merge(dayThT, growthData, by = c("Genotype"))
combineGrowthPredict$Site <- as.factor(combineGrowth$Site)
xyplot(GrowthRate ~ rfvfm | Site, data=combineGrowthPredict)
plot(GrowthRate~ rfvfm, data=combineGrowthPredict)
basicGrowth <- lm(GrowthRate~ rfvfm + Site, data=combineGrowthPredict)
anova(basicGrowth)
TukeyHSD(aov(basicGrowth))
cld(lsmeans(basicGrowth, "rfvfm"))
abline(basicGrowth)
library(effects)
plot(Effect("rfvfm", basicGrowth, col="red"), main="Growth rate vs. relative Fv/Fm", ylab = "Growth Rate cm/day", xlab="rFv/Fm")

abline(basicGrowth)
modGrowth <- lmer(GrowthRate ~ rfvfm+(1|Site), data=combineGrowthPredict)
plot(modGrowth, specs = "rfvfm")

#Subset data based on sites by nursery
nursery <- subset(combineGrowthPredict,subset = combineGrowthPredict$Site=="Nursery")
plot(GrowthRate~ rfvfm, data=nursery)
nurserylm <- lm(GrowthRate~ rfvfm, data=nursery)
anova(nurserylm)
abline(nurserylm)

nursery <- subset(combineGrowth,subset = combineGrowth$Site=="Nursery")
plot(GrowthRate~ Result, data=nursery)
nurserylm <- lm(GrowthRate~ Result, data=nursery)
anova(nurserylm)
cld(lsmeans(nurserylm, "Result"))
abline(nurserylm)

#Subset data based on sites by stephs
stephs <- subset(combineGrowthPredict,subset = combineGrowthPredict$Site=="Steph's")
plot(GrowthRate~ rfvfm, data=stephs)
stephslm <- lm(GrowthRate~ rfvfm, data=stephs)
anova(stephslm)
abline(stephslm)

nursery <- subset(combineGrowth,subset = combineGrowth$Site=="Steph's")
plot(GrowthRate~ Result, data=nursery)
nurserylm <- lm(GrowthRate~ Result, data=nursery)
anova(nurserylm)
cld(lsmeans(nurserylm, "Result"))
abline(nurserylm)

#Subset data based on sites by inshore
inshore <- subset(combineGrowthPredict,subset = combineGrowthPredict$Site=="Inshore")
plot(GrowthRate~ rfvfm, data=inshore)
inshorelm <- lm(GrowthRate~ rfvfm, data=inshore)
anova(inshorelm)
abline(inshorelm)

nursery <- subset(combineGrowth,subset = combineGrowth$Site=="Inshore")
plot(GrowthRate~ Result, data=nursery)
nurserylm <- lm(GrowthRate~ Result, data=nursery)
anova(nurserylm)
cld(lsmeans(nurserylm, "Result"))
abline(nurserylm)

#Subset data based on sites by struggle bus
strug <- subset(combineGrowthPredict,subset = combineGrowthPredict$Site=="Struggle Bus")
plot(GrowthRate~ rfvfm, data=strug)
struglm <- lm(GrowthRate~ rfvfm, data=strug)
anova(struglm)
abline(struglm)

#Subset data based on sites by Jon's
jon <- subset(combineGrowthPredict,subset = combineGrowthPredict$Site=="Jon's")
plot(GrowthRate~ rfvfm, data=jon)
jonlm <- lm(GrowthRate~ rfvfm, data=jon)
anova(jonlm)
abline(jonlm)

#Subset data based on sites by Miami Beach
mb <- subset(combineGrowthPredict,subset = combineGrowthPredict$Site=="Miami Beach")
plot(GrowthRate~ rfvfm, data=mb)
mblm <- lm(GrowthRate~ rfvfm, data=mb)
anova(mblm)
abline(mblm)

#Subset data based on sites by COoper
coop <- subset(combineGrowthPredict,subset = combineGrowthPredict$Site=="Cooper's")
plot(GrowthRate~ rfvfm, data=coop)
cooplm <- lm(GrowthRate~ rfvfm, data=coop)
anova(cooplm)
abline(cooplm)

#Subset data based on sites by CVFD
cvfd <- subset(combineGrowthPredict,subset = combineGrowthPredict$Site=="CVFD")
plot(GrowthRate~ rfvfm, data=cvfd)
cvfdlm <- lm(GrowthRate~ rfvfm, data=cvfd)
anova(cvfdlm)
abline(cvfdlm)

#Subset data based on sites by Grounding
ground <- subset(combineGrowthPredict,subset = combineGrowthPredict$Site=="Grounding")
plot(GrowthRate~ rfvfm, data=ground)
groundlm <- lm(GrowthRate~ rfvfm, data=ground)
anova(groundlm)
abline(groundlm)

nursery <- subset(combineGrowth,subset = combineGrowth$Site=="Grounding")
plot(GrowthRate~ Result, data=nursery)
nurserylm <- lm(GrowthRate~ Result, data=nursery)
anova(nurserylm)
cld(lsmeans(nurserylm, "Result"))
abline(nurserylm)

#Show graph of 3 designations for Ford genotypes
plotfitsnormal(graphData, c("1731", "1735", "1727"), NA, T)
lines(rfvfm~Days, data=newdatamean, lwd=2, col="red")
abline(v=32, col= "red")
text(x= 32, y = 1.0, adj = c(.5), "Statistically Different")

legend("bottomleft", inset=c(.1, 0), legend=c("Not", "Below", "Above"), bty="n",
       y.intersp=0.3)

#Phil time point graph
sampleInfo <- data.frame(NA)
sampleDays <- c(2, 11, 18, 25, 31, 33, 36, 38, 40)
for (i in 1:length(sampleDays)) {
  sampleInfo[i, "Day"] <- sampleDays[i]
  sampleInfo[i, "TimePoint"] <- paste("T_", i-1)   
}
plotfitsnormal(graphData, c("1731", "1735", "1727"), NA, T)
lines(rfvfm~Days, data=newdatamean, lwd=2, col="red")
abline(v=sampleDays, lwd=2, col="blue")
text(x=sampleInfo$Day, y=0.9, labels = sampleInfo$TimePoint, pos = 2, cex=.7, offset = .1) 

legend(9, .1, legend=c("Not", "Below", "Above"), bty="n",
       y.intersp=0.3)

#For plotting ipam and sampling days
ipamdays <- unique(all.f$Days)
plot(temp$Projected ~ temp$Days, type = "l", ylim=c(25, 35), main="", xaxs="i", lwd = 2, xlim=c(0,40), xlab="Days", ylab="Temperature (oC)", yaxs="i")
lines(temp$Thres~temp$Days, col = 'orange', lty = 2)
lines(averageTemp$mean ~ averageTemp$Days, lwd = (.5), xlim=c(0, 39))
lines(averageTemp$smooth ~ averageTemp$Days, col = 'red')
abline(v=ipamdays, lwd=2, col="blue")
text(x=sampleInfo$Day, y=33, labels = sampleInfo$TimePoint, pos = 2, cex=.7, offset = .1) 

legend("bottomright", legend = c("Projected Temperature", "Mean Temperature", "I-Pam Sample"), lwd = c(2, 2, 2), col = c("black", "red", "blue"))

#Plot normal showing differences 

plotfitsnormal(graphData, c("1721", "1735", "1748"), "Model of normalized photochemical efficiency for 3 genotypes", T)
lines(rfvfm~Days, data=newdatamean, lwd=2, col="red")
abline(v=32, col= "red")
text(x= 32, y = 0.9, "Statistically Different", pos = 4)
legend("bottomleft", inset=c(.5, 0), legend=c("Below Mean", "Above Mean", "Not Significantly Different"), bty="n",
       y.intersp=0.3)
