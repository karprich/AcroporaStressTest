#load packages
library(devtools)
library(plyr)
library(stringr)
library(reshape2)
library(lattice)
library(lsmeans)
library(RColorBrewer)

#load qpcr script
devtools::source_url("https://raw.githubusercontent.com/jrcunning/steponeR/master/steponeR.R")

#load files 
plates<-list.files(path="Data/SampleExtraction/PCR/rdata",pattern="csv$", full.names=TRUE)
Acerv<-steponeR(files=plates,target.ratios=c("A.Ac"),
                fluor.norm=list(A=-.064,Ac=5),
                copy.number=list(A=1,Ac=1),
                ploidy=list(A=1,Ac=2),
                extract=list(A=.813,Ac=.982))
  


result<-Acerv$result



#organize data 

result2<-result[!grepl(pattern="*p5",result$Sample.Name), ]
result3<-result2[!grepl(pattern="*R",result$Sample.Name), ]
result4<-result3[!grepl(pattern="*58C_0",result$Sample.Name), ]
result5<-result4[!grepl(pattern="*55E_0",result$Sample.Name), ]

result3<-result2[!grepl(pattern="*NTC", result2$Sample.Name), ]

result4<-result3[!grepl(pattern="*dilute", result3$Sample.Name), ]

#result5<-result4[!grepl(pattern="*.2", result4$Sample.Name), ]
result6<-result4[which(result4$A.Ac!="NaN"),]
result6<-result6[which(result6$Sample.Name!="58C_0"),]
result6<-result6[which(result6$Sample.Name!="55E_0_2"),]
result6<-result6[which(result6$Sample.Name!="55E_0"),]

#assign columns

result6$TimePoint<-str_extract(result6$Sample.Name,pattern="_.")
result6$TimePoint<-str_sub(result6$TimePoint,start = 2, end = 3)
result6$Genotype<-paste("17", str_extract(result6$Sample.Name,pattern=".."), sep = "")
result6$FragID<-str_extract(result6$Sample.Name,pattern="[A-Z]")

#merge grouping

merged.results<-merge(result6,compareGeno,by=c("Genotype"))



#make the timepoints the right order 
merged.results$TimePoint<-factor(merged.results$TimePoint,levels=c("0","1","2","3","4","5","6","7", "8"))
merged.results$Result<-factor(merged.results$Result,levels=c("above","not","below"))

#take log
merged.results$logSH<-log10(merged.results$A.Ac)

xyplot(logSH ~ TimePoint | Result, data=merged.results)
xyplot(logSH ~ TimePoint | Genotype, data=merged.results)
boxplot(logSH~TimePoint, data=merged.results)
boxplot(logSH~TimePoint, data=merged.results[which(merged.results$Result=="above"),])
boxplot(logSH~TimePoint, data=merged.results[which(merged.results$Result=="below"),])
boxplot(logSH~TimePoint, data=merged.results[which(merged.results$Result=="not"),])

#Merge blasting data
qPCR_Results <- merge(blast, merged.results, by = c("TimePoint", "Genotype")) #With Frag ID as factor. Misses 19 obs
qPCR_Results$Days <- as.integer(qPCR_Results$DateSampled-as.POSIXct("2017-05-22"))
qPCR_Results$Genotype <- as.factor(qPCR_Results$Genotype)
xyplot(logSH ~ Days | Result, data=qPCR_Results)
xyplot(logSH ~ Days | Genotype, data=qPCR_Results)

plot(logSH ~ Days, data=qPCR_Results[which(qPCR_Results$Result=="not"),], col="green")
points(logSH ~ Days, data=qPCR_Results[which(qPCR_Results$Result=="above"),], col="blue")
points(logSH ~ Days, data=qPCR_Results[which(qPCR_Results$Result=="below"),], col="red")

xyplot(A.Ac~ Days | Result, data = qPCR_Results, ylim=c(0,500))

qpcr_mod<- lm(logSH~Days * Result, data=qPCR_Results)
anova(qpcr_mod)
lsmeans(qpcr_mod, specs="Result")

#Try ddply
dataqpcr_summary <- ddply(qPCR_Results, c("Result", "Days"), summarise, N    = length(logSH),
                          mean = mean(logSH),
                          sd   = sd(logSH),
                          se   = sd / sqrt(N)
)

plot(mean~Days, dataqpcr_summary[which(dataqpcr_summary$Result=="above"),], col="blue", ylim=c(-4,4))
points(mean~Days, dataqpcr_summary[which(dataqpcr_summary$Result=="not"),], col="green")
points(mean~Days, dataqpcr_summary[which(dataqpcr_summary$Result=="below"),], col="red")

xyplot(mean~Days | Result, dataqpcr_summary)

#Form relative boxplots of each qpcr
above_qpcr <- qPCR_Results[which(qPCR_Results$Result=="above"),]

above_qpcr$logRel <- above_qpcr$logSH/dataqpcr_summary[which(dataqpcr_summary$Result=="above" & dataqpcr_summary$Days==2), "mean"]
aboveqpcr_summary <- ddply(above_qpcr, c("Days"), summarise, N    = length(logSH),
                           mean = mean(logSH),
                           sd   = sd(logSH),
                           se   = sd / sqrt(N),
                           relmean=mean(logRel),
                           rsd=sd(logRel),
                           rse=rsd/sqrt(N)
)
aboveqpcr_summary$Result <- as.factor("above")
mod_abovePCRsum <- lm(mean ~ Days, aboveqpcr_summary)
plot(mean ~ Days, aboveqpcr_summary, ylim=c(-5, 5), ylab="log(S:H)", xlab="Days", main="Above average genotypes", yaxs="i", col=alpha(colors[1], .9), pch=19)
abline(mod_abovePCRsum)
arrows(aboveqpcr_summary$Days, aboveqpcr_summary$mean-aboveqpcr_summary$se, aboveqpcr_summary$Days, aboveqpcr_summary$mean+aboveqpcr_summary$se, length=0.05, angle=90, code=3)
#Relative above
mod_abovePCRsumrel <- lm(relmean ~ Days, aboveqpcr_summary)
plot(relmean ~ Days, aboveqpcr_summary, ylim=c(-2, 2), ylab="log(S:H)", xlab="Days", main="Above average genotypes", yaxs="i", col=alpha(colors[1], .9), pch=19)
abline(mod_abovePCRsumrel)
arrows(aboveqpcr_summary$Days, aboveqpcr_summary$relmean-aboveqpcr_summary$rse, aboveqpcr_summary$Days, aboveqpcr_summary$relmean+aboveqpcr_summary$rse, length=0.05, angle=90, code=3)


#Below summary
below_qpcr <- qPCR_Results[which(qPCR_Results$Result=="below"),]
below_qpcr$Result <- "below"
below_qpcr$logRel <- below_qpcr$logSH/dataqpcr_summary[which(dataqpcr_summary$Result=="below" & dataqpcr_summary$Days==2), "mean"]
belowqpcr_summary <- ddply(below_qpcr, c("Days"), summarise, N    = length(logSH),
                           mean = mean(logSH),
                           sd   = sd(logSH),
                           se   = sd / sqrt(N),
                           relmean=mean(logRel),
                           rsd=sd(logRel),
                           rse=rsd/sqrt(N)
)
belowqpcr_summary$Result <- as.factor("below")
mod_belowPCRsum <- lm(mean ~ Days, belowqpcr_summary)
plot(mean ~ Days, belowqpcr_summary)
abline(mod_belowPCRsum)
arrows(belowqpcr_summary$Days, belowqpcr_summary$mean-belowqpcr_summary$se, belowqpcr_summary$Days, belowqpcr_summary$mean+belowqpcr_summary$se, length=0.05, angle=90, code=3)

mod_belowPCRsumrel <- lm(relmean ~ Days, belowqpcr_summary)
points(relmean ~ Days, belowqpcr_summary, ylim=c(-1.5, 1.5), col="red")
abline(mod_belowPCRsumrel)
arrows(belowqpcr_summary$Days, belowqpcr_summary$relmean-belowqpcr_summary$rse, belowqpcr_summary$Days, belowqpcr_summary$relmean+belowqpcr_summary$rse, length=0.05, angle=90, code=3)


#Not
not_qpcr <- qPCR_Results[which(qPCR_Results$Result=="not"),]
not_qpcr$Result <- as.factor("not")
not_qpcr$logRel <- not_qpcr$logSH/dataqpcr_summary[which(dataqpcr_summary$Result=="not" & dataqpcr_summary$Days==2), "mean"]

notqpcr_summary <- ddply(not_qpcr, c("Days"), summarise, N    = length(logSH),
                         mean = mean(logSH),
                         sd   = sd(logSH),
                         se   = sd / sqrt(N),
                         relmean=mean(logRel),
                         rsd=sd(logRel),
                         rse=rsd/sqrt(N)
)
notqpcr_summary$Result <- as.factor("not")
mod_notPCRsum <- lm(mean ~ Days, notqpcr_summary)
plot(mean ~ Days, notqpcr_summary, ylim=c(-5, 3.5))
abline(mod_notPCRsum)
arrows(notqpcr_summary$Days, notqpcr_summary$mean-notqpcr_summary$se, notqpcr_summary$Days, notqpcr_summary$mean+notqpcr_summary$se, length=0.05, angle=90, code=3)

mod_notPCRsumrel <- lm(relmean ~ Days, notqpcr_summary)
points(relmean ~ Days, notqpcr_summary, ylim=c(-1.5, 1.5))
abline(mod_notPCRsumrel)
arrows(notqpcr_summary$Days, notqpcr_summary$relmean-notqpcr_summary$rse, notqpcr_summary$Days, notqpcr_summary$relmean+notqpcr_summary$rse, length=0.05, angle=90, code=3)



# 3 plot with brewer color pallete
getColor <- colorRampPalette(brewer.pal(5, "Dark2"))
colors <- getColor(5)

#par(mfrow = c(1, 3))
plot(mean ~ Days, data=aboveqpcr_summary, ylim=c(-2.5, 2.5), ylab="log(S:H)", xlab="Days", main="Above average genotypes", yaxs="i", col=alpha(colors[1], .5), pch=19)
#abline(mod_abovePCRsum)
arrows(aboveqpcr_summary$Days, aboveqpcr_summary$mean-aboveqpcr_summary$se, aboveqpcr_summary$Days, aboveqpcr_summary$mean+aboveqpcr_summary$se, length=0.05, angle=90, code=3)
#text(x = 32, y = 2, labels = expression("Adj. R"^2*"=0.4419"))
#text(x=20, y= -1, labels = expression("Slope: -0.0364"))
plot(mean ~ Days, notqpcr_summary, ylim=c(-2.5, 2.5), ylab="log(S:H)", xlab="Days", main="Not different genotypes", yaxs="i", col=alpha(colors[2], .5), pch=19)
#abline(mod_notPCRsum)
arrows(notqpcr_summary$Days, notqpcr_summary$mean-notqpcr_summary$se, notqpcr_summary$Days, notqpcr_summary$mean+notqpcr_summary$se, length=0.05, angle=90, code=3)
#text(x = 32, y = 2, labels = expression("Adj. R"^2*"=0.4918"))
plot(mean ~ Days, belowqpcr_summary, ylim=c(-2.5, 2.5), ylab="log(S:H)", xlab="Days", main="Below average genotypes", yaxs="i", col=alpha(colors[3], .5), pch=19)
abline(mod_belowPCRsum)
arrows(belowqpcr_summary$Days, belowqpcr_summary$mean-belowqpcr_summary$se, belowqpcr_summary$Days, belowqpcr_summary$mean+belowqpcr_summary$se, length=0.05, angle=90, code=3)
text(x = 32, y = 2, labels = expression("Adj. R"^2*"=0.6375"))
#par(mfrow = c(1, 1))

#relative qpcr
plot(relmean ~ Days, data=aboveqpcr_summary, ylim=c(-2.5, 2.5), ylab="log(S:H)", xlab="Days", main="Above average genotypes", yaxs="i", col=alpha(colors[1], .5), pch=19)
arrows(aboveqpcr_summary$Days, aboveqpcr_summary$relmean-aboveqpcr_summary$rse, aboveqpcr_summary$Days, aboveqpcr_summary$relmean+aboveqpcr_summary$rse, length=0.05, angle=90, code=3)
#text(x = 32, y = 2, labels = expression("Adj. R"^2*"=0.4419"))
#text(x=20, y= -1, labels = expression("Slope: -0.0364"))
points(relmean ~ Days, notqpcr_summary, ylim=c(-2.5, 2.5), ylab="log(S:H)", xlab="Days", main="Not different genotypes", yaxs="i", col=alpha(colors[2], .5), pch=19)
#abline(mod_notPCRsum)
arrows(notqpcr_summary$Days, notqpcr_summary$relmean-notqpcr_summary$rse, notqpcr_summary$Days, notqpcr_summary$relmean+notqpcr_summary$rse, length=0.05, angle=90, code=3)
#text(x = 32, y = 2, labels = expression("Adj. R"^2*"=0.4918"))
points(relmean ~ Days, belowqpcr_summary, ylim=c(-2.5, 2.5), ylab="log(S:H)", xlab="Days", main="Below average genotypes", yaxs="i", col=alpha(colors[3], .3), pch=19)

arrows(belowqpcr_summary$Days, belowqpcr_summary$relmean-belowqpcr_summary$rse, belowqpcr_summary$Days, belowqpcr_summary$relmean+belowqpcr_summary$rse, length=0.05, angle=90, code=3)

#Analysis of whole

pcr_summary <- ddply(qPCR_Results, c("Days", "Result"), summarise, N    = length(logSH),
                     mean = mean(logSH),
                     sd   = sd(logSH),
                     se   = sd / sqrt(N)
)
mod_pcr_summary <- lm(logSH ~ Days *Result, data=qPCR_Results)
anova(mod_pcr_summary)
pcr_lst <- lstrends(mod_pcr_summary, "Result", var="Days")
pairs(pcr_lst)
cld(lstrends(mod_pcr_summary, "Result", var="Days"))
TukeyHSD(aov(mod_pcr_summary))
cld(lsmeans(mod_pcr_summary, "Result"))

xyplot(mean~Days, data=pcr_summary, groups = Result)

relcombineqpcr<- rbind(aboveqpcr_summary, notqpcr_summary, belowqpcr_summary)
mod_pcr_summaryrel <- lm(relmean ~ Days *Result, data=relcombineqpcr)
anova(mod_pcr_summaryrel)
pcr_rellst <- lstrends(mod_pcr_summaryrel, "Result", var="Days")
pairs(pcr_rellst)
cld(lstrends(mod_pcr_summaryrel, "Result", var="Days"))
xyplot(relmean~Days, data=relcombineqpcr, groups = Result, pch=19)
#Run GAMM
modMean_pcr <- gamm4(logSH ~ s(Days, k=5), random=~(1|Genotype)+(1|FragID.x), data=qPCR_Results)
modResult_pcr <- gamm4(logSH ~ Result + s(Days, k=5, by=Result), random=~(1|Genotype)+(1|FragID.x), data=qPCR_Results)
modGeno_pcr <- gamm4(logSH ~ Genotype + s(Days, k=5, by=Genotype), random=~(1|FragID.x), data=qPCR_Results)

data_pcrRes <- expand.grid(Days=seq(0, 40), Result=levels(qPCR_Results$Result))
data_pcrMean <- expand.grid(Days=seq(0, 40))
data_pcrGeno<- expand.grid(Days=seq(0, 40), Genotype=levels(qPCR_Results$Genotype))

data_pcrMean$fit <-predict(modMean_pcr$gam, data_pcrMean)
data_pcrRes$fit <-predict(modResult_pcr$gam, data_pcrRes)
data_pcrGeno$fit <- predict(modGeno_pcr$gam, data_pcrGeno)

xyplot(fit~Days, data=data_pcrRes, groups = Result, type="l")
plot(fit ~ Days, data=data_pcrRes[which(data_pcrRes$Result=="above"),], type="l", ylim=c(-3, 3))
lines(fit ~ Days, data=data_pcrRes[which(data_pcrRes$Result=="below"),], col="red")
lines(fit ~ Days, data=data_pcrRes[which(data_pcrRes$Result=="not"),], col="green")
lines(fit~ Days, data=data_pcrMean)
xyplot(fit ~ Days | Genotype, data=data_pcrGeno, type="l")
#For fitted model of res
set.seed(198) 
Rbeta <- mvrnorm(n = 10000, coef(modResult_pcr$gam), vcov(modResult_pcr$gam))
Xp <- predict(modResult_pcr$gam, newdata = data_pcrRes, type = "lpmatrix")
sim_pcr <- Xp %*% t(Rbeta)

#For fitted model of genos
set.seed(516) 
Rbetag <- mvrnorm(n = 10000, coef(modGeno_pcr$gam), vcov(modGeno_pcr$gam))
Xpg <- predict(modGeno_pcr$gam, newdata = data_pcrGeno, type = "lpmatrix")
sim_pcrgeno <- Xpg %*% t(Rbetag)

#For mean not confinced this is necassary... for the derivs
set.seed(897)
RbetaMean <- mvrnorm(n = 10000, coef(modMean_pcr$gam), vcov(modMean_pcr$gam))
XpMean <- predict(modMean_pcr$gam, newdata = data_pcrMean, type = "lpmatrix")
simMean_pcr <- XpMean %*% t(RbetaMean)




# Extract 90% confidence intervals from simulation
data_pcrRes$lci <- apply(sim_pcr, 1, quantile, 0.05)
data_pcrRes$uci <- apply(sim_pcr, 1, quantile, 0.95)
data_pcrMean$lci <- apply(simMean_pcr, 1, quantile, 0.05)
data_pcrMean$uci <- apply(simMean_pcr, 1, quantile, 0.95)
data_pcrGeno$lci <- apply(sim_pcrgeno, 1, quantile, 0.05)
data_pcrGeno$uci <- apply(sim_pcrgeno, 1, quantile, 0.95)


#Plot fits chl function
plotfitspcrres <- function(data, res, title, key) {
  df <- droplevels(subset(data, Result %in% res))
  getColor <- colorRampPalette(brewer.pal(length(res), "Dark2"))
  colors <- getColor(length(res))
  dfsplit <- split(df, f=df$Result)
  plot(NA, xlim=c(0, 40), ylim=range(df$fit), xaxs="i", yaxs="i", xlab = "Days", ylab = "log(S:H)", main=title)
  if(key==T) {
    par(new=T)
    plot(temp$Projected ~ temp$Days, type = "l", ylim=c(25, 35), xaxs="i", axes=F, lwd = 2, xlim=c(0,39), xlab="", ylab="")
    lines(temp$Thres~temp$Days, col = 'orange', lty = 2)
    lines(averageTemp$mean ~ averageTemp$Days, lwd = (.5), xlim=c(0, 39))
    lines(averageTemp$smooth ~ averageTemp$Days, col = 'red')
    axis(side=4)
    mtext("Tempearture (oC)", side=4)
    par(new=T)
    plot(NA, xlim=c(0, 40), ylim=range(df$fit), xaxs="i", yaxs="i", xlab = "Days", ylab = "log(S:H)", main=title)
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

plotfitspcrres(data_pcrRes, c("above", "below", "not"), NA, F)

#Normalize plotting
plotfitspcrresnormal <- function(data, res, title, key) {
  df <- droplevels(subset(data, Result %in% res))
  getColor <- colorRampPalette(brewer.pal(length(res), "Dark2"))
  colors <- getColor(length(res))
  dfsplit <- split(df, f=df$Result)
  plot(NA, xlim=c(0, 40), ylim=range(df$rlog), xaxs="i", yaxs="i", xlab = "Days", ylab = "Relative log(S:H)", main=title)
  if(key==T) {
    par(new=T)
    plot(temp$Projected ~ temp$Days, type = "l", ylim=c(25, 35), xaxs="i", axes=F, lwd = 2, xlim=c(0,39), xlab="", ylab="")
    lines(temp$Thres~temp$Days, col = 'orange', lty = 2)
    lines(averageTemp$mean ~ averageTemp$Days, lwd = (.5), xlim=c(0, 39))
    lines(averageTemp$smooth ~ averageTemp$Days, col = 'red')
    axis(side=4)
    mtext("Tempearture (oC)", side=4)
    par(new=T)
    plot(NA, xlim=c(0, 40), ylim=range(df$rlog), xaxs="i", yaxs="i", xlab = "Days", ylab = "Relative log(S:H)", main=title)
  }
  for (i in 1:length(dfsplit)) {
    with(dfsplit[[i]], {
      lines(Days, rlog)
      addpoly(Days, rlci, ruci, col=alpha(colors[i], 0.3))
    })
  }
  
  legend("bottom", legend=levels(df$Result), pch=15, bty="n",
         col=alpha(colors, 0.4), pt.cex=1.5,
         y.intersp=0.3)
}

plotfitspcrGeno <- function(data, geno, title) {
  df <- droplevels(subset(data, Genotype %in% geno))
  getColor <- colorRampPalette(brewer.pal(length(geno), "Dark2"))
  colors <- getColor(length(geno))
  dfsplit <- split(df, f=df$geno)
  plot(NA, xlim=c(0, 40), ylim=range(df$fit), xaxs="i", yaxs="i", xlab = "Days", ylab = "log(S:H)", main=title)
  for (i in 1:length(dfsplit)) {
    with(dfsplit[[i]], {
      lines(Days, fit)
      addpoly(Days, lci, uci, col=alpha(colors[i], 0.5))
    })
  }
  
  legend("bottom", legend=levels(df$Result), pch=15, bty="n",
         col=alpha(colors, 0.4), pt.cex=1.5,
         y.intersp=0.3)
}


#Normalize qPCR value
for(res in data_pcrRes) {
  sdata <- NA
  sdata <- data_pcrRes[which(data_pcrRes$Result==res),]
  for(day in sdata$Days) {
    sdata[which(sdata$Days ==day), "rfit"] <- sdata[which(sdata$Days==day), "fit"] / sdata[which(sdata$Days==0), "fit"]
    sdata[which(sdata$Days ==day), "rlci"] <- sdata[which(sdata$Days==day), "lci"] / sdata[which(sdata$Days==0), "fit"]
    sdata[which(sdata$Days ==day), "ruci"] <- sdata[which(sdata$Days==day), "uci"] / sdata[which(sdata$Days==0), "fit"]
  }
  data_pcrRes[which(data_pcrRes$Result==res),"rfit"] <- sdata[, "rfit"]
  data_pcrRes[which(data_pcrRes$Result==res), "rlci"] <- sdata[, "rlci"]
  data_pcrRes[which(data_pcrRes$Result==res), "ruci"] <- sdata[, "ruci"]
}

plotfitspcrresnormal(data_pcrRes, c("above", "below", "not"), NA, F)

#qPCR to rfvfm at day 32... should normalize
dayThT_qPCR <- qPCR_Results[which(qPCR_Results$Days==33),]

summary_genopcr <- ddply(qPCR_Results, c("Genotype", "Days"), summarise, N    = length(logSH),
                          mean = mean(logSH),
                          sd   = sd(logSH),
                          se   = sd / sqrt(N)
)


combineDayThT<- merge(dayThT, dayThT_qPCR, by="Genotype")
plot(logSH~ rfvfm, data=combineDayThT)
pcrrfvfm <- lm(logSH~ rfvfm, data=combineDayThT)
abline(pcrrfvfm)
anova(pcrrfvfm)


mod_pcrdayThT <-lm(logSH ~ Result, data=qPCR_Results)
anova(mod_pcrdayThT)
cld(lsmeans(mod_pcrdayThT, "Result"))

#qpcr to symbiont density... should make relative?


combineDayThT_cell <- merge(dayThT_qPCR, combineCells, by=c("Genotype", "TimePoint"))
plot(totalcellcm ~logSH, data=combineDayThT_cell)
symqpcr <- lm(totalcellcm~logSH,data=combineDayThT_cell)
anova(symqpcr)
pcrcell <- merge(qPCR_Results, combineCells, by=c("Genotype", "TimePoint"))
plot(totalcellcm ~ logSH, data=pcrcell)
xyplot(totalcellcm ~ logSH, data=pcrcell, groups=Result.x)
symqpcr <- lm(totalcellcm~logSH, data=pcrcell)
anova(symqpcr)
abline(symqpcr)

#Test day differences
mod_pcrResulttime <- lm(logSH ~factor(Days) * Result, data=qPCR_Results )
anova(mod_pcrResulttime)
TukeyHSD(aov(mod_pcrResulttime))
View(cld(lsmeans(mod_pcrResulttime, c("Result", "Days"))))

#Relative S:H~ relative cells plot
#Normalize pcr geno fits
data_pcrGeno

for(geno in data_pcrGeno$Genotype) {
  sdata <- NA
  sdata <- data_pcrGeno[which(data_pcrGeno$Genotype==geno),]
  for(day in sdata$Days) {
    sdata[which(sdata$Days ==day), "rlog"] <- sdata[which(sdata$Days==day), "fit"] / sdata[which(sdata$Days==0), "fit"]
    sdata[which(sdata$Days ==day), "rlci"] <- sdata[which(sdata$Days==day), "lci"] / sdata[which(sdata$Days==0), "fit"]
    sdata[which(sdata$Days ==day), "ruci"] <- sdata[which(sdata$Days==day), "uci"] / sdata[which(sdata$Days==0), "fit"]
  }
  data_pcrGeno[which(data_pcrGeno$Genotype==geno),"rlog"] <- sdata[, "rlog"]
  data_pcrGeno[which(data_pcrGeno$Genotype==geno), "rlci"] <- sdata[, "rlci"]
  data_pcrGeno[which(data_pcrGeno$Genotype==geno), "ruci"] <- sdata[, "ruci"]
}

for(geno in data_symGeno$Genotype) {
  sdata <- NA
  sdata <- data_symGeno[which(data_symGeno$Genotype==geno),]
  for(day in sdata$Days) {
    sdata[which(sdata$Days ==day), "rsym"] <- sdata[which(sdata$Days==day), "fit"] / sdata[which(sdata$Days==0), "fit"]
    sdata[which(sdata$Days ==day), "rlci"] <- sdata[which(sdata$Days==day), "lci"] / sdata[which(sdata$Days==0), "fit"]
    sdata[which(sdata$Days ==day), "ruci"] <- sdata[which(sdata$Days==day), "uci"] / sdata[which(sdata$Days==0), "fit"]
  }
  data_symGeno[which(data_symGeno$Genotype==geno),"rsym"] <- sdata[, "rsym"]
  data_symGeno[which(data_symGeno$Genotype==geno), "rlci"] <- sdata[, "rlci"]
  data_symGeno[which(data_symGeno$Genotype==geno), "ruci"] <- sdata[, "ruci"]
}
relcellsh <- merge(data_symGeno, data_pcrGeno, by=c("Genotype", "Days"))
plot(rsym~rlog, relcellsh[which(relcellsh$Days==32),], type="p")
rfvfmcellsh <- merge(relcellsh, newdata, by=c("Genotype", "Days"))
plot(fit.x~rfvfm, rfvfmcellsh[which(rfvfmcellsh$Days==32),], type="p")
mod_symvcf <- lm(fit.x~rfvfm, rfvfmcellsh[which(rfvfmcellsh$Days==32),])
anova(mod_symvcf)
summary(mod_symvcf)
