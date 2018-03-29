#Import Genotype CSV
options(stringsAsFactors = FALSE)
genotypes1 <- read.csv("Data/22052017_initialIpam/Processed/R_genotypes.csv")
genotypes2 <- read.csv("Data/25052017_2ndIpam/Processed/rGenotypes.csv")
genotypes3 <- read.csv("Data/29052017_3rdIpam/Processed/3rdgenotypes.csv")
genotypes4 <- read.csv("Data/01062017_4thIpam/Processed/4rgenotypes.csv")
genotypes5 <- read.csv("Data/05062017_5thIpam/Processed/rGenotypes.csv")
genotypes5wbox <- read.csv("Data/05062017_5thIpam/Processed/rGenotypeswboxtype.csv")
genotypes6 <- read.csv("Data/08062017_6thIpam/Processed/rgenotypes.csv")
genotypes6wbox <- read.csv("Data/08062017_6thIpam/Processed/rgenotypeswbox.csv")
genotypes7 <- read.csv("Data/12062017_7thIpam/Processed/rgenotypes.csv")
genotypes8 <- read.csv("Data/15062017_8thIpam/Processed/rgenotypes.csv")
genotypes9 <- read.csv("Data/19062017_9thIpam/Processed/rgenotypes.csv")
genotypes10 <- read.csv("Data/21062017_10thIpam/Processed/rgenotypes.csv")
genotypes11 <- read.csv("Data/23062017_11thIpam/Processed/rgenotypes.csv")
genotypes12 <- read.csv("Data/26062017_12thIpam/Processed/rgenotypes.csv")
genotypes13 <- read.csv("Data/28062017_13thIpam/Processed/rgenotypes.csv")
genotypes14 <- read.csv("Data/30062017_14thIpam/Processed/rgenotypes.csv")
genotypes15 <- read.csv("Data/AcroporaRecovery/13072017_InitialIpam/Processed/rgenotypes.csv")
genotypes16 <- read.csv("Data/AcroporaRecovery/20072017_2ndIpam/Processed/rgenotypes.csv")
genotypes17 <- read.csv("Data/AcroporaRecovery/15082017_3rdIpam/Processed/rgenotypes.csv")
genotypes18_own <- read.csv("Data/AcroporaRecovery/23082017_4thIpam/OwnHead/Processed/rgenotypes.csv")
genotypes18_loan <- read.csv("Data/AcroporaRecovery/23082017_4thIpam/LoanerHead/Processed/rgenotypes.csv")

#Import iPAM data
devtools::source_url("https://raw.githubusercontent.com/jrcunning/IPAM2R/master/R/import_ipam.R")
ipam1_values <- import_ipam("Data/22052017_initialIpam/Raw Data/iPAMinitials", 
                            info.pattern = NULL)
ipam2_values <- import_ipam("Data/25052017_2ndIpam/RawData/25052017_2ndIpam", 
                            info.pattern = NULL)
ipam3_values <- import_ipam("Data/29052017_3rdIpam/Raw Data", 
                            info.pattern = NULL)
ipam4_values <- import_ipam("Data/01062017_4thIpam/RawData", 
                            info.pattern = NULL)
ipam5_values <- import_ipam("Data/05062017_5thIpam/RawData", 
                            info.pattern = NULL)
ipam6_values <- import_ipam("Data/08062017_6thIpam/RawData", 
                            info.pattern = NULL)
ipam7_values <- import_ipam("Data/12062017_7thIpam/RawData", 
                            info.pattern = NULL)
ipam8_values <- import_ipam("Data/15062017_8thIpam/RawData", 
                            info.pattern = NULL)
ipam9_values <- import_ipam("Data/19062017_9thIpam/RawData",  info.pattern = NULL)  
ipam10_values <- import_ipam("Data/21062017_10thIpam/RawData", info.pattern = NULL)
ipam11_values <- import_ipam("Data/23062017_11thIpam/RawData", info.pattern = NULL)
ipam12_values <- import_ipam("Data/26062017_12thIpam/RawData", info.pattern = NULL)
ipam13_values <- import_ipam("Data/28062017_13thIpam/RawData", info.pattern = NULL)
ipam14_values <- import_ipam("Data/30062017_14thIpam/RawData", info.pattern = NULL)
ipam15_values <- import_ipam("Data/AcroporaRecovery/13072017_InitialIpam/RawData", info.pattern = NULL)
ipam16_values <- import_ipam("Data/AcroporaRecovery/20072017_2ndIpam/RawData", info.pattern = NULL)
ipam17_values <- import_ipam("Data/AcroporaRecovery/15082017_3rdIpam/RawData", info.pattern = NULL)
ipam18_own_values <- import_ipam("Data/AcroporaRecovery/23082017_4thIpam/OwnHead/RawData", info.pattern = NULL)
ipam18_loan_values <- import_ipam("Data/AcroporaRecovery/23082017_4thIpam/LoanerHead/RawData", info.pattern = NULL)


# Check Photos and Files for Same AOI
View(table(ipam15_values$file))
View(table(genotypes15$Picture))



combine_ipam1 <- merge(genotypes1, ipam1_values, by = c("file", "AOI"))
combine_ipam1$Date <- as.Date("2017-05-22")
combine_ipam1$Tank <- "Initial"
combine_ipam2 <- merge(genotypes2, ipam2_values, by = c("file", "AOI"))
combine_ipam2$Date <- as.Date("2017-05-25") # make new column with date field
combine_ipam3 <- merge(genotypes3, ipam3_values, by = c("file", "AOI"))
combine_ipam3$Date <- as.Date("2017-05-29")
combine_ipam4 <- merge(genotypes4, ipam4_values, by = c("file", "AOI"))
combine_ipam4$Date <- as.Date("2017-06-01")
combine_ipam5 <- merge(genotypes5, ipam5_values, by = c("file", "AOI"))
combine_ipam5$Date <- as.Date("2017-06-05")
combine_ipam5wbox <- merge(genotypes5wbox, ipam5_values, by = c("file", "AOI"))
combine_ipam5wbox$Date <- as.Date("2017-06-05")
combine_ipam6 <- merge(genotypes6, ipam6_values, by = c("file", "AOI"))
combine_ipam6$Date <- as.Date("2017-06-08")
combine_ipam6wbox <- merge(genotypes6wbox, ipam6_values, by = c("file", "AOI"))
combine_ipam6wbox$Date <- as.Date("2017-06-08")
combine_ipam7 <- merge(genotypes7, ipam7_values, by = c("file", "AOI"))
combine_ipam7$Date <- as.Date("2017-06-12")
combine_ipam8 <- merge(genotypes8, ipam8_values, by = c("file", "AOI"))
combine_ipam8$Date <- as.Date("2017-06-15")
combine_ipam9 <- merge(genotypes9, ipam9_values, by = c("file", "AOI"))
combine_ipam9$Date <- as.Date("2017-06-19")
combine_ipam10 <- merge(genotypes10, ipam10_values, by = c("file", "AOI"))
combine_ipam10$Date <- as.Date("2017-06-21")
combine_ipam11 <- merge(genotypes11, ipam11_values, by = c("file", "AOI"))
combine_ipam11$Date <- as.Date("2017-06-23")
combine_ipam12 <- merge(genotypes12, ipam12_values, by = c("file", "AOI"))
combine_ipam12$Date <- as.Date("2017-06-26")
combine_ipam13 <- merge(genotypes13, ipam13_values, by = c("file", "AOI"))
combine_ipam13$Date <- as.Date("2017-06-28")
combine_ipam14 <- merge(genotypes14, ipam14_values, by = c("file", "AOI"))
combine_ipam14$Date <- as.Date("2017-06-30")
combine_ipam15 <- merge(genotypes15, ipam15_values, by = c("file", "AOI"))
combine_ipam15$Date <- as.Date("2017-07-13")
combine_ipam16 <- merge(genotypes16, ipam16_values, by = c("file", "AOI"))
combine_ipam16$Date <- as.Date("2017-07-20")
combine_ipam17 <- merge(genotypes17, ipam17_values, by = c("file", "AOI"))
combine_ipam17$Date <- as.Date("2017-08-15")
combine_ipam18_loan <- merge(genotypes18_loan, ipam18_loan_values, by = c("file", "AOI"))
combine_ipam18_loan$Date <- as.Date("2017-08-22")
combine_ipam18_own <- merge(genotypes18_own, ipam18_own_values, by = c("file", "AOI"))
combine_ipam18_own$Date <- as.Date("2017-08-22")


#correction factor
cor_combine_ipam5 <- combine_ipam5
cor_combine_ipam4 <- combine_ipam4
#cor_combine_ipam4$F <- 0.999 * combine_ipam4$F
cor_combine_ipam4$Fm <- 1.03 * combine_ipam4$Fm
cor_combine_ipam4$Y <- (cor_combine_ipam4$Fm-cor_combine_ipam4$F)/cor_combine_ipam4$Fm

#cor_combine_ipam5$F <- 0.999 * combine_ipam5$F
cor_combine_ipam5$Fm <- 1.03 * combine_ipam5$Fm
cor_combine_ipam5$Y <- (cor_combine_ipam5$Fm-cor_combine_ipam5$F)/cor_combine_ipam5$Fm
#par(mfrow= c(1, 2))
#boxplot(combine_ipam1$Y, combine_ipam2$Y, combine_ipam3$Y, cor_combine_ipam4$Y, cor_combine_ipam5$Y, combine_ipam6$Y, names = c("1st IPam", "2nd IPam", "3rd IPam", "4th IPam", "5th IPam", "6th IPam"), ylab = "Y", ylim = c(0.3, 0.7))
#boxplot(combine_ipam1$Y, combine_ipam2$Y, combine_ipam3$Y, combine_ipam4$Y, combine_ipam5$Y, combine_ipam6$Y, names = c("1st IPam", "2nd IPam", "3rd IPam", "4th IPam", "5th IPam", "6th IPam"), ylab = "Y", ylim = c(0.3, 0.7))

# combine time points
# reorder columns in ipam1
combine_ipam1 <- combine_ipam1[,c(1:5,10,6:9)]
# attach all rows together
allRec <- rbind(combine_ipam1, combine_ipam2, combine_ipam3, cor_combine_ipam4, cor_combine_ipam5, combine_ipam6, combine_ipam7, combine_ipam8, combine_ipam9, combine_ipam10, combine_ipam11, combine_ipam12, combine_ipam13, combine_ipam14, combine_ipam15, combine_ipam16, combine_ipam17, combine_ipam18_loan)
# remove rows with "BL" or "BACK" in them
all.rec <- subset(allRec, !allRec$Genotype %in% c("17BL", "BACK", "Blan"))
# remove rows with 1732 and 1734 from 12th Ipam as they are dead
all.rec[which(all.f$Date=="2017-06-26" & all.f$Genotype=="1732"),"Y"]<- NA
all.rec[which(all.f$Date=="2017-06-26" & all.f$Genotype=="1734"),"Y"] <- NA
#omit Data that were positively selected for:
#1721 Ipam 13+14
all.rec[which(all.f$Date=="2017-06-28" || all.f$Date=="2017-06-30" & all.f$Genotype=="1721"),"Y"]<- NA
#1722 Ipam 14
all.rec[which(all.f$Date=="2017-06-30" & all.f$Genotype=="1722"),"Y"]<- NA
#1727 Ipam 14
all.rec[which(all.f$Date=="2017-06-30" & all.f$Genotype=="1727"),"Y"]<- NA
#1728 Ipam 14
all.rec[which(all.f$Date=="2017-06-30" & all.f$Genotype=="1728"),"Y"]<- NA
#1731 Ipam 13+14
all.rec[which(all.f$Date=="2017-06-28" || all.f$Date=="2017-06-30" & all.f$Genotype=="1731"),"Y"]<- NA
#1732 Ipam 12, 13+14 Previously deleted Ipam 12
all.[which(all.f$Date=="2017-06-26" || all.f$Date=="2017-06-28" || all.f$Date=="2017-06-30" & all.f$Genotype=="1732"),"Y"]<- NA
#1734 Ipam 12, 13+14 Previously deleted Ipam 12
all.f[which(all.f$Date=="2017-06-26" || all.f$Date=="2017-06-28" || all.f$Date=="2017-06-30" & all.f$Genotype=="1734"),"Y"]<- NA
#1736 Ipam 14
all.f[which(all.f$Date=="2017-06-30" & all.f$Genotype=="1736"),"Y"]<- NA
#1737 Ipam 14
all.f[which(all.f$Date=="2017-06-30" & all.f$Genotype=="1737"),"Y"]<- NA
#1738 Ipam 13+14
all.f[which(all.f$Date=="2017-06-28" || all.f$Date=="2017-06-30" & all.f$Genotype=="1738"),"Y"]<- NA
#1739 Ipam 13+14
all.f[which(all.f$Date=="2017-06-28" || all.f$Date=="2017-06-30" & all.f$Genotype=="1739"),"Y"]<- NA
#1747 Ipam 14
all.f[which(all.f$Date=="2017-06-30" & all.f$Genotype=="1747"),"Y"]<- NA
#1750 Ipam 14
all.f[which(all.f$Date=="2017-06-30" & all.f$Genotype=="1750"),"Y"]<- NA
#1752 Ipam 14
all.f[which(all.f$Date=="2017-06-30" & all.f$Genotype=="1752"),"Y"]<- NA
#1753 Ipam 14
all.f[which(all.f$Date=="2017-06-30" & all.f$Genotype=="1753"),"Y"]<- NA
#1755 Ipam 13+14
all.f[which(all.f$Date=="2017-06-28" || all.f$Date=="2017-06-30" & all.f$Genotype=="1755"),"Y"]<- NA
#1757 Ipam 14
all.f[which(all.f$Date=="2017-06-30" & all.f$Genotype=="1757"),"Y"]<- NA
#1758 Ipam 14
all.f[which(all.f$Date=="2017-06-30" & all.f$Genotype=="1758"),"Y"]<- NA
#1759 Ipam 14
all.f[which(all.f$Date=="2017-06-30" & all.f$Genotype=="1759"),"Y"]<- NA

#Plot

library(lattice)
png(filename = "Output/Figures/2ndrecovery.png", width=7, height = 7, units = "in", res = 300)
xyplot(Y ~ Date | Genotype, data=all.rec, type=c("p"), ylim=c(0.0, 0.7))
dev.off()
png(filename = "Output/Figures/1strecoveryboxplot.png", width=7, height = 7, units = "in", res = 300)
boxplot(combine_ipam1$Y, combine_ipam2$Y, combine_ipam3$Y, cor_combine_ipam4$Y, cor_combine_ipam5$Y, combine_ipam6$Y, combine_ipam7$Y, combine_ipam8$Y, combine_ipam9$Y, combine_ipam10$Y, combine_ipam11$Y, combine_ipam12$Y, combine_ipam13$Y, combine_ipam14$Y, combine_ipam15$Y, names = c("1st IPam", "2nd IPam", "3rd IPam", "4th IPam", "5th IPam", "6th IPam", "7th IPam", "8th IPam", "9th IPam", "10th IPam", "11th IPam", "12th IPam", "13th Ipam", "14th Ipam", "15th Ipam"), ylab = "Y", ylim = c(0.0, 0.7))
boxplot(Y~Date, data=all.rec, names = c("1st IPam", "2nd IPam", "3rd IPam", "4th IPam", "5th IPam", "6th IPam", "7th IPam", "8th IPam", "9th IPam", "10th IPam", "11th IPam", "12th IPam", "13th Ipam", "14th Ipam", "15th Ipam", "16th Ipam", "17th Ipam", "18th Ipam"), ylab = "Y", ylim = c(0.0, 0.7))
dev.off()

xyplot(Y ~ factor(Tank) | factor(Genotype), data=combine_ipam15, type=c("p"), ylim = c(0.0, 0.7), xlab = "Tank", main = "Ipam 9 by Tank")

#Compare Ipam arrays for Ipam 18
own <- combine_ipam18_own
own$array <- "own"
loan <- combine_ipam18_loan
loan$array <- "loan"
array <- rbind(own, loan)
boxplot(Y~array, data=array)
mod_array <- lm(Y~array, data=array)
anova(mod_array)
TukeyHSD(aov(mod_array))
xyplot(Y~factor(array)|factor(FragID), data = array, type=c("p"))
t.test(own$Y, loan$Y)
