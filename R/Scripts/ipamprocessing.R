#Import Genotype CSV
options(stringsAsFactors = FALSE)
genotypes1 <- read.csv("/Users/Rich/Documents/GradSchool/Internship/Lab/Data/AcroporaStressTest/Data/22052017_initialIpam/Processed/R_genotypes.csv")
genotypes2 <- read.csv("/Users/Rich/Documents/GradSchool/Internship/Lab/Data/AcroporaStressTest/Data/25052017_2ndIpam/Processed/rGenotypes.csv")
genotypes3 <- read.csv("/Users/Rich/Documents/GradSchool/Internship/Lab/Data/AcroporaStressTest/Data/29052017_3rdIpam/Processed/3rdgenotypes.csv")
genotypes4 <- read.csv("/Users/Rich/Documents/GradSchool/Internship/Lab/Data/AcroporaStressTest/Data/01062017_4thIpam/Processed/4rgenotypes.csv")
genotypes5 <- read.csv("/Users/Rich/Documents/GradSchool/Internship/Lab/Data/AcroporaStressTest/Data/05062017_5thIpam/Processed/rGenotypes.csv")
genotypes5wbox <- read.csv("/Users/Rich/Documents/GradSchool/Internship/Lab/Data/AcroporaStressTest/Data/05062017_5thIpam/Processed/rGenotypeswboxtype.csv")
genotypes6 <- read.csv("/Users/Rich/Documents/GradSchool/Internship/Lab/Data/AcroporaStressTest/Data/08062017_6thIpam/Processed/rgenotypes.csv")
genotypes6wbox <- read.csv("/Users/Rich/Documents/GradSchool/Internship/Lab/Data/AcroporaStressTest/Data/08062017_6thIpam/Processed/rgenotypeswbox.csv")
genotypes7 <- read.csv("/Users/Rich/Documents/GradSchool/Internship/Lab/Data/AcroporaStressTest/Data/120517_7thIpam/Processed/rgenotypes.csv")
genotypes8 <- read.csv("/Users/Rich/Documents/GradSchool/Internship/Lab/Data/AcroporaStressTest/Data/15062017_8thIpam/Processed/rgenotypes.csv")
genotypes9 <- read.csv("/Users/Rich/Documents/GradSchool/Internship/Lab/Data/AcroporaStressTest/Data/19062017_9thIpam/Processed/rgenotypes.csv")

#Import iPAM data
devtools::source_url("https://raw.githubusercontent.com/jrcunning/IPAM2R/master/R/import_ipam.R")
ipam1_values <- import_ipam("/Users/Rich/Documents/GradSchool/Internship/Lab/Data/AcroporaStressTest/Data/22052017_initialIpam/Raw Data/iPAMinitials", 
                            info.pattern = NULL)
ipam2_values <- import_ipam("/Users/Rich/Documents/GradSchool/Internship/Lab/Data/AcroporaStressTest/Data/25052017_2ndIpam/RawData/25052017_2ndIpam", 
                            info.pattern = NULL)
ipam3_values <- import_ipam("/Users/Rich/Documents/GradSchool/Internship/Lab/Data/AcroporaStressTest/Data/29052017_3rdIpam/Raw Data", 
                            info.pattern = NULL)
ipam4_values <- import_ipam("/Users/Rich/Documents/GradSchool/Internship/Lab/Data/AcroporaStressTest/Data/01062017_4thIpam/RawData", 
                            info.pattern = NULL)
ipam5_values <- import_ipam("/Users/Rich/Documents/GradSchool/Internship/Lab/Data/AcroporaStressTest/Data/05062017_5thIpam/RawData", 
                            info.pattern = NULL)
ipam6_values <- import_ipam("/Users/Rich/Documents/GradSchool/Internship/Lab/Data/AcroporaStressTest/Data/08062017_6thIpam/RawData", 
                            info.pattern = NULL)
ipam7_values <- import_ipam("/Users/Rich/Documents/GradSchool/Internship/Lab/Data/AcroporaStressTest/Data/120517_7thIpam/RawData", 
                            info.pattern = NULL)
ipam8_values <- import_ipam("/Users/Rich/Documents/GradSchool/Internship/Lab/Data/AcroporaStressTest/Data/15062017_8thIpam/RawData", 
                            info.pattern = NULL)
ipam9_values <- import_ipam("/Users/Rich/Documents/GradSchool/Internship/Lab/Data/AcroporaStressTest/Data/19062017_9thIpam/RawData",  info.pattern = NULL)                 

combine_ipam1 <- merge(genotypes1, ipam1_values, by = c("file", "AOI"))
combine_ipam1$Date <- as.Date("2017-05-22")
combine_ipam1$Tank <- NA
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

#correction factor
cor_combine_ipam5 <- combine_ipam5
cor_combine_ipam4 <- combine_ipam4
#cor_combine_ipam4$F <- 0.999 * combine_ipam4$F
cor_combine_ipam4$Fm <- 1.03 * combine_ipam4$Fm
cor_combine_ipam4$Y <- (cor_combine_ipam4$Fm-cor_combine_ipam4$F)/cor_combine_ipam4$Fm

#cor_combine_ipam5$F <- 0.999 * combine_ipam5$F
cor_combine_ipam5$Fm <- 1.03 * combine_ipam5$Fm
cor_combine_ipam5$Y <- (cor_combine_ipam5$Fm-cor_combine_ipam5$F)/cor_combine_ipam5$Fm
par(mfrow= c(1, 2))
boxplot(combine_ipam1$Y, combine_ipam2$Y, combine_ipam3$Y, cor_combine_ipam4$Y, cor_combine_ipam5$Y, combine_ipam6$Y, names = c("1st IPam", "2nd IPam", "3rd IPam", "4th IPam", "5th IPam", "6th IPam"), ylab = "Y", ylim = c(0.3, 0.7))
boxplot(combine_ipam1$Y, combine_ipam2$Y, combine_ipam3$Y, combine_ipam4$Y, combine_ipam5$Y, combine_ipam6$Y, names = c("1st IPam", "2nd IPam", "3rd IPam", "4th IPam", "5th IPam", "6th IPam"), ylab = "Y", ylim = c(0.3, 0.7))

# combine time points
# reorder columns in ipam1
combine_ipam1 <- combine_ipam1[,c(1:5,10,6:9)]
# attach all rows together
all <- rbind(combine_ipam1, combine_ipam2, combine_ipam3, cor_combine_ipam4, cor_combine_ipam5, combine_ipam6, combine_ipam7, combine_ipam8, combine_ipam9)
# remove rows with "BL" or "BACK" in them
all.f <- subset(all, !all$Genotype %in% c("17BL", "BACK", "Blan"))


#correct all
#par(mfrow= c(1, 1))

library(lattice)
xyplot(Y ~ Date | Genotype, data=all.f, type=c("p", "r"), ylim=c(0.3, 0.7))

boxplot(combine_ipam1$Y, combine_ipam2$Y, combine_ipam3$Y, cor_combine_ipam4$Y, cor_combine_ipam5$Y, combine_ipam6$Y, combine_ipam7$Y, combine_ipam8$Y, combine_ipam9$Y, names = c("1st IPam", "2nd IPam", "3rd IPam", "4th IPam", "5th IPam", "6th IPam", "7th IPam", "8th IPam", "9th IPam"), ylab = "Y", ylim = c(0.3, 0.7))
#To look at numbers by genotype and decline
#summ <- aggregate(all.f$FragID, by = list(all.f$Genotype, all.f$Date), FUN = length)
#xyplot(x~ Group.2|Group.1, data = summ, type ="o", ylim = c(0, 17))

xyplot(Y ~ factor(Tank) | factor(Genotype), data=combine_ipam9, type=c("p"), ylim = c(0.3, 0.7), xlab ="Tank")
xyplot(Y ~ factor(Box) | factor(Genotype), data=combine_ipam6wbox, type=c("p"), ylim = c(0.3, 0.7))
View(table(combine_ipam3$Tank))
? xyplot 

View(combine_ipam3)
#Write figures
#Write Genotype Plot
png(filename = "Output/Figures/first9ipamsbygeno.png", width=7, height = 7, units = "in", res = 300)
xyplot(Y ~ Date | Genotype, data=all.f, type=c("p", "r"), ylim=c(0.4, 0.7), main = "Ipam 9 by Genotype")
dev.off()
#Write Boxplots
png(filename = "Output/Figures/first9ipamsboxplot.png", width=3, height = 5, units = "in", res = 300)
boxplot(combine_ipam1$Y, combine_ipam2$Y, combine_ipam3$Y, cor_combine_ipam4$Y, cor_combine_ipam5$Y, combine_ipam6$Y, combine_ipam7$Y, combine_ipam8$Y, names = c("1st IPam", "2nd IPam", "3rd IPam", "4th IPam", "5th IPam", "6th IPam", "7th IPam", "8th IPam"), ylab = "Y", ylim = c(0.3, 0.7))
dev.off()
#Write Tanks plot
png(filename = "Output/Figures/ipam9bytank.png", width=5, height = 5, units = "in", res = 300)
xyplot(Y ~ factor(Tank) | factor(Genotype), data=combine_ipam8, type=c("p"), ylim = c(0.4, 0.7), xlab = "Tank", main = "Ipam 9 by Tank")
dev.off()

#mean Y Table
y_meangeno <-aggregate(combine_ipam9$Y, by = list(combine_ipam9$Genotype), FUN = mean)
write.csv(y_meangeno, "Output/8thIpamYmean.csv")

#Number of Genotypes with x number of genotypes
length(which(geno2_count$x >= 8))
View(which(combine_ipam5$Genotype < 8))
View(table(combine_ipam6$Genotype))

#Number of geno types per tank
View(table(combine_ipam5$Genotype, combine_ipam5$Tank))
View(table(combine_ipam5$Genotype))
library(reshape2)
ipam3_tanks <- dcast(combine_ipam3, combine_ipam3$Genotype ~ combine_ipam3$Tank)
ipam6_tanks <- dcast(combine_ipam6, combine_ipam6$Genotype ~ combine_ipam6$Tank)
View(ipam6_tanks)
#Write figures
png(filename = "Output/Figures/first3ipams.png", width=5, height = 5, units = "in", res = 300)
xyplot(Y ~ Date | Genotype, data=all.f, type=c("p", "r"), ylim=c(0.4, 0.7))
dev.off()
# NEW LINE OF CODE!!!!!!!!!

#find a value
all[which(all$Genotype=="17BL"),]
all[which(all$Y > 1.0),]
View(all.f[which(all.f$Y < 0.5),])

#basic analysis

#boxplot(combine_ipam6wbox$Y ~ combine_ipam6wbox$Box, ylim=c(0.3,0.7))
#mod_box <- lm(combine_ipam6wbox$Y ~ combine_ipam6wbox$Box)
boxplot(combine_ipam6$Y ~ combine_ipam6$Tank, ylim=c(0.3,0.7))
mod_Tank <- lm(combine_ipam8$Y ~ combine_ipam8$Tank)
mod_Geno <- lm(combine_ipam8$Y ~ factor(combine_ipam8$Genotype))
mod_Date <-lm(all.f$Y ~ factor(all.f$Date))
mod_all <- lm(all.f$Y ~ factor(all.f$Date) * factor(all.f$Genotype))
anova(mod_Geno)
TukeyHSD(aov(mod_Geno))

anova(mod_Tank)
TukeyHSD(aov(mod_Tank))
anova(mod_Date)
TukeyHSD(aov(mod_Date))
anova(mod_all)
TukeyHSD(aov(mod_all))
table(combine_ipam2$Genotype)
View(ipam3_values)

