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
genotypes7 <- read.csv("Data/120517_7thIpam/Processed/rgenotypes.csv")
genotypes8 <- read.csv("Data/15062017_8thIpam/Processed/rgenotypes.csv")
genotypes9 <- read.csv("Data/19062017_9thIpam/Processed/rgenotypes.csv")
genotypes10 <- read.csv("Data/21062017_10thIpam/Processed/rgenotypes.csv")
genotypes11 <- read.csv("Data/23062017_11thIpam/Processed/rgenotypes.csv")
genotypes12 <- read.csv("Data/26062017_12thIpam/Processed/rgenotypes.csv")
genotypes13 <- read.csv("Data/28062017_13thIpam/Processed/rgenotypes.csv")
genotypes14 <- read.csv("Data/30062017_14thIpam/Processed/rgenotypes.csv")

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
ipam7_values <- import_ipam("Data/120517_7thIpam/RawData", 
                            info.pattern = NULL)
ipam8_values <- import_ipam("Data/15062017_8thIpam/RawData", 
                            info.pattern = NULL)
ipam9_values <- import_ipam("Data/19062017_9thIpam/RawData",  info.pattern = NULL)  
ipam10_values <- import_ipam("Data/21062017_10thIpam/RawData", info.pattern = NULL)
ipam11_values <- import_ipam("Data/23062017_11thIpam/RawData", info.pattern = NULL)
ipam12_values <- import_ipam("Data/26062017_12thIpam/RawData", info.pattern = NULL)
ipam13_values <- import_ipam("Data/28062017_13thIpam/RawData", info.pattern = NULL)
ipam14_values <- import_ipam("Data/30062017_14thIpam/RawData", info.pattern = NULL)


# Check Photos and Files for Same AOI
View(table(ipam14_values$file))
View(table(genotypes14$Picture))



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
all <- rbind(combine_ipam1, combine_ipam2, combine_ipam3, cor_combine_ipam4, cor_combine_ipam5, combine_ipam6, combine_ipam7, combine_ipam8, combine_ipam9, combine_ipam10, combine_ipam11, combine_ipam12, combine_ipam13, combine_ipam14)
# remove rows with "BL" or "BACK" in them
all.f <- subset(all, !all$Genotype %in% c("17BL", "BACK", "Blan"))
# remove rows with 1732 and 1734 from 12th Ipam as they are dead
all.f[which(all.f$Date=="2017-06-26" & all.f$Genotype=="1732"),"Y"]<- NA
all.f[which(all.f$Date=="2017-06-26" & all.f$Genotype=="1734"),"Y"] <- NA
#omit Data that were positively selected for:
#1721 Ipam 13+14
all.f[which(all.f$Date=="2017-06-28" || all.f$Date=="2017-06-30" & all.f$Genotype=="1721"),"Y"]<- NA
#1722 Ipam 14
all.f[which(all.f$Date=="2017-06-30" & all.f$Genotype=="1722"),"Y"]<- NA
#1727 Ipam 14
all.f[which(all.f$Date=="2017-06-30" & all.f$Genotype=="1727"),"Y"]<- NA
#1728 Ipam 14
all.f[which(all.f$Date=="2017-06-30" & all.f$Genotype=="1728"),"Y"]<- NA
#1731 Ipam 13+14
all.f[which(all.f$Date=="2017-06-28" || all.f$Date=="2017-06-30" & all.f$Genotype=="1731"),"Y"]<- NA
#1732 Ipam 12, 13+14 Previously deleted Ipam 12
all.f[which(all.f$Date=="2017-06-26" || all.f$Date=="2017-06-28" || all.f$Date=="2017-06-30" & all.f$Genotype=="1732"),"Y"]<- NA
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


# Save data as .RData
save(all.f, file="Output/all.f.RData")


#correct all
#par(mfrow= c(1, 1))

library(lattice)
xyplot(Y ~ Date | Genotype, data=all.f, type=c("p"), ylim=c(0.0, 0.7))

boxplot(combine_ipam1$Y, combine_ipam2$Y, combine_ipam3$Y, cor_combine_ipam4$Y, cor_combine_ipam5$Y, combine_ipam6$Y, combine_ipam7$Y, combine_ipam8$Y, combine_ipam9$Y, combine_ipam10$Y, combine_ipam11$Y, combine_ipam12$Y, combine_ipam13$Y, combine_ipam14$Y, names = c("1st IPam", "2nd IPam", "3rd IPam", "4th IPam", "5th IPam", "6th IPam", "7th IPam", "8th IPam", "9th IPam", "10th IPam", "11th IPam", "12th IPam", "13th Ipam", "14th Ipam"), ylab = "Y", ylim = c(0.0, 0.7))
#To look at numbers by genotype and decline
#summ <- aggregate(all.f$FragID, by = list(all.f$Genotype, all.f$Date), FUN = length)
#xyplot(x~ Group.2|Group.1, data = summ, type ="o", ylim = c(0, 17))

xyplot(Y ~ factor(Tank) | factor(Genotype), data=combine_ipam14, type=c("p"), ylim = c(0.0, 0.7), xlab ="Tank")
#xyplot(Y ~ factor(Box) | factor(Genotype), data=combine_ipam6wbox, type=c("p"), ylim = c(0.3, 0.7))
View(table(combine_ipam3$Tank))
? xyplot 

View(combine_ipam3)
#Write figures
#Write Genotype Plot
png(filename = "Output/Figures/first14ipamsbygeno.png", width=7, height = 7, units = "in", res = 300)
xyplot(Y ~ Date | Genotype, data=all.f, type=c("p"), ylim=c(0.0, 0.7), main = "Ipam 14 by Genotype")
dev.off()
#Write Boxplots
png(filename = "Output/Figures/first14ipamsboxplot.png", width=3, height = 5, units = "in", res = 300)
boxplot(combine_ipam1$Y, combine_ipam2$Y, combine_ipam3$Y, cor_combine_ipam4$Y, cor_combine_ipam5$Y, combine_ipam6$Y, combine_ipam7$Y, combine_ipam8$Y, combine_ipam9$Y, combine_ipam10$Y, combine_ipam11$Y, combine_ipam12$Y, combine_ipam13$Y, combine_ipam14$Y, names = c("1st IPam", "2nd IPam", "3rd IPam", "4th IPam", "5th IPam", "6th IPam", "7th IPam", "8th IPam", "9th IPam", "10th IPam", "11th IPam", "12th IPam", "13th IPam", "14th IPam"), ylab = "Y", ylim = c(0.0, 0.7))
dev.off()
#Write Tanks plot
png(filename = "Output/Figures/ipam9bytank.png", width=5, height = 5, units = "in", res = 300)
xyplot(Y ~ factor(Tank) | factor(Genotype), data=combine_ipam14, type=c("p"), ylim = c(0.0, 0.7), xlab = "Tank", main = "Ipam 9 by Tank")
dev.off()

#mean Y Table
y_meangeno_10 <-aggregate(combine_ipam10$Y, by = list(combine_ipam10$Genotype), FUN = mean)
y_meangeno_11 <-aggregate(combine_ipam11$Y, by = list(combine_ipam11$Genotype), FUN = mean)
y_meangeno_12 <-aggregate(combine_ipam12$Y, by = list(combine_ipam12$Genotype), FUN = mean)
y_meangeno_13 <-aggregate(combine_ipam13$Y, by = list(combine_ipam13$Genotype), FUN = mean)
y_meangeno_14 <-aggregate(combine_ipam14$Y, by = list(combine_ipam14$Genotype), FUN = mean)
ipam_count <- table(combine_ipam12$Genotype)
ipam_count_13 <- table(combine_ipam13$Genotype)
ipam_count_14 <- table(combine_ipam14$Genotype)
View(aggregate(ipam_count_13, y_meangeno_13, FUN = paste))
View(aggregate(ipam_count_14, y_meangeno_14, FUN = paste))
View(cbind(ipam_count, y_meangeno_10, y_meangeno_11, y_meangeno_12, y_meangeno_13))
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
all.f[which(all.f$Genotype=="1732"),]
all.f[which(all.f$Date=="2017-06-26"),]
all.f[which(all.f$Date=="2017-06-26" & all.f$Genotype=="1732"),]

#basic analysis

#boxplot(combine_ipam6wbox$Y ~ combine_ipam6wbox$Box, ylim=c(0.3,0.7))
#mod_box <- lm(combine_ipam6wbox$Y ~ combine_ipam6wbox$Box)
library(multcomp)
library(multcompView)
library(lsmeans)


mod_Tank <- lm(Y ~ Tank, data= combine_ipam13)
combine_ipam13$Genotype <- factor(combine_ipam13$Genotype)
mod_Geno <- lm(Y ~ Genotype, data = combine_ipam13)

all.f$Date <- factor(all.f$Date)
all.f$Genotype <- factor(all.f$Genotype)
all.f$Tank <- factor(all.f$Tank)

mod_Date <-lm(Y ~ Date, data = all.f)
mod_Tank_all <-lm(Y ~ Tank, data= all.f)
mod_all <- lm(Y ~ Date * Genotype, data=all.f)
anova(mod_Geno)

cld(lsmeans(mod_Tank_all, "Tank"))
cld(lsmeans(mod_Tank, "Tank"))
cld(lsmeans(mod_Geno, "Genotype"))
cld(lsmeans(mod_Date, "Date"))
cld(lsmeans(mod_all, specs = c("Date", "Genotype")))


anova(mod_Tank)
TukeyHSD(aov(mod_Tank))
anova(mod_Date)
TukeyHSD(aov(mod_Date))
anova(mod_all)
TukeyHSD(aov(mod_all))
table(combine_ipam2$Genotype)
View(ipam3_values)

