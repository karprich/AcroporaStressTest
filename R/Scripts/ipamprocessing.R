#Import Genotype CSV
options(stringsAsFactors = FALSE)
genotypes1 <- read.csv("/Users/Rich/Documents/GradSchool/Internship/Lab/Data/AcroporaStressTest/Data/22052017_initialIpam/Processed/R_genotypes.csv")
genotypes2 <- read.csv("/Users/Rich/Documents/GradSchool/Internship/Lab/Data/AcroporaStressTest/Data/25052017_2ndIpam/Processed/rGenotypes.csv")
genotypes3 <- read.csv("/Users/Rich/Documents/GradSchool/Internship/Lab/Data/AcroporaStressTest/Data/29052017_3rdIpam/Processed/3rdgenotypes.csv")
genotypes4 <- read.csv("/Users/Rich/Documents/GradSchool/Internship/Lab/Data/AcroporaStressTest/Data/01062017_4thIpam/Processed/4rgenotypes.csv")
genotypes5 <- read.csv("/Users/Rich/Documents/GradSchool/Internship/Lab/Data/AcroporaStressTest/Data/05062017_5thIpam/Processed/rGenotypes.csv")
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

# combine time points
# reorder columns in ipam1
combine_ipam1 <- combine_ipam1[,c(1:5,10,6:9)]
# attach all rows together
all <- rbind(combine_ipam1, combine_ipam2, combine_ipam3, combine_ipam4, combine_ipam5)
# remove rows with "BL" or "BACK" in them
all.f <- subset(all, !all$Genotype %in% c("17BL", "BACK", "Blan"))

library(lattice)
xyplot(Y ~ Date | Genotype, data=all.f, type=c("p", "r"), ylim=c(0.3, 0.7))
boxplot(combine_ipam3$Y, combine_ipam4$Y, combine_ipam5$Y, names = c("3rd IPam", "4th IPam", "5th IPam"), ylab = "Y", ylim = c(0.3, 0.7))
summ <- aggregate(all.f$Frag.ID, by = list(all.f$Genotype, all.f$Date), FUN = length)
xyplot(x~ Group.2|Group.1, data = summ, type ="o", ylim = c(0, 17))
xyplot(Y ~ factor(Tank) | factor(Genotype), data=combine_ipam5, type=c("p"), ylim = c(0.4, 0.7))
xyplot(Y ~ factor(Tank) | factor(Genotype), data=combine_ipam4, type=c("p"), ylim = c(0.4, 0.7))
View(table(combine_ipam3$Tank))
? xyplot 

View(combine_ipam3)
#Write figures
png(filename = "Output/Figures/first3ipams.png", width=5, height = 5, units = "in", res = 300)
xyplot(Y ~ Date | Genotype, data=all.f, type=c("p", "r"), ylim=c(0.4, 0.7))
dev.off()

#Number of Genotypes with x number of genotypes
length(which(geno2_count$x >= 8))


#Number of geno types per tank
View(table(combine_ipam3$Genotype, combine_ipam3$Tank))
ipam3_tanks <- dcast(combine_ipam3, combine_ipam3$Genotype ~ combine_ipam3$Tank)
View(ipam3_tanks)
#Write figures
png(filename = "Output/Figures/first3ipams.png", width=5, height = 5, units = "in", res = 300)
xyplot(Y ~ Date | Genotype, data=all.f, type=c("p", "r"), ylim=c(0.4, 0.7))
dev.off()
# NEW LINE OF CODE!!!!!!!!!

#find a value
all[which(all$Genotype=="17BL"),]
all[which(all$Y > 1.0),]


#basic analysis
head(combine_ipam3)
boxplot(combine_ipam3$Y ~ combine_ipam3$Genotype, ylim=c(0.3,0.7))
mod <- lm(combine_ipam5$Y ~ combine_ipam5$Genotype)
anova(mod)
TukeyHSD(aov(mod))
table(combine_ipam2$Genotype)
View(ipam3_values)

