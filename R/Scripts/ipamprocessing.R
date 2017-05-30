#Import Genotype CSV
options(stringsAsFactors = FALSE)
genotypes1 <- read.csv("/Users/Rich/Documents/Grad School/Internship/Lab/Data/AcroporaStressTest/Data/22052017_initialIpam/Processed/R_genotypes.csv")
genotypes2 <- read.csv("/Users/Rich/Documents/Grad School/Internship/Lab/Data/AcroporaStressTest/Data/25052017_2ndIpam/Processed/rGenotypes.csv")
genotypes3 <- read.csv("/Users/Rich/Documents/Grad School/Internship/Lab/Data/AcroporaStressTest/Data/29052017_3rdIpam/Processed/3rdgenotypes.csv")

#Import iPAM data
devtools::source_url("https://raw.githubusercontent.com/jrcunning/IPAM2R/master/R/import_ipam.R")
ipam1_values <- import_ipam("/Users/Rich/Documents/Grad School/Internship/Lab/Data/AcroporaStressTest/Data/22052017_initialIpam/Raw Data/iPAMinitials", 
                            info.pattern = NULL)
ipam2_values <- import_ipam("/Users/Rich/Documents/Grad School/Internship/Lab/Data/AcroporaStressTest/Data/25052017_2ndIpam/RawData/25052017_2ndIpam", 
                            info.pattern = NULL)
ipam3_values <- import_ipam("/Users/Rich/Documents/Grad School/Internship/Lab/Data/AcroporaStressTest/Data/29052017_3rdIpam/Raw Data", 
                            info.pattern = NULL)
combine_ipam1 <- merge(genotypes1, ipam1_values, by = c("file", "AOI"))
combine_ipam1$Date <- as.Date("2017-05-22")
combine_ipam1$Tank <- NA
combine_ipam2 <- merge(genotypes2, ipam2_values, by = c("file", "AOI"))
combine_ipam2$Date <- as.Date("2017-05-25") # make new column with date field
combine_ipam3 <- merge(genotypes3, ipam3_values, by = c("file", "AOI"))
combine_ipam3$Date <- as.Date("2017-05-29")

# combine time points
# reorder columns in ipam1
combine_ipam1 <- combine_ipam1[,c(1:5,10,6:9)]
# attach all rows together
all <- rbind(combine_ipam1, combine_ipam2, combine_ipam3)
# remove rows with "BL" or "BACK" in them
all.f <- subset(all, !all$Genotype %in% c("17BL", "BACK", "Blan"))

library(lattice)
xyplot(Y ~ Date | Genotype, data=all.f, type=c("p", "r"), ylim=c(0.4, 0.7))
summ <- aggregate(all.f$Frag.ID, by = list(all.f$Genotype, all.f$Date), FUN = length)
xyplot(x~ Group.2|Group.1, data = summ, type ="o", ylim = c(0, 17))
xyplot(Y ~ Tank | Genotype, data=all.f, type=c("p", "r"), ylim=c(0.4, 0.7))

# NEW LINE OF CODE!!!!!!!!!

#find a value
all[which(all$Genotype=="17BL"),]
all[which(all$Y > 1.0),]


#basic analysis
head(combine_ipam3)
boxplot(combine_ipam3$Y ~ combine_ipam3$Genotype, ylim=c(0.3,0.7))
mod <- lm(combine_ipam3$Y ~ combine_ipam3$Genotype)
anova(mod)
TukeyHSD(aov(mod))
table(combine_ipam2$Genotype)
View(ipam3_values)

