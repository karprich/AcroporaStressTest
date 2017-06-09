#Write Genotype Table for geno 3
genotypes4$Letter <- substr(genotypes4$FragID, start = 5, 5)
geno4_letters <-aggregate(genotypes4$Letter, by = list(genotypes4$Genotype), FUN = paste)
geno4_count <- aggregate(genotypes4$Letter, by = list(genotypes4$Genotype), FUN = length)
View(merge(geno4_count, geno4_letters, by = "Group.1"))
geno4table<-merge(geno4_count, geno4_letters, by = "Group.1")

#Write Genotype Table for geno 2
genotypes2 <- subset(genotypes2, !genotypes2$Genotype %in% c("17BL", "BACK", "Blan"))
genotypes5$Letter <- substr(genotypes5$FragID, start = 5, 5)
geno5_letters <- aggregate(as.character(genotypes5$Letter), by = list(genotypes5$Genotype), FUN = paste)
geno5_count <- aggregate(genotypes5$Letter, by = list(genotypes5$Genotype), FUN = length)
geno5table<-merge(geno5_count, geno5_letters, by = "Group.1")


genotypes6$Letter <- substr(genotypes6$FragID, start = 5, 5)
geno6_letters <- aggregate(as.character(genotypes6$Letter), by = list(genotypes6$Genotype), FUN = paste)
geno6_count <- aggregate(genotypes6$Letter, by = list(genotypes6$Genotype), FUN = length)
geno5table<-merge(geno5_count, geno5_letters, by = "Group.1")

#Number of Genotypes with x number of genotypes
length(which(geno6_count$x >= 8))

View(table(combine_ipam6$Genotype))


#Change Column Names for genotype table
library(reshape2)
head(geno2table)
colnames(geno4table) <- gsub("Group.1", "Genotype", colnames(geno4table))
colnames(geno4table) <- gsub("Frag Count", "FragCount", colnames(geno4table))
colnames(geno4table) <- gsub("Frag ID", "FragID", colnames(geno4table))
geno4table$value <- as.integer(4) #differentiates timpe points

colnames(geno5table) <- gsub("Group.1", "Genotype", colnames(geno5table))
colnames(geno5table) <- gsub("Frag Count", "FragCount", colnames(geno5table))
colnames(geno5table) <- gsub("Frag ID", "FragID", colnames(geno5table))
geno5table$value <- as.integer(5)

#determine which fragments differ
geno5table <- geno5table[order(geno5table$Genotype), ]
geno4table <- geno4table[order(geno4table$Genotype), ]
symdiff <- function(x, y) setdiff(union(x, y), intersect(x, y))

result <- data.frame(matrix(NA, nrow=nrow(geno4table), ncol=2))
for (i in 1:nrow(geno2table)) {
  frags4 <- geno4table[i, "FragID"][[1]]
  frags5 <- geno5table[i, "FragID"][[1]]
  if (setequal(frags5, frags4)) {
    result[i, 1] <- TRUE
    result[i, 2] <- NA
  } else {
    result[i, 1] <- FALSE
    result[i, 2] <- paste0(symdiff(frags4, frags5), collapse=" ")
  }
}

View(cbind(geno4table$Genotype, result))





combined$FragID <- vapply(combined$FragID, paste, collapse = ", ", character(1L))
View(combined)
write.table(combined, file = "/Users/Rich/Documents/GradSchool/Internship/Lab/Data/AcroporaStressTest/Output/Table/combinedgenotable.csv")