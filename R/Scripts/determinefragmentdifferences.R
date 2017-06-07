#Write Genotype Table for geno 3
genotypes3$Letter <- substr(genotypes3$Genotype.ID, start = 5, 5)
geno3_letters <-aggregate(genotypes3$Letter, by = list(genotypes3$Genotype), FUN = c)
geno3_count <- aggregate(genotypes3$Letter, by = list(genotypes3$Genotype), FUN = length)
View(merge(geno3_count, geno3_letters, by = "Group.1"))
geno3table<-merge(geno3_count, geno3_letters, by = "Group.1")

#Write Genotype Table for geno 2
genotypes2 <- subset(genotypes2, !genotypes2$Genotype %in% c("17BL", "BACK", "Blan"))
genotypes2$Letter <- substr(genotypes2$Genotype.ID, start = 5, 5)
geno2_letters <- aggregate(as.character(genotypes2$Letter), by = list(genotypes2$Genotype), FUN = paste)
geno2_count <- aggregate(genotypes2$Letter, by = list(genotypes2$Genotype), FUN = length)
geno2table<-merge(geno2_count, geno2_letters, by = "Group.1")

#Change Column Names for genotype table
library(reshape2)
head(geno2table)
colnames(geno2table) <- gsub("Group.1", "Genotype", colnames(geno2table))
colnames(geno2table) <- gsub("Frag Count", "FragCount", colnames(geno2table))
colnames(geno2table) <- gsub("Frag ID", "FragID", colnames(geno2table))
geno2table$value <- as.integer(2) #differentiates timpe points

colnames(geno3table) <- gsub("Group.1", "Genotype", colnames(geno3table))
colnames(geno3table) <- gsub("Frag Count", "FragCount", colnames(geno3table))
colnames(geno3table) <- gsub("Frag ID", "FragID", colnames(geno3table))
geno3table$value <- as.integer(3)

#determine which fragments differ
geno2table <- geno2table[order(geno2table$Genotype), ]
geno3table <- geno3table[order(geno3table$Genotype), ]
symdiff <- function(x, y) setdiff(union(x, y), intersect(x, y))

result <- data.frame(matrix(NA, nrow=nrow(geno2table), ncol=2))
for (i in 1:nrow(geno2table)) {
  frags2 <- geno2table[i, "FragID"][[1]]
  frags3 <- geno3table[i, "FragID"][[1]]
  if (setequal(frags2, frags3)) {
    result[i, 1] <- TRUE
    result[i, 2] <- NA
  } else {
    result[i, 1] <- FALSE
    result[i, 2] <- paste0(symdiff(frags2, frags3), collapse=" ")
  }
}

View(cbind(geno2table$Genotype, result))





combined$FragID <- vapply(combined$FragID, paste, collapse = ", ", character(1L))
View(combined)
write.table(combined, file = "/Users/Rich/Documents/GradSchool/Internship/Lab/Data/AcroporaStressTest/Output/Table/combinedgenotable.csv")