#Write Genotype Table for geno 7
genotypes7$Letter <- substr(genotypes7$FragID, start = 5, 5)
geno7_letters <-aggregate(genotypes7$Letter, by = list(genotypes7$Genotype), FUN = paste)
geno7_count <- aggregate(genotypes7$Letter, by = list(genotypes7$Genotype), FUN = length)
geno7table<-merge(geno7_count, geno7_letters, by = "Group.1")

#Write Genotype Table for geno 8
#genotypes2 <- subset(genotypes2, !genotypes2$Genotype %in% c("17BL", "BACK", "Blan"))
genotypes8$Letter <- substr(genotypes8$FragID, start = 5, 5)
geno8_letters <- aggregate(as.character(genotypes8$Letter), by = list(genotypes8$Genotype), FUN = paste)
geno8_count <- aggregate(genotypes8$Letter, by = list(genotypes8$Genotype), FUN = length)
geno8table<-merge(geno8_count, geno8_letters, by = "Group.1")


genotypes9$Letter <- substr(genotypes9$FragID, start = 5, 5)
geno9_letters <- aggregate(as.character(genotypes9$Letter), by = list(genotypes9$Genotype), FUN = paste)
geno9_count <- aggregate(genotypes9$Letter, by = list(genotypes9$Genotype), FUN = length)
geno9table<-merge(geno9_count, geno9_letters, by = "Group.1")

#Number of Genotypes with x number of genotypes
length(which(geno6_count$x >= 8))

View(table(combine_ipam6$Genotype))


#Change Column Names for genotype table
library(reshape2)

colnames(geno7table) <- gsub("Group.1", "Genotype", colnames(geno7table))
colnames(geno7table) <- gsub("x.x", "FragCount", colnames(geno7table))
colnames(geno7table) <- gsub("x.y", "FragID", colnames(geno7table))
geno7table$value <- as.integer(7) #differentiates timpe points

colnames(geno8table) <- gsub("Group.1", "Genotype", colnames(geno8table))
colnames(geno8table) <- gsub("x.x", "FragCount", colnames(geno8table))
colnames(geno8table) <- gsub("x.y", "FragID", colnames(geno8table))
geno8table$value <- as.integer(8)

colnames(geno9table) <- gsub("Group.1", "Genotype", colnames(geno9table))
colnames(geno9table) <- gsub("x.x", "FragCount", colnames(geno9table))
colnames(geno9table) <- gsub("x.y", "FragID", colnames(geno9table))
geno9table$value <- as.integer(9)

#determine which fragments differ from 8 to 9
geno9table <- geno9table[order(geno9table$Genotype), ]
geno8table <- geno8table[order(geno8table$Genotype), ]
symdiff <- function(x, y) setdiff(union(x, y), intersect(x, y))

result <- data.frame(matrix(NA, nrow=nrow(geno8table), ncol=2))
for (i in 1:nrow(geno8table)) {
  frags9 <- geno9table[i, "FragID"][[1]]
  frags8 <- geno8table[i, "FragID"][[1]]
  if (setequal(frags9, frags8)) {
    result[i, 1] <- TRUE
    result[i, 2] <- NA
  } else {
    result[i, 1] <- FALSE
    result[i, 2] <- paste0(symdiff(frags9, frags8), collapse=" ")
  }
}

View(cbind(geno7table$Genotype, result))


#Determine which fragments differ from 8 to 7
geno7table <- geno7table[order(geno7table$Genotype), ]
geno8table <- geno8table[order(geno8table$Genotype), ]
symdiff <- function(x, y) setdiff(union(x, y), intersect(x, y))

result <- data.frame(matrix(NA, nrow=nrow(geno8table), ncol=2))
for (i in 1:nrow(geno7table)) {
  frags7 <- geno7table[i, "FragID"][[1]]
  frags8 <- geno8table[i, "FragID"][[1]]
  if (setequal(frags7, frags8)) {
    result[i, 1] <- TRUE
    result[i, 2] <- NA
  } else {
    result[i, 1] <- FALSE
    result[i, 2] <- paste0(symdiff(frags7, frags8), collapse=" ")
  }
}

View(cbind(geno7table$Genotype, result))

#number of genotypes Table
number_geno_7 <- table(genotypes7$Genotype)
number_geno_8 <- table(genotypes8$Genotype)
number_geno_9 <- table(genotypes9$Genotype)

#View by numbers of fragments per
View(cbind(number_geno_7, number_geno_8))
View(cbind(number_geno_9, number_geno_8))


combined$FragID <- vapply(combined$FragID, paste, collapse = ", ", character(1L))
View(combined)
write.table(combined, file = "/Users/Rich/Documents/GradSchool/Internship/Lab/Data/AcroporaStressTest/Output/Table/combinedgenotable.csv")