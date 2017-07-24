library(reshape2)

#Write Genotype Table for geno 1
genotypes1 <- subset(genotypes1, !genotypes1$Genotype %in% c("17BL", "BACK", "Blan"))
genotypes1$Letter <- substr(genotypes1$FragID, start = 5, 5)
geno1_letters <-aggregate(genotypes1$Letter, by = list(genotypes1$Genotype), FUN = paste)
geno1_count <- aggregate(genotypes1$Letter, by = list(genotypes1$Genotype), FUN = length)
geno1table<-merge(geno1_count, geno1_letters, by = "Group.1")
colnames(geno1table) <- gsub("Group.1", "Genotype", colnames(geno1table))
colnames(geno1table) <- gsub("x.x", "FragCount", colnames(geno1table))
colnames(geno1table) <- gsub("x.y", "FragID", colnames(geno1table))
geno1table$value <- as.integer(1) #differentiates timpe points



#Write Genotype Table for geno 2
genotypes2 <- subset(genotypes2, !genotypes2$Genotype %in% c("17BL", "BACK", "Blan"))
genotypes2$Letter <- substr(genotypes2$FragID, start = 5, 5)
geno2_letters <- aggregate(as.character(genotypes2$Letter), by = list(genotypes2$Genotype), FUN = paste)
geno2_count <- aggregate(genotypes2$Letter, by = list(genotypes2$Genotype), FUN = length)
geno2table<-merge(geno2_count, geno2_letters, by = "Group.1")
colnames(geno2table) <- gsub("Group.1", "Genotype", colnames(geno2table))
colnames(geno2table) <- gsub("x.x", "FragCount", colnames(geno2table))
colnames(geno2table) <- gsub("x.y", "FragID", colnames(geno2table))
geno2table$value <- as.integer(2) #differentiates timpe points


#Write Genotype Table for geno 3
genotypes3$Letter <- substr(genotypes3$FragID, start = 5, 5)
geno3_letters <-aggregate(genotypes3$Letter, by = list(genotypes3$Genotype), FUN = paste)
geno3_count <- aggregate(genotypes3$Letter, by = list(genotypes3$Genotype), FUN = length)
geno3table<-merge(geno3_count, geno3_letters, by = "Group.1")
colnames(geno3table) <- gsub("Group.1", "Genotype", colnames(geno3table))
colnames(geno3table) <- gsub("x.x", "FragCount", colnames(geno3table))
colnames(geno3table) <- gsub("x.y", "FragID", colnames(geno3table))
geno3table$value <- as.integer(3) #differentiates timpe points

#Write Genotype Table for geno 4
genotypes4$Letter <- substr(genotypes4$FragID, start = 5, 5)
geno4_letters <-aggregate(genotypes4$Letter, by = list(genotypes4$Genotype), FUN = paste)
geno4_count <- aggregate(genotypes4$Letter, by = list(genotypes4$Genotype), FUN = length)
geno4table<-merge(geno4_count, geno4_letters, by = "Group.1")
colnames(geno4table) <- gsub("Group.1", "Genotype", colnames(geno4table))
colnames(geno4table) <- gsub("x.x", "FragCount", colnames(geno4table))
colnames(geno4table) <- gsub("x.y", "FragID", colnames(geno4table))
geno4table$value <- as.integer(4) #differentiates timpe points


#Write Genotype Table for geno 5
genotypes5$Letter <- substr(genotypes5$FragID, start = 5, 5)
geno5_letters <-aggregate(genotypes5$Letter, by = list(genotypes5$Genotype), FUN = paste)
geno5_count <- aggregate(genotypes5$Letter, by = list(genotypes5$Genotype), FUN = length)
geno5table<-merge(geno5_count, geno5_letters, by = "Group.1")
colnames(geno5table) <- gsub("Group.1", "Genotype", colnames(geno5table))
colnames(geno5table) <- gsub("x.x", "FragCount", colnames(geno5table))
colnames(geno5table) <- gsub("x.y", "FragID", colnames(geno5table))
geno5table$value <- as.integer(5) #differentiates timpe points

#Write Genotype Table for geno 6
genotypes6$Letter <- substr(genotypes6$FragID, start = 5, 5)
geno6_letters <-aggregate(genotypes6$Letter, by = list(genotypes6$Genotype), FUN = paste)
geno6_count <- aggregate(genotypes6$Letter, by = list(genotypes6$Genotype), FUN = length)
geno6table<-merge(geno6_count, geno6_letters, by = "Group.1")
colnames(geno6table) <- gsub("Group.1", "Genotype", colnames(geno6table))
colnames(geno6table) <- gsub("x.x", "FragCount", colnames(geno6table))
colnames(geno6table) <- gsub("x.y", "FragID", colnames(geno6table))
geno6table$value <- as.integer(6) #differentiates timpe points

#Write Genotype Table for geno 7
genotypes7$Letter <- substr(genotypes7$FragID, start = 5, 5)
geno7_letters <-aggregate(genotypes7$Letter, by = list(genotypes7$Genotype), FUN = paste)
geno7_count <- aggregate(genotypes7$Letter, by = list(genotypes7$Genotype), FUN = length)
geno7table<-merge(geno7_count, geno7_letters, by = "Group.1")
colnames(geno7table) <- gsub("Group.1", "Genotype", colnames(geno7table))
colnames(geno7table) <- gsub("x.x", "FragCount", colnames(geno7table))
colnames(geno7table) <- gsub("x.y", "FragID", colnames(geno7table))
geno7table$value <- as.integer(7) #differentiates timpe points

#Write Genotype Table for geno 8
genotypes8$Letter <- substr(genotypes8$FragID, start = 5, 5)
geno8_letters <-aggregate(genotypes8$Letter, by = list(genotypes8$Genotype), FUN = paste)
geno8_count <- aggregate(genotypes8$Letter, by = list(genotypes8$Genotype), FUN = length)
geno8table<-merge(geno8_count, geno8_letters, by = "Group.1")
colnames(geno8table) <- gsub("Group.1", "Genotype", colnames(geno8table))
colnames(geno8table) <- gsub("x.x", "FragCount", colnames(geno8table))
colnames(geno8table) <- gsub("x.y", "FragID", colnames(geno8table))
geno8table$value <- as.integer(8) #differentiates timpe points

#Write Genotype Table for geno 9
genotypes9$Letter <- substr(genotypes9$FragID, start = 5, 5)
geno9_letters <-aggregate(genotypes9$Letter, by = list(genotypes9$Genotype), FUN = paste)
geno9_count <- aggregate(genotypes9$Letter, by = list(genotypes9$Genotype), FUN = length)
geno9table<-merge(geno9_count, geno9_letters, by = "Group.1")
colnames(geno9table) <- gsub("Group.1", "Genotype", colnames(geno9table))
colnames(geno9table) <- gsub("x.x", "FragCount", colnames(geno9table))
colnames(geno9table) <- gsub("x.y", "FragID", colnames(geno9table))
geno9table$value <- as.integer(9) #differentiates timpe points

#Write Genotype Table for geno 10
genotypes10$Letter <- substr(genotypes10$FragID, start = 5, 5)
geno10_letters <-aggregate(genotypes10$Letter, by = list(genotypes10$Genotype), FUN = paste)
geno10_count <- aggregate(genotypes10$Letter, by = list(genotypes10$Genotype), FUN = length)
geno10table<-merge(geno10_count, geno10_letters, by = "Group.1")
colnames(geno10table) <- gsub("Group.1", "Genotype", colnames(geno10table))
colnames(geno10table) <- gsub("x.x", "FragCount", colnames(geno10table))
colnames(geno10table) <- gsub("x.y", "FragID", colnames(geno10table))
geno10table$value <- as.integer(10) #differentiates timpe points

#Write Genotype Table for geno 11
genotypes11$Letter <- substr(genotypes11$FragID, start = 5, 5)
geno11_letters <-aggregate(genotypes11$Letter, by = list(genotypes11$Genotype), FUN = paste)
geno11_count <- aggregate(genotypes11$Letter, by = list(genotypes11$Genotype), FUN = length)
geno11table<-merge(geno11_count, geno11_letters, by = "Group.1")
colnames(geno11table) <- gsub("Group.1", "Genotype", colnames(geno11table))
colnames(geno11table) <- gsub("x.x", "FragCount", colnames(geno11table))
colnames(geno11table) <- gsub("x.y", "FragID", colnames(geno11table))
geno11table$value <- as.integer(11) #differentiates timpe points

#Write Genotype Table for geno 12
genotypes12$Letter <- substr(genotypes12$FragID, start = 5, 5)
geno12_letters <-aggregate(genotypes12$Letter, by = list(genotypes12$Genotype), FUN = paste)
geno12_count <- aggregate(genotypes12$Letter, by = list(genotypes12$Genotype), FUN = length)
geno12table<-merge(geno12_count, geno12_letters, by = "Group.1")
colnames(geno12table) <- gsub("Group.1", "Genotype", colnames(geno12table))
colnames(geno12table) <- gsub("x.x", "FragCount", colnames(geno12table))
colnames(geno12table) <- gsub("x.y", "FragID", colnames(geno12table))
geno12table$value <- as.integer(12) #differentiates timpe points

#Write Genotype Table for geno 13
genotypes13$Letter <- substr(genotypes13$FragID, start = 5, 5)
geno13_letters <-aggregate(genotypes13$Letter, by = list(genotypes13$Genotype), FUN = paste)
geno13_count <- aggregate(genotypes13$Letter, by = list(genotypes13$Genotype), FUN = length)
geno13table<-merge(geno13_count, geno13_letters, by = "Group.1")
colnames(geno13table) <- gsub("Group.1", "Genotype", colnames(geno13table))
colnames(geno13table) <- gsub("x.x", "FragCount", colnames(geno13table))
colnames(geno13table) <- gsub("x.y", "FragID", colnames(geno13table))
geno13table$value <- as.integer(13) #differentiates timpe points

#Write Genotype Table for geno 14
genotypes14$Letter <- substr(genotypes14$FragID, start = 5, 5)
geno14_letters <-aggregate(genotypes14$Letter, by = list(genotypes14$Genotype), FUN = paste)
geno14_count <- aggregate(genotypes14$Letter, by = list(genotypes14$Genotype), FUN = length)
geno14table<-merge(geno14_count, geno14_letters, by = "Group.1")
colnames(geno14table) <- gsub("Group.1", "Genotype", colnames(geno14table))
colnames(geno14table) <- gsub("x.x", "FragCount", colnames(geno14table))
colnames(geno14table) <- gsub("x.y", "FragID", colnames(geno14table))
geno14table$value <- as.integer(14) #differentiates timpe points

#Number of Genotypes with x number of genotypes
#length(which(geno6_count$x >= 8))
#View(table(combine_ipam6$Genotype))


#Change Column Names for genotype table
#library(reshape2)
#colnames(geno7table) <- gsub("Group.1", "Genotype", colnames(geno7table))
#colnames(geno7table) <- gsub("x.x", "FragCount", colnames(geno7table))
#colnames(geno7table) <- gsub("x.y", "FragID", colnames(geno7table))
#geno7table$value <- as.integer(7) #differentiates timpe points

#determine which fragments differ from 1 to 2
geno1table <- geno1table[order(geno1table$Genotype), ]
geno2table <- geno2table[order(geno2table$Genotype), ]
symdiff <- function(x, y) setdiff(union(x, y), intersect(x, y))

result <- data.frame(matrix(NA, nrow=nrow(geno1table), ncol=2))
for (i in 1:nrow(geno1table)) {
  frags1 <- geno1table[i, "FragID"][[1]]
  frags2 <- geno2table[i, "FragID"][[1]]
  if (setequal(frags1, frags2)) {
    result[i, 1] <- TRUE
    result[i, 2] <- NA
  } else {
    result[i, 1] <- FALSE
    result[i, 2] <- paste0(symdiff(frags1, frags2), collapse=" ")
  }
}
View(cbind(geno1table$Genotype, result))
number_geno_1 <- table(genotypes1$Genotype)
number_geno_2 <- table(genotypes2$Genotype)
#View by numbers of fragments per
View(cbind(number_geno_1, number_geno_2))


#determine which fragments differ from 2 to 3
geno3table <- geno3table[order(geno3table$Genotype), ]
geno2table <- geno2table[order(geno2table$Genotype), ]
symdiff <- function(x, y) setdiff(union(x, y), intersect(x, y))

result <- data.frame(matrix(NA, nrow=nrow(geno2table), ncol=2))
for (i in 1:nrow(geno2table)) {
  frags3 <- geno3table[i, "FragID"][[1]]
  frags2 <- geno2table[i, "FragID"][[1]]
  if (setequal(frags3, frags2)) {
    result[i, 1] <- TRUE
    result[i, 2] <- NA
  } else {
    result[i, 1] <- FALSE
    result[i, 2] <- paste0(symdiff(frags3, frags2), collapse=" ")
  }
}
View(cbind(geno2table$Genotype, result))
number_geno_3 <- table(genotypes3$Genotype)
number_geno_2 <- table(genotypes2$Genotype)
#View by numbers of fragments per
View(cbind(number_geno_3, number_geno_2))
View(which(duplicated(genotypes3$FragID)))

#determine which fragments differ from 3 to 4
geno3table <- geno3table[order(geno3table$Genotype), ]
geno4table <- geno4table[order(geno4table$Genotype), ]
symdiff <- function(x, y) setdiff(union(x, y), intersect(x, y))

result <- data.frame(matrix(NA, nrow=nrow(geno3table), ncol=2))
for (i in 1:nrow(geno3table)) {
  frags3 <- geno3table[i, "FragID"][[1]]
  frags4 <- geno4table[i, "FragID"][[1]]
  if (setequal(frags3, frags4)) {
    result[i, 1] <- TRUE
    result[i, 2] <- NA
  } else {
    result[i, 1] <- FALSE
    result[i, 2] <- paste0(symdiff(frags3, frags4), collapse=" ")
  }
}
View(cbind(geno3table$Genotype, result))
#number of genotypes Table
number_geno_3 <- table(genotypes3$Genotype)
number_geno_4 <- table(genotypes4$Genotype)
#View by numbers of fragments per
View(cbind(number_geno_3, number_geno_4))
View(which(duplicated(genotypes3$FragID)))

#determine which fragments differ from 4 to 5
geno5table <- geno5table[order(geno5table$Genotype), ]
geno4table <- geno4table[order(geno4table$Genotype), ]
symdiff <- function(x, y) setdiff(union(x, y), intersect(x, y))

result <- data.frame(matrix(NA, nrow=nrow(geno4table), ncol=2))
for (i in 1:nrow(geno4table)) {
  frags5 <- geno5table[i, "FragID"][[1]]
  frags4 <- geno4table[i, "FragID"][[1]]
  if (setequal(frags5, frags4)) {
    result[i, 1] <- TRUE
    result[i, 2] <- NA
  } else {
    result[i, 1] <- FALSE
    result[i, 2] <- paste0(symdiff(frags5, frags4), collapse=" ")
  }
}
View(cbind(geno4table$Genotype, result))

#determine which fragments differ from 5 to 6
geno5table <- geno5table[order(geno5table$Genotype), ]
geno6table <- geno6table[order(geno6table$Genotype), ]
symdiff <- function(x, y) setdiff(union(x, y), intersect(x, y))

result <- data.frame(matrix(NA, nrow=nrow(geno5table), ncol=2))
for (i in 1:nrow(geno5table)) {
  frags5 <- geno5table[i, "FragID"][[1]]
  frags6 <- geno6table[i, "FragID"][[1]]
  if (setequal(frags5, frags6)) {
    result[i, 1] <- TRUE
    result[i, 2] <- NA
  } else {
    result[i, 1] <- FALSE
    result[i, 2] <- paste0(symdiff(frags5, frags6), collapse=" ")
  }
}
View(cbind(geno5table$Genotype, result))

#determine which fragments differ from 6 to 7
geno7table <- geno7table[order(geno7table$Genotype), ]
geno6table <- geno6table[order(geno6table$Genotype), ]
symdiff <- function(x, y) setdiff(union(x, y), intersect(x, y))

result <- data.frame(matrix(NA, nrow=nrow(geno6table), ncol=2))
for (i in 1:nrow(geno6table)) {
  frags7 <- geno7table[i, "FragID"][[1]]
  frags6 <- geno6table[i, "FragID"][[1]]
  if (setequal(frags7, frags6)) {
    result[i, 1] <- TRUE
    result[i, 2] <- NA
  } else {
    result[i, 1] <- FALSE
    result[i, 2] <- paste0(symdiff(frags7, frags6), collapse=" ")
  }
}
View(cbind(geno6table$Genotype, result))

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

View(cbind(geno8table$Genotype, result))

#determine which fragments differ from 9 to 10
geno9table <- geno9table[order(geno9table$Genotype), ]
geno10table <- geno10table[order(geno10table$Genotype), ]
symdiff <- function(x, y) setdiff(union(x, y), intersect(x, y))

result <- data.frame(matrix(NA, nrow=nrow(geno9table), ncol=2))
for (i in 1:nrow(geno9table)) {
  frags9 <- geno9table[i, "FragID"][[1]]
  frags10 <- geno10table[i, "FragID"][[1]]
  if (setequal(frags9, frags10)) {
    result[i, 1] <- TRUE
    result[i, 2] <- NA
  } else {
    result[i, 1] <- FALSE
    result[i, 2] <- paste0(symdiff(frags9, frags10), collapse=" ")
  }
}

View(cbind(geno9table$Genotype, result))


number_geno_9 <- table(genotypes9$Genotype)
number_geno_10 <- table(genotypes10$Genotype)
#View by numbers of fragments per
View(cbind(number_geno_9, number_geno_10))
genotypes10[which(duplicated(genotypes10$FragID)),]


#determine which fragments differ from 10 to 11
geno11table <- geno11table[order(geno11table$Genotype), ]
geno10table <- geno10table[order(geno10table$Genotype), ]
symdiff <- function(x, y) setdiff(union(x, y), intersect(x, y))

result <- data.frame(matrix(NA, nrow=nrow(geno10table), ncol=2))
for (i in 1:nrow(geno10table)) {
  frags11 <- geno11table[i, "FragID"][[1]]
  frags10 <- geno10table[i, "FragID"][[1]]
  if (setequal(frags11, frags10)) {
    result[i, 1] <- TRUE
    result[i, 2] <- NA
  } else {
    result[i, 1] <- FALSE
    result[i, 2] <- paste0(symdiff(frags11, frags10), collapse=" ")
  }
}

View(cbind(geno10table$Genotype, result))
number_geno_11 <- table(genotypes11$Genotype)
number_geno_10 <- table(genotypes10$Genotype)
#View by numbers of fragments per
View(cbind(number_geno_10, number_geno_11))


#determine which fragments differ from 11 to 12
geno11table <- geno11table[order(geno11table$Genotype), ]
geno12table <- geno12table[order(geno12table$Genotype), ]
symdiff <- function(x, y) setdiff(union(x, y), intersect(x, y))

result <- data.frame(matrix(NA, nrow=nrow(geno11table), ncol=2))
for (i in 1:nrow(geno11table)) {
  frags11 <- geno11table[i, "FragID"][[1]]
  frags12 <- geno12table[i, "FragID"][[1]]
  if (setequal(frags11, frags12)) {
    result[i, 1] <- TRUE
    result[i, 2] <- NA
  } else {
    result[i, 1] <- FALSE
    result[i, 2] <- paste0(symdiff(frags11, frags12), collapse=" ")
  }
}

View(cbind(geno11table$Genotype, result))

number_geno_11 <- table(genotypes11$Genotype)
number_geno_12 <- table(genotypes12$Genotype)
#View by numbers of fragments per
View(cbind(number_geno_11, number_geno_12))


#determine which fragments differ from 12 to 13
geno13table <- geno13table[order(geno13table$Genotype), ]
geno12table <- geno12table[order(geno12table$Genotype), ]
symdiff <- function(x, y) setdiff(union(x, y), intersect(x, y))

result <- data.frame(matrix(NA, nrow=nrow(geno12table), ncol=2))
for (i in 1:nrow(geno12table)) {
  for(j in 1: nrow(geno13table)) {
  if( geno12table[i, "Genotype"]==geno13table[j, "Genotype"]) {
  frags12 <- geno12table[i, "FragID"][[1]]
  frags13 <- geno13table[j, "FragID"][[1]]
  if (setequal(frags13, frags12)) {
    result[i, 1] <- TRUE
    result[i, 2] <- NA
  } else {
    result[i, 1] <- FALSE
    result[i, 2] <- paste0(symdiff(frags13, frags12), collapse=" ")
  }
  }
  }
}

View(cbind(geno12table$Genotype, result))

#determine which fragments differ from 13 to 14
geno13table <- geno13table[order(geno13table$Genotype), ]
geno14table <- geno14table[order(geno14table$Genotype), ]
symdiff <- function(x, y) setdiff(union(x, y), intersect(x, y))

result <- data.frame(matrix(NA, nrow=nrow(geno13table), ncol=2))
for (i in 1:nrow(geno13table)) {
  for(j in 1: nrow(geno14table)) {
    if(geno13table[i, "Genotype"]==geno14table[j, "Genotype"]) {
      frags14 <- geno14table[j, "FragID"][[1]]
      frags13 <- geno13table[i, "FragID"][[1]]
      if (setequal(frags13, frags14)) {
        result[i, 1] <- TRUE
        result[i, 2] <- NA
      } else {
        result[i, 1] <- FALSE
        result[i, 2] <- paste0(symdiff(frags13, frags14), collapse=" ")
      }
    }
  }
}

View(cbind(geno13table$Genotype, result))


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