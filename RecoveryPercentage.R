# Determine Percent of Genos recovered
options(stringsAsFactors = FALSE)
recovery <- read.csv("Data/AcroporaRecovery/rdata.csv", header = T)
recovery$Genotype.Initial.Count <- as.character(recovery$Genotype.Initial.Count)
recovery$GenoFirstIpam <- as.character(recovery$GenoFirstIpam)
recovery$GenosRemaining <- as.character(recovery$GenosRemaining)



quickCount <- data.frame(table(recovery$Genotype.Initial.Count))
initialPam <- data.frame(table(recovery$GenoFirstIpam))
remain <- data.frame(table(recovery$GenosRemaining))

quickCount$Var1 <- as.character(quickCount$Var1)
initialPam$Var1 <- as.character(initialPam$Var1)
remain$Var1 <- as.character(remain$Var1)
complete <- data.frame(NA)
for(i in 1:nrow(quickCount)) {
  complete[i, "Genotype"] <- quickCount[i, "Var1"]
  complete[i, "FreqInit"] <- quickCount[i, "Freq"]
}
complete<-complete[,c(2:3)]
for(i in 1:nrow(initialPam)) {
  for(j in 1:nrow(complete)){
    if(initialPam[i, "Var1"]==complete[j, "Genotype"]) {
      complete[j, "FreqPam"] <- initialPam[i, "Freq"]
    } 
  }
}
for(i in 1:nrow(complete)) {
  if(is.na(complete[i, "FreqPam"])) {
    complete[i, "FreqPam"] <- 0
  }
}


for(i in 1:nrow(remain)) {
  for(j in 1:nrow(complete)){
    if(remain[i, "Var1"]==complete[j, "Genotype"]) {
      complete[j, "FreqRemain"] <- remain[i, "Freq"]
    } 
  }
}
for(i in 1:nrow(complete)) {
  if(is.na(complete[i, "FreqRemain"])) {
    complete[i, "FreqRemain"] <- 0
  }
}

complete$PercentRecovery <- complete$FreqRemain/complete$FreqPam 
complete$PercentMortality <- 1-complete$PercentRecovery
hist(complete$PercentMortality)
sum(complete$FreqRemain)/sum(complete$FreqPam)
sum(complete$FreqRemain)/sum(complete$FreqInit)



