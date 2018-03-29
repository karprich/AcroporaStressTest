pr<- read.csv("/Users/Rich/Desktop/prac.csv", header=T)
barplot(pr$average)
arrows(pr$average, pr$std+pr$average, pr$average, pr$average-pr$std, angle = 90, length = .1)
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

barx <- barplot(pr$average, names.arg = pr$type,ylim=c(0,10), col="blue", axis.lty=1, ylab="Units")
error.bar(barx, pr$average, 1.96*pr$std)
