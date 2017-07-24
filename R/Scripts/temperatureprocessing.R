#Temperature Analysis
# ? strptime
library(zoo)
temp <- read.csv("/Users/Rich/Documents/GradSchool/Internship/Lab/Data/AcroporaStressTest/Data/TemperatureData/Rdata1.csv", header = T)
temp$Date <- as.POSIXct(temp$Date, format = "%H:%M %m/%d/%Y")

#Add Threshhold Temperature of 32
temp$Thres <- 32

#Create Days based on GAMMs
temp$Days <- difftime(temp$Date, as.POSIXct("2017-05-22 21:00:00"), units="days")

#omit sensor
temp[temp$Date>as.POSIXct("2017-06-03 11:40:00") & temp$Date<as.POSIXct("2017-06-05 16:00:00"), "Tmpx17.bottom.Right"]<-NA
temp[temp$Date>as.POSIXct("2017-06-19 12:40:00") & temp$Date<as.POSIXct("2017-06-19 16:40:00"), "Tmpx19.top.left"]<-NA
# ADD Other cleaning times


#Average Temperature of three tanks
averageTemp<-temp[,c(1:4, 6, 8, 9)]
averageTemp$mean<-rowMeans(averageTemp[,c(2:4)], na.rm=T)
averageTemp$smooth <- rollapply(data = averageTemp$mean, width = 240, FUN = mean, na.rm = T, fill = NA)
plot(temp$Projected ~ temp$Days, type = "l", ylim=c(25, 35), xaxs="i", axes=F, lwd = 2, main = "TL", xlab = "Date", ylab = "Temperature (oC)", xlim=c(0,39))
lines(temp$Thres~temp$Days, col = 'orange', lty = 2)
lines(averageTemp$mean ~ averageTemp$Days, lwd = (.5))
lines(averageTemp$smooth ~ averageTemp$Days, col = 'red')


#Write Figure 
#png(filename = "Output/Figures/tempthrough190617.png", width=5, height = 5, units = "in", res = 300)
# for multipanel plot 
library(zoo)
par(mfrow = c(2, 2))
#Plot TL
smooth_tl <- rollapply(data = temp$Tmpx19..6b.top.left, width = 240, FUN = mean, na.rm = T, fill = NA)
plot(temp$Projected ~ temp$Date, type = "l", ylim=c(25, 35), lwd = 2, main = "TL", xlab = "Date", ylab = "Temperature (oC)")
lines(temp$Thres~temp$Date, col = 'orange', lty = 2)
lines(temp$Tmpx19..6b.top.left ~ temp$Date, lwd = (.5))
lines(smooth_tl ~ temp$Date, col = 'red')


#Plot TR
smooth_tr <- rollapply(data = temp$Tmpx21...7b.top.right, width = 240, FUN = mean, na.rm = T, fill = NA)
plot(temp$Projected ~ temp$Date, type = "l", ylim=c(25, 35), lwd = 2, main = "TR", xlab = "Date", ylab = "Temperature (oC)")
lines(temp$Thres~temp$Date, col = 'orange', lty = 2)
lines(temp$Tmpx21...7b.top.right ~ temp$Date, lwd = (.5))
lines(smooth_tr ~ temp$Date, col = 'red')
#Plot BL
smooth_bl <- rollapply(data = temp$Tmpx23.6A.Bottom.Left, width = 240, FUN = mean, na.rm = T, fill = NA, xlab = "Date", ylab = "Temperature (oC)")
plot(temp$ProjectedRoss ~ temp$Date, type = "l", ylim=c(25, 35), lwd = 2, main = "BL", xlab = "Date", ylab = "Temperature (oC)")
lines(temp$Thres~temp$Date, col = 'orange', lty = 2)
lines(temp$Tmpx23.6A.Bottom.Left ~ temp$Date, lwd = (.5))
lines(smooth_bl ~ temp$Date, col = 'red')

#Plot BR
smooth_br <- rollapply(data = temp$Tmpx17.bottom.Right, width = 240, FUN = mean, na.rm = T, fill = NA, xlab = "Date", ylab = "Temperature (oC)")
plot(temp$Projected ~ temp$Date, type = "l", ylim=c(25, 35), lwd = 2, main = "BR")
lines(temp$Thres~temp$Date, col = 'orange', lty = 2)
lines(temp$Tmpx17.bottom.Right ~ temp$Date, lwd = (.5))
lines(smooth_br ~ temp$Date, col = 'red')
#dev.off()

par(mfrow = c(1, 1))
plot(temp$Projected ~ temp$Date, type = "l", ylim=c(25, 35), lwd = 2, main = "TR", xlab = "Date", ylab = "Temperature (oC)")
lines(temp$Thres~temp$Date, col = 'orange', lty = 2)
lines(smooth_tl ~ temp$Date, col = 'blue')
lines(smooth_tr ~ temp$Date, col = 'green')
lines(smooth_br ~ temp$Date, col = 'red')
abline(v=17)


