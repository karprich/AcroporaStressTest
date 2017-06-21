#Temperature Analysis
# ? strptime
temp <- read.csv("/Users/Rich/Documents/GradSchool/Internship/Lab/Data/AcroporaStressTest/Data/TemperatureData/Rdata1.csv", header = T)
temp$Date <- as.POSIXct(temp$Date, format = "%H:%M %m/%d/%Y")

#Add Threshhold Temperature of 32
temp$Thres <- 32

#omit sensor
temp[temp$Date>as.POSIXct("2017-06-03 11:40:00") & temp$Date<as.POSIXct("2017-06-05 16:00:00"), "Tmpx17.bottom.Right"]<-NA
temp[temp$Date>as.POSIXct("2017-06-19 12:40:00") & temp$Date<as.POSIXct("2017-06-19 16:40:00"), "Tmpx19.top.left"]<-NA



#Write Figure 
#png(filename = "Output/Figures/tempthrough190617.png", width=5, height = 5, units = "in", res = 300)
# for multipanel plot 
library(zoo)
par(mfrow = c(2, 2))
#Plot TL
smooth_tl <- rollapply(data = temp$Tmpx19..6b.top.left, width = 100, FUN = mean, na.rm = T, fill = NA)
plot(temp$Projected ~ temp$Date, type = "l", ylim=c(25, 35), lwd = 2, main = "TL", xlab = "Date", ylab = "Temperature (oC)")
lines(temp$Thres~temp$Date, col = 'orange', lty = 2)
lines(temp$Tmpx19..6b.top.left ~ temp$Date, lwd = (.5))
lines(smooth_tl ~ temp$Date, col = 'red')
#Plot TR
smooth_tr <- rollapply(data = temp$Tmpx21...7b.top.right, width = 100, FUN = mean, na.rm = T, fill = NA)
plot(temp$Projected ~ temp$Date, type = "l", ylim=c(25, 35), lwd = 2, main = "TR", xlab = "Date", ylab = "Temperature (oC)")
lines(temp$Thres~temp$Date, col = 'orange', lty = 2)
lines(temp$Tmpx21...7b.top.right ~ temp$Date, lwd = (.5))
lines(smooth_tr ~ temp$Date, col = 'red')
#Plot BL
smooth_bl <- rollapply(data = temp$Tmpx23.6A.Bottom.Left, width = 100, FUN = mean, na.rm = T, fill = NA, xlab = "Date", ylab = "Temperature (oC)")
plot(temp$Projected ~ temp$Date, type = "l", ylim=c(25, 35), lwd = 2, main = "BL", xlab = "Date", ylab = "Temperature (oC)")
lines(temp$Thres~temp$Date, col = 'orange', lty = 2)
lines(temp$Tmpx23.6A.Bottom.Left ~ temp$Date, lwd = (.5))
lines(smooth_bl ~ temp$Date, col = 'red')
#Plot BR
smooth_br <- rollapply(data = temp$Tmpx17.bottom.Right, width = 100, FUN = mean, na.rm = T, fill = NA, xlab = "Date", ylab = "Temperature (oC)")
plot(temp$Projected ~ temp$Date, type = "l", ylim=c(25, 35), lwd = 2, main = "BR")
lines(temp$Thres~temp$Date, col = 'orange', lty = 2)
lines(temp$Tmpx17.bottom.Right ~ temp$Date, lwd = (.5))
lines(smooth_br ~ temp$Date, col = 'red')
#dev.off()
