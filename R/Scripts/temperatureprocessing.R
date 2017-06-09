#Temperature Analysis
# ? strptime
temp <- read.csv("/Users/Rich/Documents/GradSchool/Internship/Lab/Data/AcroporaStressTest/Data/TemperatureData/Rdata1.csv", header = T)
temp$Date <- as.POSIXct(temp$Date, format = "%H:%M %m/%d/%Y")


#omit sensor
temp[temp$Date>as.POSIXct("2017-06-04 02:00:00") & temp$Date<as.POSIXct("2017-06-05 16:00:00"), "Tmpx17.bottom.Right"]<-NA


# for multipanel plot 
library(zoo)
par(mfrow = c(2, 2))
smooth_tl <- rollapply(data = temp$Tmpx19..6b.top.left, width = 100, FUN = mean, na.rm = T, fill = NA)
plot(temp$Projected ~ temp$Date, type = "l", ylim=c(25, 35), lwd = 2, main = "TL", xlab = "Date", ylab = "Temperature (oC)")
lines(temp$Tmpx19..6b.top.left ~ temp$Date, lwd = (.5))
lines(smooth_tl ~ temp$Date, col = 'red')
smooth_tr <- rollapply(data = temp$Tmpx21...7b.top.right, width = 100, FUN = mean, na.rm = T, fill = NA)
plot(temp$Projected ~ temp$Date, type = "l", ylim=c(25, 35), lwd = 2, main = "TR", xlab = "Date", ylab = "Temperature (oC)")
lines(temp$Tmpx21...7b.top.right ~ temp$Date, lwd = (.5))
lines(smooth_tr ~ temp$Date, col = 'red')
smooth_bl <- rollapply(data = temp$Tmpx23.6A.Bottom.Left, width = 100, FUN = mean, na.rm = T, fill = NA, xlab = "Date", ylab = "Temperature (oC)")
plot(temp$Projected ~ temp$Date, type = "l", ylim=c(25, 35), lwd = 2, main = "BL", xlab = "Date", ylab = "Temperature (oC)")
lines(temp$Tmpx23.6A.Bottom.Left ~ temp$Date, lwd = (.5))
lines(smooth_bl ~ temp$Date, col = 'red')

smooth_br <- rollapply(data = temp$Tmpx17.bottom.Right, width = 100, FUN = mean, na.rm = T, fill = NA)
plot(temp$Projected ~ temp$Date, type = "l", ylim=c(25, 35), lwd = 2, main = "BR")
lines(temp$Tmpx17.bottom.Right ~ temp$Date, lwd = (.5))
lines(smooth_br ~ temp$Date, col = 'red')

