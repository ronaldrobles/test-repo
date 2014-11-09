bigtable <- read.csv("household_power_consumption.txt", header =T, sep=";", na.string= "?")
table1<- bigtable[bigtable$Date == "1/2/2007",]
table2<- bigtable[bigtable$Date == "2/2/2007",]
table<- rbind(table1, table2)

paste_date_time <- paste(table$Date, table$Time)
nwDate <- strptime(paste_date_time, format= "%d/%m/%Y %H:%M:%S")
table2 <- cbind(nwDate, table)

########  PLOT 3  #############

png(file = "plot3.png", bg = "transparent")

plot(table2$nwDate, table2$Sub_metering_1, type = "l", ylab="Energy sub metering", xlab="")
points(table2$nwDate, table2$Sub_metering_2, type = "l",col= "red")
points(table2$nwDate, table2$Sub_metering_3, type = "l", col= "blue")
legend("topright", legend = colnames(table2[8:10]),lty=1, col= c("black", "red", "blue"))

dev.off()