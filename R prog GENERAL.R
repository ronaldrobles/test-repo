library(installr)
updateR(F, T, T, F, T, F, T) 

devtools::install_github('rstudio/shinyapps')

######### R Programming clase1 ###################

### Objects

?vector()
vector("numeric", length = 10)

as.vector(x = 1:3, mode = "numeric")
as.vector(x = 1:3, mode = "character")
as.vector(x = 1:3, mode = "logical")

x <- c(a = 1, b = 2)
is.vector(x)
as.vector(x)


x <- c(a = c(1,2,3), b = 2)

b<- 3L
b
is.numeric(b)
str(b)
is.integer(b)

NaN

#### Attributes

x <- c(a = 1, b = 2)
names(x)
dimnames(x)
dim(x)
class(x)
length(x)


x <- c(0.5, 0.6)       ## numeric
x <- c(TRUE, FALSE)    ## logical
x <- c(T, F)           ## logical
x <- c("a", "b", "c")  ## character
x <- 9:29              ## integer
x <- c(1+0i, 2+4i)     ## complex


## Coercion : character<numeric<logical

y <- c(1.7, "a")   ## character
y <- c(TRUE, 2)    ## numeric
y <- c("a", TRUE)  ## character
class(y)

### Explicit Coercion

x <- 0:6
class(x)
as.numeric(x)
as.logical(x)
as.character(x)
as.complex(x)
x

###### Matrices (vectors with dimension attribute)

m <- matrix(nrow = 2, ncol = 3)
m
dim(m)
attributes(m)

## Creating matrices

m <- matrix(1:20, nrow = 5, ncol = 4)
m

## Creating matrices

m <- 1:20
m
dim(m)
class(m)
dim(m) <- c(4, 5)
class(m)

## Creating matrices

x <- 1:3
y <- 10:12

cbind(x, y)
rbind(x, y)

######### Lists

x <- list(1, "a", c(1:25), TRUE, 1 + 4i) 
x

#### Factors

x <- factor(c("yes", "yes", "no", "yes", "no"), levels = c("yes", "no"))
x
unclass(x)

y <- factor(c(1, 2, 1, 1, 2, 2, 1, 2), levels = c(1,2),labels = c("red", "blue"))

y

##### Missing values

x <- list(1, "a", c(1:25), TRUE, 1 + 4i, NaN) 
is.na(x)
is.nan(x)

#### Data Frames

x <- data.frame(foo = 1:4, bar = c(T, T, F, F))
nrow(x)
ncol(x)
rownames(x)
colnames(x)

########## Names###

##dataframes
x <- 1:3
names(x)
#NULL
names(x) <- c("foo", "bar", "norf")
x
names(x)

##lists
x <- list(a = 1, b = 2, c = 3)
x

##matrices

m <- matrix(1:4, nrow = 2, ncol = 2)

dimnames(m)
dim(m)
class(m)
dimnames(m) <- list(c("a", "b"))

m
dimnames(m) <- list(c("a", "b"), c("c", "d"))



######## R Programming Quiz 1 #################

tabla <- read.table(file = "clipboard", sep = "\t", header=TRUE)
View(tabla)

colnames(tabla)

tabla[1:2,]

str(tabla)

tabla[152:153,]

tabla[47,]
tabla$Ozone[47]

is.na(tabla$Ozone)
sum(is.na(tabla$Ozone))

mean(tabla$Ozone, na.rm=TRUE)



newdata <- na.omit(tabla)
above<- newdata[newdata$Ozone>31,]
above2<- above[above$Temp>90,]
summary(above2)

above<- tabla[tabla$Ozone>31,]
above2<- above[above$Temp>90,]



summary(tabla[tabla$Month==6,])

summary(tabla[tabla$Month==5,])


## http://www.statmethods.net/input/missingdata.html


######### R Programming clase 2 ###################


x <- c("a", "b", "c", "c", "d", "a")
x[1]
#[1] "a"

x[2]
#[1] "b"

x[1:4]
#[1] "a" "b" "c" "c"

x[x > "a"]
#[1] "b" "c" "c" "d"

###Logical index
u <- x > "a"
u
#[1] FALSE TRUE TRUE TRUE TRUE FALSE

x[u]
#[1] "b" "c" "c" "d"

x[!u]
# [1] "a" "a"

### Subsetting a Matrix

x <- matrix(1:6, 2, 3)

##       [,1] [,2] [,3]
## [1,]    1    3    5
## [2,]    2    4    6

#x[row, col]

x[1, 2]
#[1] 3

x[2, 1]
#[1] 2

x[1, ]
#[1] 1 3 5

x[, 2]
#[1] 3 4


### drop (mantiene la clase )

x[1, 2, drop = FALSE]

##       [,1]
## [1,]    3

x[1, , drop = FALSE]

#       [,1] [,2] [,3]
# [1,]    1    3    5

### Subsetting a List

x <- list(foo = 1:4, bar = 0.6)

x[1]
# $foo
# [1] 1 2 3 4

x[[1]]
# [1] 1 2 3 4

x$foo


x$bar
# [1] 0.6

x[["bar"]]
# [1] 0.6

x[[2]]
# [1] 0.6

x["bar"]
# $bar
# [1] 0.6

x[[x$bar]]
# [1] 1 2 3 4


x <- list(foo = 1:4, bar = 0.6, baz = "hello")
x[c(1, 3)]
# $foo
# [1] 1 2 3 4
# $baz
# [1] "hello"


x <- list(foo = 1:4, bar = 0.6, baz = "hello")
name <- "foo"
x[name]
x[[name]]
x$name


x <- list(a = list(10, 12, 14), b = c(3.14, 2.81))

x[[c(1, 3)]]
# [1] 14

x[[1]][[3]]
# [1] 14

x[[c(2, 1)]]
# [1] 3.14


### Partial marching  [[ and $

x <- list(aardvark = 1:5)
x$a
# [1] 1 2 3 4 5

x[["a"]]
# NULL

x[["a", exact = FALSE]]
# [1] 1 2 3 4 5

### Removing NA values

x <- c(1, 2, NA, 4, NA, 5)
is.na(x)
bad <- is.na(x)
bad
x[!bad]
# [1] 1 2 4 5


# complete.cases (argumentos con la misma longitud)

x <- c(1, 2, NA, 4, NA, 6)
y <- c("a", "b", NA, "d", NA, "f")
good <- complete.cases(x, y)
good
# [1]  TRUE  TRUE FALSE  TRUE FALSE  TRUE
x[good]
# [1] 1 2 4 5
y[good]
# [1] "a" "b" "d" "f"
x[!good]
# [1] NA NA


airquality[1:6, ]

good <- complete.cases(airquality)
good
airquality[good, ]
airquality[good, ][1:6, ]




### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
## http://www.statmethods.net/input/

# recode 99 to missing for variable v1
# select rows where v1 is 99 and recode column v1 
mydata$v1[mydata$v1==99] <- NA 

x <- c(1,2,NA,3)
mean(x) # returns NA
mean(x, na.rm=TRUE) # returns 2

complete.cases(airquality)
airquality[complete.cases(airquality),]
airquality[!complete.cases(airquality),]

na.omit(airquality) 
newdata <- na.omit(airquality) 

#You can go beyond pairwise of listwise deletion of missing values
#through methods such as multiple imputation. Good implementations
#that can be accessed through R include Amelia II, Mice, and mitools. 

# enter data using editor 
mydata <- data.frame(age=numeric(0), gender=character(0), weight=numeric(0))
mydata <- edit(mydata)

mydat$v1 <- factor(mydat$v1,
                    levels = c(1,2,3),
                    labels = c("red", "blue", "green")) 

mydata$v1 <- ordered(mydata$y,
                     levels = c(1,3, 5),
                     labels = c("Low", "Medium", "High")) 

#################   EXPLORATORY DATA ANALYSIS    ################

eq <- rnorm(100)
hist(eq)
y <- rnorm(100)
plot(eq,y)
z <- rnorm(100)
plot(eq,z)
par(mar=c(6,6,2,2))
plot(eq,y, pch=20)
plot(eq,y, pch=2)
plot(eq,y, pch=20)
plot(eq,y, pch=15)
plot(eq,y, pch=20)
plot(eq,y, pch=20)

plot(eq,z, pch=20)
title("Scaterrplot")
text(-2,-2,"Label")
legend("topleft", legend = "Data", pch=20)
fit <- lm(y~eq)
abline(fit)
abline(fit, lwd=3, col="blue")

plot(eq,y, xlab= "Weight", ylab = "Height", main ="Scatterplot", pch=20)
legend("topright", legend = "Data", pch=20)
fit <- lm(y~eq)
abline(fit, lwd=3, col="red")

z<- rpois(100,2)
par(mfrow = c(2,1))
plot(eq, y, pch=19)
plot(eq, z, pch=20)
par("mar")
par(mar= c(2,2,1,1))
plot(eq, y, pch=19)
plot(eq, z, pch=20)
par(mfrow = c(1,2))
plot(eq, y, pch=19)
plot(eq, z, pch=20)
par(mar= c(4,4,1,1))
plot(eq, y, pch=19)
plot(eq, z, pch=20)


###### PROJECT 1 #######
getwd()
setwd("C:/Users/Ronald/GitHub/ExData_Plotting1")

bigtable <- read.csv("household_power_consumption.txt", header =T, sep=";", na.string= "?")
table1<- bigtable[bigtable$Date == "1/2/2007",]
table2<- bigtable[bigtable$Date == "2/2/2007",]
table<- rbind(table1, table2)
str(table)


#str(strptime(table$Date, format= "%d/%m/%Y"))
#table$Date<- strptime(table$Date, format= "%d/%m/%Y")
#str(table)
#str(strptime(table$Time, format= "%H:%M:%S"))

paste_date_time <- paste(table$Date, table$Time)
str(paste_date_time)

str(strptime(paste_date_time, format= "%d/%m/%Y %H:%M:%S"))
nwDate <- strptime(paste_date_time, format= "%d/%m/%Y %H:%M:%S")

table2 <- cbind(nwDate, table)
str(table2)

##  for(i in c(4:7)) {  table2[,i] <- as.numeric(table2[,i]) }


Sys.getlocale("LC_TIME")
# [1] "Spanish_Peru.1252"

Sys.setlocale("LC_TIME", "English")


Sys.setlocale("LC_TIME", "Spanish_Peru.1252")

########  PLOT 1  #############

png(file = "plot1.png", bg = "transparent")
hist(table2$Global_active_power, col="red", main="Global Active Power", xlab="Global Active Power (kilowatts)")
dev.off()


########  PLOT 2  #############

png(file = "plot2.png", bg = "transparent")
plot(table2$nwDate, table2$Global_active_power, type = "l", ylab="Global Active Power (kilowatts)", xlab="")
dev.off()

########  PLOT 3  #############

png(file = "plot3.png", bg = "transparent")

plot(table2$nwDate, table2$Sub_metering_1, type = "l", ylab="Energy sub metering", xlab="")
points(table2$nwDate, table2$Sub_metering_2, type = "l",col= "red")
points(table2$nwDate, table2$Sub_metering_3, type = "l", col= "blue")
legend("topright", legend = colnames(table2[8:10]),lty=1, col= c("black", "red", "blue"))

dev.off()

########  PLOT 4  #############

png(file = "plot4.png", bg = "transparent")

par(mfrow= c(2,2))

plot(table2$nwDate, table2$Global_active_power, type = "l", ylab="Global Active Power", xlab="")
plot(table2$nwDate, table2$Voltage, type = "l", ylab="Voltage", xlab="datetime")

plot(table2$nwDate, table2$Sub_metering_1, type = "l", ylab="Energy sub metering", xlab="")
points(table2$nwDate, table2$Sub_metering_2, type = "l",col= "red")
points(table2$nwDate, table2$Sub_metering_3, type = "l", col= "blue")
legend("topright", legend = colnames(table2[8:10]),lty=1, col= c("black", "red", "blue"), bty="n")

plot(table2$nwDate, table2$Global_reactive_power, type = "l", ylab="Global Reactive Power", xlab="datetime")

dev.off()










elecdata <- fread("http://d396qusza40orc.cloudfront.net/exdata%2Fdata%2Fhousehold_power_consumption.zip", sep=";",
                  nrows=2880, na.string="?", skip=68076)








dtPOSIXct <- as.POSIXct(dtstring)

# extract time of 'date+time' (POSIXct) in hours as numeric
dtTime <- as.numeric(dtPOSIXct - trunc(dtPOSIXct, "days"))


########## Data PRoducts Quiz1 ##############

library(manipulate)
myPlot <- function(s) {
  plot(cars$dist - mean(cars$dist), cars$speed - mean(cars$speed))
  abline(0, s)
}

manipulate(myPlot, s = slider(0, 2, step = 0.1))  

manipulate(myPlot(s), s = slider(0, 2, step = 0.1))  

manipulate(myPlot(s), x.s = slider(0, 2, step = 0.1))  

manipulate(myPlot(s), slider = x(0, 2, step = 0.1))  

shinyapps::setAccountInfo(
  name="rrobles", 
  token="8F8882AAD490C299507703FAA765A051", 
  secret="gREsHhBuNkyEB26Go7Z8sS0OvPvzUZ3aXAOfadEF")






library(shiny)
library(ggplot2)

shinyServer(function(input, output) {
  
  dataset <- reactive(function() {
    diamonds[sample(nrow(diamonds), input$sampleSize),]
  })
  
  output$plot <- reactivePlot(function() {
    
    p <- ggplot(dataset(), aes_string(x=input$x, y=input$y)) + geom_point()
    
    if (input$color != 'None')
      p <- p + aes_string(color=input$color)
    
    facets <- paste(input$facet_row, '~', input$facet_col)
    if (facets != '. ~ .')
      p <- p + facet_grid(facets)
    
    if (input$jitter)
      p <- p + geom_jitter()
    if (input$smooth)
      p <- p + geom_smooth()
    
    print(p)
    
  }, height=700)
  
})



