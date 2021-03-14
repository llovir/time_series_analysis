podaci <- read.csv("putnici.csv")
str(podaci)
attach(podaci)
dim(podaci) #144

install.packages("tseries")
install.packages("fGarch")
install.packages("TSA")
install.packages("aTSA")
install.packages("FinTS")
install.packages("astsa")
install.packages("forecast")

library("tseries")
library("TSA")
library("fGarch")
library("aTSA")
library("FinTS")
library("astsa")
library("forecast")
acf <- stats::acf
arima <- stats::arima

putnik <- podaci$Airlinepassengers
putnik
u <- ts(putnik, start=c(1949,1), end=c(1960,12), frequency = 12)
u
length(u)
v<-window(u, end=c(1956,12))
length(v)
v
summary(v)
sd(v)
par(mfrow = c(1,1))
plot(u, typecol="blue", xlab = "Date", ylab = "Number of passengers",main="The number of airline passengers per month")
plot(v,type="o", xlab="Date",ylab="Number of passengers",col="darkblue",main="The number of airline passengers per month")
##increasing trend
lines(v)
summary(v)
(minimalan<-Month[putnik==104])
(maximalan<-Month[putnik==413])
str(podaci)

#boxplot by month:

#comparison of months:
v7 <- window(v, start = c(1949,7), freq = T)
v11 <- window(v, start = c(1949,11), freq = T)
mean(v7)
mean(v11)

boxplot(v7,v11, col="darkblue",names = c("July","November"))

acf(v, main="Autocorrelation function", xlab='Lag', ylab='ACF', lag.max=72)
plot(v)

boxplot(v~cycle(v), col='darkblue', names=c('January','February','March','April','May','June','July','August',
                                            'September','October','November','December'), xlab='Month', ylab='Number of passengers',main='Number of passengers per month')

putnici7mj<-window(v, start=c(1949,7), freq=T)
putnici11mj<-window(v,start=c(1949,11),freq=T)

plot(putnici7mj, ylim=c(50,460), col='darkblue', lwd=2, xlab='Date', ylab='Number of passengers', main='Number of passengers in july and november')
lines(putnici11mj, col='red', lwd=2)
legend("topleft", legend=c('July','November'), col=c('darkblue','red'), lty=c(1,1))

# log
putnik <- ts(podaci$Airlinepassengers, start = c(1949, 1), freq = 12)
putnikM <- window(putnik, end = c(1956, 12))
length(putnikM)
logputnik<-log(putnik)
plot(putnikM, xlab = "Date", ylab = "Number of passengers", col = "darkblue")
points(putnikM, col = "red", pch = 16, cex = 0.5)
adf.test(putnikM)
plot(diff(putnikM))
logputnikM <- log(putnikM)
plot(logputnikM,xlab = "Date", ylab = "log (Number of passengers)", col = "darkblue")
points(logputnikM, col = "red", pch = 16, cex = 0.5)
acf(logputnikM, main="Autocorrelation  function", xlab="Lag",lag.max = 72)
acf(logputnikM, lag.max = 72)
#Differencing
DlogputnikM <- diff(logputnikM)
plot(DlogputnikM , xlab = "Date", ylab = "log- Number of passengers", col = "darkblue")
points(DlogputnikM, col = "red", pch = 16, cex = 0.5)
acf(DlogputnikM, main="Autocorrelation  function", xlab="Lag",lag.max = 72)
adf.test(DlogputnikM)
#Differencing lag 12
diff_ss12 <- diff(DlogputnikM,lag=12)
adf.test(diff_ss12) 

plot(diff_ss12, main="", xlab="Date", ylab = "log- Number of passengers")
acf(diff_ss12, main="Autocorrelation function", lag.max=72, xlab="Lag")
plot(diff_ss12, col='darkblue', main='', xlab='Date', ylab='Number of passengers')

################################

library("TSA")

#AIC 
najboljiAICs <- function(armax, mamax, podaci) {
  redAR <- c()
  redMA <- c()
  redSAR <- c()
  redSMA <- c()
  AICvrijednost <- c()
  for (i in 0:armax) {
    for (j in 0:mamax) {
      for (k in 0:armax) {
        for (l in 0:mamax) {
          fitaic <- tryCatch(AIC(arima(podaci, c(i, 0, j), seasonal = c(k, 0, l))), error = function(e) NaN) 
          redAR <- c(redAR,i)
          redMA <- c(redMA, j)
          redSAR <- c(redSAR, k)
          redSMA <- c(redSMA, l)
          AICvrijednost <- c(AICvrijednost, fitaic)
        }
      }
    }
  }
  rez <- data.frame(p = redAR, q = redMA, P = redSAR, Q = redSMA, AICvrijednosti = AICvrijednost)
  rez <- rez[order(rez$AICvrijednosti), ]
  return(rez[1:min(10,length(redAR)), ])
}

najboljiAICs(2,2,diff_ss12)


# p q P Q AICvrijednosti
# 0 1 0 1      -289.5586
# 2 1 0 1      -288.3886
# 1 0 0 1      -288.3126
# 1 2 0 1      -288.0699
# 0 2 0 1      -287.9115
# 0 1 1 1      -287.6364
# 0 1 0 2      -287.6330
# 1 1 0 1      -286.7063
# 2 0 0 1      -286.4467
#2 1 0 2      -286.3992
#log
#BIC
najboljiBICs <- function(armax, mamax, podaci) {
  redAR <- c()
  redMA <- c()
  redSAR <- c()
  redSMA <- c()
  BICvrijednost <- c()
  for (i in 0:armax) {
    for (j in 0:mamax) {
      for (k in 0:armax) {
        for (l in 0:mamax) {
          fitbic <- tryCatch(BIC(arima(podaci, c(i, 0, j), seasonal = c(k, 0, l))), error = function(e) NaN) 
          redAR <- c(redAR,i)
          redMA <- c(redMA, j)
          redSAR <- c(redSAR, k)
          redSMA <- c(redSMA, l)
          BICvrijednost <- c(BICvrijednost, fitbic)
        }
      }
    }
  }
  rez <- data.frame(p = redAR, q = redMA, P = redSAR, Q = redSMA, BICvrijednosti = BICvrijednost)
  rez <- rez[order(rez$BICvrijednosti), ]
  return(rez[1:min(10,length(redAR)), ])
}

najboljiBICs(2, 2, diff_ss12)

#log
# p q P Q BICvrijednosti
# 11 0 1 0 1      -279.8832
# 29 1 0 0 1      -278.6372
# 20 0 2 0 1      -275.8173
# 14 0 1 1 1      -275.5422
# 12 0 1 0 2      -275.5388
# 13 0 1 1 0      -274.7572
# 38 1 1 0 1      -274.6121
# 56 2 0 0 1      -274.3525
# 32 1 0 1 1      -274.2356
# 30 1 0 0 2      -274.2347

#candidates for model

model<- sarima(logputnikM, 0,1,1 ,0, 1, 1, 12)
model$fit
BIC(model$fit)
coef(model$fit)
confint(model$fit)
Box.test(model$fit$residuals, type = "Ljung-Box")
#uncorrelated
# var= 0.001515:  log likelihood = 148.76,  aic = -291.53
# BIC(model$fit)= -284.2695
 
model0<- sarima(logputnikM, 2,1,1,0,1,1, 12)
model0$fit
BIC(model0$fit)
confint(model0$fit)
Box.test(model0$fit$residuals, type = "Ljung-Box")

model1<- sarima(logputnikM, 1,1,0,0,1,1, 12)
model1$fit
BIC(model1$fit)
confint(model1$fit)
Box.test(model1$fit$residuals, type = "Ljung-Box")

model2<- sarima(logputnikM, 1,1,2,0,1,1, 12)
confint(model2$fit)
model2$fit
BIC(model2$fit)
Box.test(model2$fit$residuals, type = "Ljung-Box")

model3<- sarima(logputnikM,0,1,2,0,1,1, 12) 
confint(model3$fit)
model3$fit
BIC(model3$fit)
Box.test(model3$fit$residuals, type = "Ljung-Box")

model4<- sarima(logputnikM,0,1,1,1,1,1, 12)
confint(model4$fit)
model4$fit
BIC(model4$fit)
Box.test(model4$fit$residuals, type = "Ljung-Box")

model5<- sarima(logputnikM,0,1,1,0,1,2, 12) 
confint(model5$fit)
model5$fit
BIC(model5$fit)
Box.test(model5$fit$residuals, type = "Ljung-Box")

model6<- sarima(logputnikM,1,1,1,0,1,1, 12) ###final model
confint(model6$fit)
model6$fit
BIC(model6$fit)
Box.test(model6$fit$residuals, type = "Ljung-Box") #model with best variance, AIC and BIC

model7<- sarima(logputnikM,2,1,0,0,1,1, 12) 
model7$fit
BIC(model7$fit)
Box.test(model7$fit$residuals, type = "Ljung-Box")

model8<- sarima(logputnikM,2,1,1,0,1,2, 12) 
model8$fit
BIC(model8$fit)
Box.test(model8$fit$residuals, type = "Ljung-Box")

model9<- sarima(logputnikM,1,1,0,1,1,1, 12) 
model9$fit
BIC(model9$fit)
Box.test(model9$fit$residuals, type = "Ljung-Box")

model11<- sarima(logputnikM,1,1,0,0,1,2, 12) 
model11$fit
BIC(model11$fit)
Box.test(model11$fit$residuals, type = "Ljung-Box")
#################################################

####residuals
reziduali<-model6$fit$residuals
plot(reziduali, main='Residuals', xlab='Date', ylab='Residuals')
par(mfrow=c(1,1))
hist(reziduali, main='Histogram of residulas ', xlab='Residuals', ylab='Frequence')
qqnorm(reziduali,xlab="Theoretical quantile",ylab="Quantile of empirical distribution", main="QQ-Plot")
qqline(reziduali,col="red",lwd=1)
#test for normality
shapiro.test(reziduali)
acf(reziduali,main='Autocorrelation function of residuals', xlab='Lag', ylab='ACRF',lag.max=36)
Box.test(reziduali, type = "Ljung-Box")

modelfix <- arima(x = logputnikM, order = c(1, 1, 1), seasonal = list(order = c(0,1,1), period = 12), fixed = c(0, NA, NA))
modelfix
confint(modelfix) 
reziduali<-modelfix$residuals
plot(reziduali, main='Residuals', xlab='Date', ylab='Residuals')
par(mfrow=c(1,1))
hist(rez, main='Histogram of residuals ', xlab='Residuals', ylab='Frequence')
qqnorm(reziduali,xlab="Theoretical quantile",ylab="Quantile of empirical distribution", main="QQ-Plot")
qqline(reziduali,col="red",lwd=1)
shapiro.test(reziduali) 
acf(reziduali,main='Autocorrelation function of residuals', xlab='Lag', ylab='ACRF',lag.max=36)
Box.test(rez, type = "Ljung-Box") 


# Prediction
razlika <- length(u) - length(v)
razlika 

library(forecast)
par(mfrow=c(1,1))
plot(forecast(modelfix,48), main="Prediction of model", xlab="Year", ylab="Number of passengers")
plot(forecast(modelfix,48),xlim=c(1949,1960),type="l",main="Prediction od model",ylab="Number of passengers",xlab="Year") 
pred <- forecast(modelfix,h = 48)
str(pred)

plot(logputnik, xlab="Date", ylab="Airline passengers",xlim=c(1949,1960))
polygon(c(time(pred$lower), rev(time(pred$upper))), c(pred$lower[, 2], rev(pred$upper[, 2])), col = rgb(0,0,1,0.1), border = FALSE)
polygon(c(time(pred$lower), rev(time(pred$upper))), c(pred$lower[, 1], rev(pred$upper[, 1])), col = rgb(0,0,1,0.2), border = FALSE)
points(pred$mean,col="red",pch=16,cex=0.8)
points(logputnik,col="blue",pch=16,cex=0.8)
legend(cex=1,"topleft", c("real", "prediction"), col = c("blue", "red"), pch = c(16, 16))