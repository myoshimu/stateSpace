df <- read.table("data/UKdriversKSI.txt", skip=1)
d <- ts(log(df[,1]))

# linear regression
tslm <- lm(d ~ time(d))
print(tslm)

plot(d, type = "p",pch  = 3, xlab = NA, ylab = NA, main = "Figure 1.1")
abline(tslm, col="red")
#legend("top", legend = c("log UK drivers KSI against time (in months)", "regression line"),  lty = c(-1, 1),    pch = c(3, -1), 
#col = c("black", "red"),cex = .8)

# time series
plot(d, type = "n", xlab = NA, ylab = NA,  main = "Figure 1.2")
lines(d, col = "black")
legend("top", legend = "log UK drivers KSI", lty= 1, col= c("black", "red"),cex= .8)


plot(d - predict(tslm),type="l",lty= 2,col= "green",xlab= NA, ylab  = NA, main  = "Figure 1.3")
abline(h=0, col="gray")
legend("top", lty= 2, legend = "residuals",  col= "green",cex = .8  )
  

#ACF
acf(ts(as.ts(rnorm(length(d))), frequency=1), lag.max = 15,lwd = 12, col = "green"
    , main = "Figure 1.4", xlab = NA, ylab = NA , xlim = c(1,15))
legend("top", legend = "ACF-random residuals", lty= 1, col = "green", cex = .8) 


acf(ts(d - predict(tslm), frequency=1), lag.max = 15,lwd = 12, col = "green"
    , main = "Figure 1.5", xlab = NA, ylab = NA , xlim = c(1,15))
legend("top", legend = "ACF-regression residuals", lty= 1, col = "green", cex = .8) 
