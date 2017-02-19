# - - - - - - - - - - - - - - - - - 
exercise.Chapter9.4a <- function() {
  
  library(dlm)
  library(KFAS)
  
  # constants - - - - - -
  sFILENAME <- "UKfrontrearseatKSI.txt"
  # - - - - - - - - - - - 
  
  # read data
  dfData <- read.table(sFILENAME, skip=1)
  colnames(dfData) <- c("agDriver", "agFront", "agRear", "agKilo", "agPetrol")
  dfData$abSeatbelt <- 0
  dfData$abSeatbelt[170:nrow(dfData)] <- 1
  tsData <- ts(dfData, frequency=12, start=c(1969, 1))
  n      <- length(tsData)
  
  cat("\n")
  
  cat("estimated by KFAS package:\n")
  
  oModel.KFAS <- SSModel(
    cbind(agFront, agRear) ~ SSMtrend(degree=1, type="distinct", Q=matrix(NA, nrow=2, ncol=2))
     + SSMseasonal(12, sea.type="dummy")
     -1 + agKilo + agPetrol + abSeatbelt, 
    H=matrix(NA, nrow=2, ncol=2),
    data = tsData
  )
  ## print(oModel.KFAS)   # -> states are (level*2, seasondummy12*2 ...)
  ## print(t(oModel.KFAS$Z[,,1]))
  
  # Fitting
  oFitted.KFAS   <- fitSSM(
    oModel.KFAS, 
    inits  = rep(0, 8), 
    method = "BFGS"
  )
  stopifnot(oFitted.KFAS$optim.out$convergence == 0)
  
  cat("estimated H matrix:\n")
  print(oFitted.KFAS$model$H[,,1])
  cat("\n")
  
  cat("estimated Q matrix:\n")
  print(oFitted.KFAS$model$Q[,,1])
  cat("\n")
  
  # Filtering & Smoothing
  oEstimated.KFAS <- KFS(oFitted.KFAS$model, smoothing=c("state", "mean", "disturbance"))
  
  cat("correlation between level disturbance:\n")
  print(cor(oEstimated.KFAS$etahat[,1], oEstimated.KFAS$etahat[,2]))
  cat("\n")
  
  cat("regresssion coefficients:\n")
  print(as.matrix(oEstimated.KFAS$alphahat[1,1:6]))
  cat("\n")
  
  # plot 
  cat("making Figure 9.1 ... \n")
  plot(tsData, main = "Figure 9.1")
  
  cat("making Figure 9.2 ... \n")
  plot(
    as.numeric(oEstimated.KFAS$etahat[,2]), 
    as.numeric(oEstimated.KFAS$etahat[,1]), 
    xlab = "eta_Rear", ylab = "eta_Front", 
    main = "Figure 9.2"
  )
  abline(h=0)
  abline(lm(oEstimated.KFAS$etahat[,1] ~ oEstimated.KFAS$etahat[,2]))
  
  cat("making Figure 9.3 ... \n")
  plot(
    oEstimated.KFAS$alphahat[,7:8], 
    main = "Figure 9.3"
  )
  
  cat("making Figure 9.4 ... \n")
  plot(
    as.numeric(oEstimated.KFAS$alphahat[,8]), 
    as.numeric(oEstimated.KFAS$alphahat[,7]), 
    xlab = "level_Rear", ylab = "level_Front", 
    main = "Figure 9.4"
  )
  abline(lm(oEstimated.KFAS$alphahat[,7] ~ oEstimated.KFAS$alphahat[,8]))
}
# - - - - - - - - - -
exercise.Chapter9.4b <- function() {
  
  library(dlm)
  library(KFAS)
  
  # constants - - - - - -
  sFILENAME <- "UKfrontrearseatKSI.txt"
  # - - - - - - - - - - - 
  
  # read data
  dfData <- read.table(sFILENAME, skip=1)
  colnames(dfData) <- c("agDriver", "agFront", "agRear", "agKilo", "agPetrol")
  dfData$abSeatbelt <- 0
  dfData$abSeatbelt[170:nrow(dfData)] <- 1
  tsData <- ts(dfData, frequency=12, start=c(1969, 1))
  n      <- length(tsData)
  
  cat("\n")
  
  cat("estimated by KFAS package:\n")
  
  oModel.KFAS <- SSModel(
    cbind(agFront,agRear) ~ -1 
       + agPetrol 
       + agKilo
       + SSMregression(~-1+abSeatbelt, data=tsData, index=1)
       + SSMseasonal(period=12, sea.type='dummy')
       + SSMcustom(
           Z     = diag(2),
           T     = diag(2),
           R     = matrix(1,nrow=2,ncol=1),
           Q     = matrix(1),
           P1inf = diag(2)
         )
       ,
    data = tsData,
    H = matrix(NA, ncol=2, nrow=2)
  )
  # see manual of KFS()
  funUpdate <- function(pars, model, ...){
    model$H[,,1]   <- exp(0.5*pars[1:2])
    model$H[1,2,1] <- model$H[2,1,1] <- tanh(pars[3])*prod(sqrt(exp(0.5*pars[1:2])))
    model$R[28:29] <- exp(pars[4:5])
    ## print(model$H[,,1])
    ## print(model$R[28:29])
    model
  }
  # fitting
  oFitted.KFAS   <- fitSSM(
    oModel.KFAS, 
    inits      = c(-7, -7, 1, -1, -3),  # i'm not sure how we can get it ...
    updatefn  = funUpdate, 
    method    = "BFGS"
  )
  stopifnot(oFitted.KFAS$optim.out$convergence == 0)
  
  cat("estimated H matrix:\n")
  print(oFitted.KFAS$model$H[,,1])
  cat("\n")
  
  cat("estimated R matrix for custom disturbance:\n")
  print(oFitted.KFAS$model$R[28:29,,1])
  agQMatrix <- oFitted.KFAS$model$R[28:29,,1]%*%t(oFitted.KFAS$model$R[28:29,,1])
  cat("\n")
  
  cat("shown as Q matrix:\n")
  print(agQMatrix)
  cat("\n")
  
  # Filtering & Smoothing
  oEstimated.KFAS <- KFS(oFitted.KFAS$model, smoothing=c("state", "mean", "disturbance"))
  
  cat("regression coefficients:\n")
  print(as.matrix(oEstimated.KFAS$alphahat[1,1:5]))
  cat("\n")
  
  # plot
  cat("making Figure 9.6 ... \n")
  plot(
    as.numeric(oEstimated.KFAS$alphahat[,"custom2"]), 
    as.numeric(oEstimated.KFAS$alphahat[,"custom1"]), 
    main = "Figure 9.6", 
    xlab = "custom2", ylab="custom1"
  )
  cat("making Figure 9.7 ... \n")
  plot(
    oEstimated.KFAS$alphahat[,c("custom1", "custom2")], 
    main = "Figure 9.7"
  )
  cat("making Figure 9.8 ... \n")
  plot(
    oEstimated.KFAS$alphahat[,"custom1"] 
    + oEstimated.KFAS$alphahat[,"abSeatbelt.agFront"] * tsData[,"abSeatbelt"], 
    main = "Figure 9.8a", 
    ylab = NULL, xlab=NULL
  )
  legend("top", legend="level+intervention front", lty=1, cex=0.8)
  plot(
    oEstimated.KFAS$alphahat[,"custom2"], 
    main = "Figure 9.8b", 
    ylab = NULL, xlab=NULL
  )
  legend("top", legend="level rear", lty=1, cex=0.8)
  cat("making Figure 9.9 ... \n")
  plot(
    oEstimated.KFAS$alphahat[,"sea_dummy1.agFront"], 
    main = "Figure 9.9a", 
    ylab = NULL, xlab=NULL
  )
  legend("top", legend="seasonal front", lty=1, cex=0.8)
  
  plot(
    oEstimated.KFAS$alphahat[,"sea_dummy1.agRear"], 
    main = "Figure 9.9b", 
    ylab = NULL, xlab=NULL
  )
  legend("top", legend="seasonal rear", lty=1, cex=0.8)
}
# - - - - - - - - - -
exercise.Chapter10 <- function() {
  
  set.seed(123123)
  
  tsRandom <- ts(rnorm(200))
  sub.plotFigure_twoseries(
    tsRandom, 
    NULL, 
    "Figure10.1", 
    "random process"
  )
  sub.plotFigure_ACF (
    tsRandom, 
    12, 
    "Figure 10.2", 
    "ACF - random process"
  )
  
  tsRW <- ts(cumsum(rnorm(200)))
  sub.plotFigure_twoseries(
    tsRW, 
    NULL, 
    "Figure10.3", 
    "random walk"
  )
  sub.plotFigure_ACF (
    tsRW, 
    12, 
    "Figure 10.4", 
    "ACF - random walk"
  )
  
  tsMA <- ts(arima.sim(model=list(ma=0.5), n=200))
  sub.plotFigure_twoseries(
    tsMA, 
    NULL, 
    "Figure10.5", 
    "MA(1) process"
  )
  sub.plotFigure_ACF (
    tsMA, 
    12, 
    "Figure 10.6", 
    "ACF - MA(1) process"
  )
  
  tsAR <- ts(arima.sim(model=list(ar=0.5), n=200))
  sub.plotFigure_twoseries(
    tsAR, 
    NULL, 
    "Figure10.7", 
    "AR(1) process"
  )
  sub.plotFigure_ACF (
    tsAR, 
    12, 
    "Figure 10.8", 
    "ACF - AR(1) process"
  )
  
  tsARMA <- ts(arima.sim(model=list(ar=0.5, ma=0.5), n=200))
  sub.plotFigure_twoseries(
    tsARMA, 
    NULL, 
    "Figure10.9", 
    "ARMA(1,1) process"
  )
  sub.plotFigure_ACF (
    tsARMA, 
    12, 
    "Figure 10.10", 
    "ACF - ARMA(1,1) process"
  )
}
# = = = = = = = = = - - = = = 

# png(width = 480, height = 300)

# exercise.Chapter1 ()
#  
# exercise.Chapter2.1 ()
# exercise.Chapter2.2 ()
# exercise.Chapter2.3 () 
# 
# exercise.Chapter3.1 ()
# exercise.Chapter3.2 ()
# exercise.Chapter3.3 ()
# exercise.Chapter3.4a ()
# exercise.Chapter3.4b ()  
# 
# exercise.Chapter4.1 ()
exercise.Chapter4.2 ()
# exercise.Chapter4.3 () 
# exercise.Chapter4.4 ()
#  
# exercise.Chapter5.1a ()
# exercise.Chapter5.1b ()
# exercise.Chapter5.2 ()
#  
# exercise.Chapter6.1 ()
# exercise.Chapter6.2 ()
# 
# exercise.Chapter7.1 ()
# exercise.Chapter7.2 ()
# exercise.Chapter7.3 ()
# exercise.Chapter7.4 ()
# 
# exercise.Chapter8.6a ()
# exercise.Chapter8.6b () 
# exercise.Chapter8.6c () 
# exercise.Chapter8.7 () 
# 
# exercise.Chapter9.4a () 
# exercise.Chapter9.4b () 
# 
# exercise.Chapter10 () 

# dev.off()
