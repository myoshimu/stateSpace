# - - - - - - - - - - - - - - - - - 
exercise.Chapter5.1a <- function(){
  
  library(dlm)
  library(KFAS)
  
  # constants - - - - - -
  sFILENAME <- "UKdriversKSI.txt"
  gCKLogLik <- 0.4140728
  # - - - - - - - - - - - 
  
  # reading data
  dfData <- read.table(sFILENAME, skip=1)
  tsData <- ts(log(dfData[,1]), frequency=12, start=c(1969, 1))
  n      <- length(tsData)
  
  x      <- ts(as.numeric(1:n), frequency=12, start=c(1969, 1))
  
  # model specs
  nStateVar   <- 2
  nHyperParam <- 1
  
  # init value for variances (epsilon)
  gInit     <- c(var(tsData))  
  
  cat("\n")
  
  cat("estimated by dlm package:\n")
  funModel <- function(parm){
    m1 <- dlmModPoly(order = 1, dV=exp(parm[1]), dW=0)
    m2 <- dlmModReg(x, dV=0, addInt = FALSE)
    return( m1 + m2 )
  }
  # Fitting
  oMLE.DLM <- dlmMLE(
    tsData, 
    parm  = log(gInit), 
    build = funModel
  )
  stopifnot(oMLE.DLM$convergence == 0)
  oFitted.DLM    <- funModel(oMLE.DLM$par)
  # Filtering
  oFiltered.DLM <- dlmFilter(tsData, oFitted.DLM)
  agStdPredErr.DLM  <- residuals(oFiltered.DLM, type = c("standardized"), sd=FALSE)
  lPredErr.DLM      <- residuals(oFiltered.DLM, type = c("raw"), sd=TRUE)  
  # Smoothing
  oSmoothed.DLM <- dlmSmooth(tsData, oFitted.DLM)
  # LL
  gLLbyFun.DLM    <- -oMLE.DLM$value - 0.5*n*log(2*pi)
  gLLbyErr.DLM    <- sub.LogLik(lPredErr.DLM$res, lPredErr.DLM$sd^2, nStateVar)
  # report
  cat(sprintf("hat(sigma^2_epsilon) = %f\n", drop(oFitted.DLM$V)))
  cat(sprintf("hat(mu_1)            = %f\n", oSmoothed.DLM$s[2,1]))
  cat(sprintf("hat(beta_1)          = %f\n", oSmoothed.DLM$s[2,2]))
  cat("\n") 
  
  cat("estimated by KFAS package:\n")
  oModel.KFAS <- SSModel(
    tsData ~ SSMtrend(degree=1, Q=list(matrix(0))) + 
             SSMregression(~ -1 + x, Q=matrix(0)), 
    H=matrix(NA)
  )
  ## print(oModel.KFAS) -> states are (x, level)
  ## print(oModel.KFAS$Q) -> Q is (2,2) matrix of 0
  # Fitting
  oFitted.KFAS   <- fitSSM(
    oModel.KFAS, 
    inits  = log(gInit), 
    method = "BFGS"
  )
  stopifnot(oFitted.KFAS$optim.out$convergence == 0)
  # Filtering & Smoothing
  oEstimated.KFAS <- KFS(oFitted.KFAS$model)
  agStdPredErr.KFAS  <- rstandard(oEstimated.KFAS, type="recursive")
  # LL
  gLLbyFun.KFAS    <- oEstimated.KFAS$logLik
  gLLbyErr.KFAS    <- sub.LogLik(oEstimated.KFAS$v[,1], t(oEstimated.KFAS$F), nStateVar)
  # report
  cat(sprintf("hat(sigma^2_epsilon) = %f\n", drop(oFitted.KFAS$model$H)))
  cat(sprintf("hat(mu_1)            = %f\n", coef(oEstimated.KFAS)[1,2]))
  cat(sprintf("hat(beta_1)          = %f\n", coef(oEstimated.KFAS)[1,1]))
  cat("\n") 
  
  cat("max loglikelihood and AIC:\n")
  sub.ShowLogLik (
    n, nStateVar, nHyperParam, 
    gLLbyFun.DLM, gLLbyErr.DLM, gLLbyFun.KFAS, gLLbyErr.KFAS, gCKLogLik*n
  )
  cat("\n") 
}
# - - - - - - - - - - - - - - - - - 
exercise.Chapter5.1b <- function(){
  
  library(dlm)
  library(KFAS)
  
  # constants - - - - - -
  sFILENAME    <- "UKdriversKSI.txt"
  sXFILENAME   <- "logUKpetrolprice.txt"
  gCKLogLik    <- 0.4457201
  # - - - - - - - - - - - 
  
  # reading data
  dfData <- read.table(sFILENAME, skip=1)
  tsData <- ts(log(dfData[,1]), frequency=12, start=c(1969, 1))
  n      <- length(tsData)
  
  dfXData <- read.table(sXFILENAME, skip=1)
  x <- ts(dfXData[,1], frequency=12, start=c(1969, 1))
  stopifnot (length(x) == n)
  
  # model specs
  nStateVar   <- 2
  nHyperParam <- 1
  
  # init value for variances (epsilon)
  gInit     <- var(tsData)
  
  cat("\n")
  
  cat("estimated by dlm package:\n")
  funModel <- function(parm){
    m1 <- dlmModPoly(order=1, dV=exp(parm[1]), dW=0)
    m2 <- dlmModReg(x, dV=0, addInt=FALSE)
    return( m1 + m2 )
  }
  # Fitting
  oMLE.DLM <- dlmMLE(
    tsData, 
    parm  = log(gInit), 
    build = funModel
  )
  stopifnot(oMLE.DLM$convergence == 0)
  oFitted.DLM    <- funModel(oMLE.DLM$par)
  # Filtering
  oFiltered.DLM <- dlmFilter(tsData, oFitted.DLM)
  agStdPredErr.DLM  <- residuals(oFiltered.DLM, type = c("standardized"), sd=FALSE)
  lPredErr.DLM      <- residuals(oFiltered.DLM, type = c("raw"), sd=TRUE)  
  # Smoothing
  oSmoothed.DLM <- dlmSmooth(tsData, oFitted.DLM)
  tsSmoState.DLM   <- dropFirst(oSmoothed.DLM$s)
  # LL
  gLLbyFun.DLM    <- -oMLE.DLM$value - 0.5*n*log(2*pi)
  gLLbyErr.DLM    <- sub.LogLik(lPredErr.DLM$res, lPredErr.DLM$sd^2, nStateVar)
  # report
  cat(sprintf("hat(sigma^2_epsilon) = %f\n", drop(oFitted.DLM$V)))
  cat(sprintf("hat(mu_1)            = %f\n", oSmoothed.DLM$s[2,1]))
  cat(sprintf("hat(beta_1)          = %f\n", oSmoothed.DLM$s[2,2]))
  cat("\n") 
  
  cat("estimated by KFAS package:\n")
  oModel.KFAS <- SSModel(
    tsData ~ SSMtrend(degree=1, Q=list(matrix(0))) + 
             SSMregression(~ -1 + x, Q=matrix(0)), 
    H=matrix(NA)
  )
  ## print(oModel.KFAS) -> states are (x, level)
  ## print(oModel.KFAS$Q) -> Q is (2,2) matrix of 0
  # Fitting
  oFitted.KFAS   <- fitSSM(
    oModel.KFAS, 
    inits  = log(gInit), 
    method = "BFGS"
  )
  stopifnot(oFitted.KFAS$optim.out$convergence == 0)
  # Filtering & Smoothing
  oEstimated.KFAS <- KFS(oFitted.KFAS$model)
  agStdPredErr.KFAS  <- rstandard(oEstimated.KFAS, type="recursive")
  # LL
  gLLbyFun.KFAS    <- oEstimated.KFAS$logLik
  gLLbyErr.KFAS    <- sub.LogLik(oEstimated.KFAS$v[,1], t(oEstimated.KFAS$F), nStateVar)
  # report
  cat(sprintf("hat(sigma^2_epsilon) = %f\n", drop(oFitted.KFAS$model$H)))
  cat(sprintf("hat(mu_1)            = %f\n", coef(oEstimated.KFAS)[1,2]))
  cat(sprintf("hat(beta_1)          = %f\n", coef(oEstimated.KFAS)[1,1]))
  cat("\n") 
  
  cat("max loglikelihood and AIC:\n")
  sub.ShowLogLik (
    n, nStateVar, nHyperParam, 
    gLLbyFun.DLM, gLLbyErr.DLM, gLLbyFun.KFAS, gLLbyErr.KFAS, gCKLogLik*n
  )
  cat("\n") 

  # plot
  sub.plotFigure_twoseries(
    tsData, 
    tsSmoState.DLM[,1] + tsSmoState.DLM[,2] * x, 
    "Figure 5.1",
    c("log UK driver KSI", "deterministic level + beta*log(PETROL PRICE)")
  )
  sub.plotFigure_xyscatterreg(
    tsData, 
    x, 
    tsData ~ x,
    "Figure 5.2", 
    c(
      "log UK driver KSI againsg log PETROL PRICE", 
      "deterministic level + beta*log(PETROL PRICE)"
    )
  )
  sub.plotFigure_residual(
    tsData - (tsSmoState.DLM[,1] + tsSmoState.DLM[,2] * x), 
    "Figure 5.3",
    "irregular"
  )
}
# - - - - - - - - - - - - - - - - - 
exercise.Chapter5.2 <- function(){
  
  library(dlm)
  library(KFAS)
  
  # constants - - - - - -
  sFILENAME    <- "UKdriversKSI.txt"
  sXFILENAME   <- "logUKpetrolprice.txt"
  gCKLogLik    <- 0.6456361
  # - - - - - - - - - - - 
  
  # reading data
  dfData <- read.table(sFILENAME, skip=1)
  tsData <- ts(log(dfData[,1]), frequency=12, start=c(1969, 1))
  n      <- length(tsData)
  
  dfXData <- read.table(sXFILENAME, skip=1)
  x <- ts(dfXData[,1], frequency=12, start=c(1969, 1))
  stopifnot (length(x) == n)
  
  # model specs
  nStateVar   <- 2
  nHyperParam <- 2
  
  # init value for variances (epsilon, xi)
  gInit     <- c(var(tsData), 0.001)  
  
  cat("\n")
  
  cat("estimated by dlm package:\n")
  funModel <- function(parm){
    m1 <- dlmModPoly(order = 1, dV=exp(parm[1]), dW=exp(parm[2]))
    m2 <- dlmModReg(x, dV=0, addInt=FALSE)
    return( m1 + m2 )
  }
  # Fitting
  oMLE.DLM <- dlmMLE(
    tsData, 
    parm  = log(gInit), 
    build = funModel
  )
  stopifnot(oMLE.DLM$convergence == 0)
  oFitted.DLM    <- funModel(oMLE.DLM$par)
  # Filtering
  oFiltered.DLM <- dlmFilter(tsData, oFitted.DLM)
  agStdPredErr.DLM  <- residuals(oFiltered.DLM, type = c("standardized"), sd=FALSE)
  lPredErr.DLM      <- residuals(oFiltered.DLM, type = c("raw"), sd=TRUE)  
  # Smoothing
  oSmoothed.DLM <- dlmSmooth(tsData, oFitted.DLM)
  tsSmoState.DLM   <- dropFirst(oSmoothed.DLM$s)
  # LL
  gLLbyFun.DLM    <- -oMLE.DLM$value - 0.5*n*log(2*pi)
  gLLbyErr.DLM    <- sub.LogLik(lPredErr.DLM$res, lPredErr.DLM$sd^2, nStateVar)
  # report
  agW           <- diag(oFitted.DLM$W)
  cat(sprintf("hat(sigma^2_epsilon) = %f\n", drop(oFitted.DLM$V)))
  cat(sprintf("hat(sigma^2_xi)      = %f\n", agW[1]))
  cat(sprintf("hat(mu_1)            = %f\n", oSmoothed.DLM$s[2,1]))
  cat(sprintf("hat(beta_1)          = %f\n", oSmoothed.DLM$s[2,2]))
  cat("\n") 
  
  cat("estimated by KFAS package:\n")
  oModel.KFAS <- SSModel(
    tsData ~ SSMtrend(degree=1, Q=list(matrix(NA))) + 
      SSMregression(~ -1 + x, Q=matrix(0)), 
    H=matrix(NA)
  )
  ## print(oModel.KFAS)   -> states are (x, level)
  ## print(oModel.KFAS$Q) -> Q is (2,2) matrix where diag is (0, NA)
  # Fitting
  oFitted.KFAS   <- fitSSM(
    oModel.KFAS, 
    inits  = log(gInit)
  )
  stopifnot(oFitted.KFAS$optim.out$convergence == 0)
  # Filtering & Smoothing
  oEstimated.KFAS <- KFS(oFitted.KFAS$model)
  agStdPredErr.KFAS  <- rstandard(oEstimated.KFAS, type="recursive")
  # LL
  gLLbyFun.KFAS    <- oEstimated.KFAS$logLik
  gLLbyErr.KFAS    <- sub.LogLik(oEstimated.KFAS$v[,1], t(oEstimated.KFAS$F), nStateVar)
  # report
  agQ            <- diag(oFitted.KFAS$model$Q[,,1])
  cat(sprintf("hat(sigma^2_epsilon) = %f\n", drop(oFitted.KFAS$model$H)))
  cat(sprintf("hat(sigma^2_xi)      = %f\n", agQ[2]))
  cat(sprintf("hat(mu_1)            = %f\n", coef(oEstimated.KFAS)[1,2]))
  cat(sprintf("hat(beta_1)          = %f\n", coef(oEstimated.KFAS)[1,1]))
  cat("\n") 
  
  cat("max loglikelihood and AIC:\n")
  sub.ShowLogLik (
    n, nStateVar, nHyperParam, 
    gLLbyFun.DLM, gLLbyErr.DLM, gLLbyFun.KFAS, gLLbyErr.KFAS, gCKLogLik*n
  )
  cat("\n") 
  
  # plot
  sub.plotFigure_twoseries(
    tsData, 
    tsSmoState.DLM[,1] + tsSmoState.DLM[,2] * x, 
    "Figure 5.4",
    c("log UK driver KSI", "stochastic level + beta*log(PETROL PRICE)")
  )
  sub.plotFigure_residual(
    tsData - (tsSmoState.DLM[,1] + tsSmoState.DLM[,2] * x), 
    "Figure 5.5",
    "irregular"
  )
}
