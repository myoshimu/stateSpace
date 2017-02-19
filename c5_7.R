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
# - - - - - - - - - - - - - - - - - 
exercise.Chapter6.1 <- function(){
  
  library(dlm)
  library(KFAS)
  
  # constants - - - - - -
  sFILENAME    <- "UKdriversKSI.txt"
  gCKLogLik    <- 0.4573681
  # - - - - - - - - - - - 
  
  # reading data
  dfData <- read.table(sFILENAME, skip=1)
  tsData <- ts(log(dfData[,1]), frequency=12, start=c(1969, 1))
  n      <- length(tsData)
  
  x <- ts(rep(0, n), frequency=12, start=c(1969, 1))
  x[170:n] <- 1
  
  # model specs
  nStateVar   <- 2
  nHyperParam <- 1
  
  # init value for variances (epsilon)
  gInit <- var(tsData)
  
  cat("\n")
  
  cat("estimated by dlm package:\n")
  funModel <- function(parm){
    m1 <- dlmModPoly(order = 1, dV=exp(parm), dW=0)
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
  tsSmoState.DLM   <- dropFirst(oSmoothed.DLM$s)
  # LL
  gLLbyFun.DLM    <- -oMLE.DLM$value - 0.5*n*log(2*pi)
  gLLbyErr.DLM    <- sub.LogLik(lPredErr.DLM$res, lPredErr.DLM$sd^2, nStateVar)
  # report
  cat(sprintf("hat(sigma^2_epsilon) = %f\n", drop(oFitted.DLM$V)))
  cat(sprintf("hat(mu_1)            = %f\n", oSmoothed.DLM$s[2,1]))
  cat(sprintf("hat(lambda_1)        = %f\n", oSmoothed.DLM$s[2,2]))
  cat("\n") 
  
  cat("estimated by KFAS package:\n")
  oModel.KFAS <- SSModel(
    tsData ~ SSMtrend(degree=1, Q=list(matrix(0))) + 
             SSMregression(~ -1 + x, Q=matrix(0)), 
    H=matrix(NA)
  )
  ## print(oModel.KFAS)   -> states are (x, level)
  ## print(oModel.KFAS$Q) -> Q is (2,2) matrix of 0
  # Fitting
  oFitted.KFAS   <- fitSSM(
    oModel.KFAS, 
    inits  = log(gInit), 
    method ="BFGS"
  )
  stopifnot(oFitted.KFAS$optim.out$convergence == 0)
  # Filtering & Smoothing
  oEstimated.KFAS <- KFS(oFitted.KFAS$model, simplify=FALSE)
  agStdPredErr.KFAS  <- rstandard(oEstimated.KFAS, type="recursive")
  # LL
  gLLbyFun.KFAS    <- oEstimated.KFAS$logLik
  gLLbyErr.KFAS    <- sub.LogLik(oEstimated.KFAS$v[,1], t(oEstimated.KFAS$F), nStateVar)
  # report
  cat(sprintf("hat(sigma^2_epsilon) = %f\n", drop(oFitted.KFAS$model$H)))
  cat(sprintf("hat(mu_1)            = %f\n", coef(oEstimated.KFAS)[1,2]))
  cat(sprintf("hat(lambda_1)        = %f\n", coef(oEstimated.KFAS)[1,1]))
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
    "Figure 6.1",
    c("log UK driver KSI", "deterministic level + lambda*(SEATBELT LAW)")
  )
  sub.plotFigure_xyscatterreg(
    tsData, 
    x, 
    tsData ~ x,
    "Figure 6.2", 
    c(
      "log UK driver KSI against SEATBELT LAW", 
      "deterministic level + lambda*(SEATBELT LAW)"
    )
  )
  sub.plotFigure_residual(
    tsData - (tsSmoState.DLM[,1] + tsSmoState.DLM[,2] * x), 
    "Figure 6.3",
    "irregular"
  )
  
  # comparison of prediction error variances between ldm and KFAS
  cat("making additional chart ... \n")
  agPEV.DLM <- as.numeric(lPredErr.DLM$sd^2)
  agPEV.DLM[agPEV.DLM > 1] <- NA
  agPEV.KFAS <- as.numeric(t(oEstimated.KFAS$F))
  tsPEV <- ts(cbind(agPEV.DLM, agPEV.KFAS), frequency=12, start=c(1969, 1))
  colnames(tsPEV) <- c("ldm", "KFAS")
  plot(tsPEV, main="prediction error variance")
  
}
# - - - - - - - - - - - - - - - - - 
exercise.Chapter6.2 <- function(){
  
  library(dlm)
  library(KFAS)
  
  # constants - - - - - -
  sFILENAME    <- "UKdriversKSI.txt"
  gCKLogLik    <- 0.6630851
  # - - - - - - - - - - - 
  
  # reading data
  dfData <- read.table(sFILENAME, skip=1)
  tsData <- ts(log(dfData[,1]), frequency=12, start=c(1969, 1))
  n      <- length(tsData)
  
  x <- ts(rep(0, n), frequency=12, start=c(1969, 1))
  x[170:n] <- 1
  
  # model specs
  nStateVar   <- 2
  nHyperParam <- 2
  
  # init value for variances (epsilon, xi)
  agInit     <- c(var(tsData), 0.001)  
  
  cat("\n")
  
  cat("estimated by dlm package:\n")
  funModel <- function(parm){
    m1 <- dlmModPoly(order = 1, dV=exp(parm[1]), dW=exp(parm[2]))
    m2 <- dlmModReg(x, dV=0, addInt = FALSE)
    return( m1 + m2 )
  }
  # Fitting
  oMLE.DLM <- dlmMLE(
    tsData, 
    parm  = log(agInit), 
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
  cat(sprintf("hat(lambda_1)        = %f\n", oSmoothed.DLM$s[2,2]))
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
    inits = log(agInit)
  )
  stopifnot(oFitted.KFAS$optim.out$convergence == 0)
  # Filtering & Smoothing
  oEstimated.KFAS <- KFS(oFitted.KFAS$model, simplify=FALSE)
  agStdPredErr.KFAS  <- rstandard(oEstimated.KFAS, type="recursive")
  # LL
  gLLbyFun.KFAS    <- oEstimated.KFAS$logLik
  gLLbyErr.KFAS    <- sub.LogLik(oEstimated.KFAS$v[,1], t(oEstimated.KFAS$F), nStateVar)
  # report
  agQ            <- diag(oFitted.KFAS$model$Q[,,1])
  cat(sprintf("hat(sigma^2_epsilon) = %f\n", drop(oFitted.KFAS$model$H)))
  cat(sprintf("hat(sigma^2_xi)      = %f\n", agQ[2]))
  cat(sprintf("hat(mu_1)            = %f\n", coef(oEstimated.KFAS)[1,2]))
  cat(sprintf("hat(lambda_1)        = %f\n", coef(oEstimated.KFAS)[1,1]))
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
    "Figure 6.4",
    c("log UK driver KSI", "stochastic level + lambda*(SEATBELT LAW)")
  )
  sub.plotFigure_residual(
    tsData - (tsSmoState.DLM[,1] + tsSmoState.DLM[,2] * x), 
    "Figure 6.5",
    "irregular"
  )
}
# - - - - - - - - - - - - - - - - - 
exercise.Chapter7.1 <- function(){
  
  library(dlm)
  library(KFAS)
  
  # constants - - - - - -
  sFILENAME    <- "UKdriversKSI.txt"
  sXFILENAME   <- "logUKpetrolprice.txt"
  gCKLogLik    <- 0.8023777
  # - - - - - - - - - - - 
  
  # reading data
  dfData <- read.table(sFILENAME, skip=1)
  tsData <- ts(log(dfData[,1]), frequency=12, start=c(1969, 1))
  n      <- length(tsData)
  
  dfX1Data <- read.table(sXFILENAME, skip=1)
  x2 <- rep(0, n)
  x2[170:n] <- 1
  X <- ts(
    cbind(dfX1Data[,1], x2),
    frequency=12, start=c(1969, 1)
  )
  colnames(X) <- c("x1", "x2")
  
  # model specs
  nStateVar   <- 14
  nHyperParam <- 1
  
  # init value for variances (epsilon)
  gInit  <- var(tsData)
  
  cat("\n")
  
  cat("estimated by dlm package:\n")
  funModel <- function(parm){
    m1 <- dlmModPoly(order = 1, dV=exp(parm), dW=0)
    m2 <- dlmModSeas(12, dV=0, dW=rep(0,11))
    m3 <- dlmModReg(X, dV=0, addInt=FALSE)
    return( m1 + m2 + m3)
  }
  ## print(funModel(log(123)))
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
  tsSmoState.DLM  <- dropFirst(oSmoothed.DLM$s)
  # LL
  gLLbyFun.DLM    <- -oMLE.DLM$value - 0.5*n*log(2*pi)
  gLLbyErr.DLM    <- sub.LogLik(lPredErr.DLM$res, lPredErr.DLM$sd^2, nStateVar)
  # report
  cat(sprintf("hat(sigma^2_epsilon) = %f\n", drop(oFitted.DLM$V)))
  cat(sprintf("hat(mu_1)            = %f\n", oSmoothed.DLM$s[2,1]))
  cat(sprintf("hat(gamma_1)         = %f\n", oSmoothed.DLM$s[2,2]))
  cat(sprintf("hat(beta_1)          = %f\n", oSmoothed.DLM$s[2,13]))
  cat(sprintf("hat(lambda_1)        = %f\n", oSmoothed.DLM$s[2,14]))
  cat("\n") 
  
  cat("estimated by KFAS package:\n")
  oModel.KFAS <- SSModel(
    tsData ~ SSMtrend(degree=1, Q=list(matrix(0))) + 
             SSMseasonal(12, sea.type="dummy") + 
             SSMregression(~ -1 + X, Q=matrix(0)), 
    H=matrix(NA)
  )
  ## print(oModel.KFAS)   -> states are (x1, x2, level, seasondummies)
  ## print(oModel.KFAS$Q) -> Q is (2,2) matrix where diag is (0, 0)
  # fitting
  oFitted.KFAS   <- fitSSM(
    oModel.KFAS, 
    inits  = log(gInit), 
    method = "BFGS"
  )
  stopifnot(oFitted.KFAS$optim.out$convergence == 0)
  # Filtering & Smoothing
  oEstimated.KFAS <- KFS(oFitted.KFAS$model, simplify=FALSE)
  agStdPredErr.KFAS  <- rstandard(oEstimated.KFAS, type="recursive")
  # LL
  gLLbyFun.KFAS    <- oEstimated.KFAS$logLik
  gLLbyErr.KFAS    <- sub.LogLik(oEstimated.KFAS$v[,1], t(oEstimated.KFAS$F), nStateVar)
  # report
  cat(sprintf("hat(sigma^2_epsilon) = %f\n", drop(oFitted.KFAS$model$H)))
  cat(sprintf("hat(mu_1)            = %f\n", coef(oEstimated.KFAS)[1,3]))
  cat(sprintf("hat(gamma_1)         = %f\n", coef(oEstimated.KFAS)[1,4]))
  cat(sprintf("hat(beta_1)          = %f\n", coef(oEstimated.KFAS)[1,1]))
  cat(sprintf("hat(lambda_1)        = %f\n", coef(oEstimated.KFAS)[1,2]))
  cat("\n") 
  
  cat("max loglikelihood and AIC:\n")
  sub.ShowLogLik (
    n, nStateVar, nHyperParam, 
    gLLbyFun.DLM, gLLbyErr.DLM, gLLbyFun.KFAS, gLLbyErr.KFAS, gCKLogLik*n
  )
  cat("\n") 
  
  cat("diagnostics of standardized prediction error (estimated by dlm):\n")
  sub.ShowDiagnostics (
    agStdPredErr.DLM, nStateVar, nHyperParam, 15, c(1,12)
  )
  
  # plot
  sub.plotFigure_twoseries(
    tsData, 
    tsSmoState.DLM[,1] + (tsSmoState.DLM[,13] * X[,1]) + (tsSmoState.DLM[,14] * X[,2]), 
    "Figure 7.1",
    c("log UK driver KSI", "det.level+beta*log(PETROL PRICE)+lambda*(SEATBELT LAW)")
  )
  
  # plot for section 7.3
  sub.plotFigure_ACF(
    tsData - (tsSmoState.DLM[,1] + tsSmoState.DLM[,2] + (tsSmoState.DLM[,13] * X[,1]) + (tsSmoState.DLM[,14] * X[,2])), 
    15, 
    "Figure 7.5", 
    "ACF-det.level and seasonal model residuals"
  )
}
# - - - - - - - - - - - - - - - - - 
exercise.Chapter7.2 <- function(){
  
  library(dlm)
  library(KFAS)
  
  # constants - - - - - -
  sFILENAME    <- "UKdriversKSI.txt"
  sXFILENAME   <- "logUKpetrolprice.txt"
  gCKLogLik    <- 0.9825225
  # - - - - - - - - - - - 
  
  # reading data
  dfData <- read.table(sFILENAME, skip=1)
  tsData <- ts(log(dfData[,1]), frequency=12, start=c(1969, 1))
  n      <- length(tsData)
  
  dfX1Data <- read.table(sXFILENAME, skip=1)
  x2 <- rep(0, n)
  x2[170:n] <- 1
  X <- ts(
    cbind(dfX1Data[,1], x2),
    frequency=12, start=c(1969, 1)
  )
  colnames(X) <- c("x1", "x2")
  
  # model specs
  nStateVar   <- 14
  nHyperParam <- 3
  
  # init value for variances (epsilon, xi, omega)
  agInit     <- c(var(tsData), 0.0001, 0.0001)
  
  cat("\n")
  
  cat("estimated by dlm package:\n")
  funModel <- function(parm){
    m1 <- dlmModPoly(order = 1, dV=exp(parm[1]), dW=exp(parm[2]))
    m2 <- dlmModSeas(12, dV=0, dW=c(exp(parm[3]), rep(0,10)))
    m3 <- dlmModReg(X, dV=0, addInt=FALSE)
    return( m1 + m2 + m3 )
  }
  ## print(funModel(log(c(123, 456, 789))))
  # Fitting
  oMLE.DLM <- dlmMLE(
    tsData, 
    parm  = log(agInit), 
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
  cat(sprintf("hat(sigma^2_xi)      = %e\n", agW[1]))
  cat(sprintf("hat(sigma^2_omega)   = %e\n", agW[2]))
  cat(sprintf("hat(mu_1)            = %f\n", oSmoothed.DLM$s[2,1]))
  cat(sprintf("hat(gamma_1)         = %f\n", oSmoothed.DLM$s[2,2]))
  cat(sprintf("hat(beta_1)          = %f\n", oSmoothed.DLM$s[1,13]))
  cat(sprintf("hat(lambda_1)        = %f\n", oSmoothed.DLM$s[1,14]))
  cat("\n") 
  
  cat("estimated by KFAS package:\n")
  oModel.KFAS <- SSModel(
    tsData ~ SSMtrend(degree=1, Q=list(matrix(NA))) + 
             SSMseasonal(12, sea.type="dummy", Q=matrix(NA), n=n) + 
             SSMregression(~ -1 + X, Q=matrix(0, ncol=2, nrow=2)), 
    H=matrix(NA)
  )
  ## print(oModel.KFAS)   -> states are (x1, x2, level, seasondummies)
  ## print(oModel.KFAS$Q) -> Q is (4,4) matrix where diag is (0,0,NA,NA)
  # Fitting
  oFitted.KFAS   <- fitSSM(
    oModel.KFAS, 
    inits  = log(agInit)
  )
  stopifnot(oFitted.KFAS$optim.out$convergence == 0)
  # Filtering & Smoothing
  oEstimated.KFAS <- KFS(oFitted.KFAS$model, simplify=FALSE)
  agStdPredErr.KFAS  <- rstandard(oEstimated.KFAS, type="recursive")
  # LL
  gLLbyFun.KFAS    <- oEstimated.KFAS$logLik
  gLLbyErr.KFAS    <- sub.LogLik(oEstimated.KFAS$v[,1], t(oEstimated.KFAS$F), nStateVar)
  # report
  agQ            <- diag(oFitted.KFAS$model$Q[,,1])
  cat(sprintf("hat(sigma^2_epsilon) = %f\n", drop(oFitted.KFAS$model$H)))
  cat(sprintf("hat(sigma^2_xi)      = %e\n", agQ[3]))
  cat(sprintf("hat(sigma^2_omega)   = %e\n", agQ[4]))
  cat(sprintf("hat(mu_1)            = %f\n", coef(oEstimated.KFAS)[1,3]))
  cat(sprintf("hat(gamma_1)         = %f\n", coef(oEstimated.KFAS)[1,4]))
  cat(sprintf("hat(beta_1)          = %f\n", coef(oEstimated.KFAS)[1,1]))
  cat(sprintf("hat(lambda_1)        = %f\n", coef(oEstimated.KFAS)[1,2]))
  cat("\n") 
  
  cat("max loglikelihood and AIC:\n")
  sub.ShowLogLik (
    n, nStateVar, nHyperParam, 
    gLLbyFun.DLM, gLLbyErr.DLM, gLLbyFun.KFAS, gLLbyErr.KFAS, gCKLogLik*n
  )
  cat("\n") 
  
  cat("diagnostics of standardized prediction error (estimated by dlm):\n")
  sub.ShowDiagnostics (
    agStdPredErr.DLM, nStateVar, nHyperParam, 15, c(1,12)
  )
  
  # plot
  sub.plotFigure_twoseries(
    tsData, 
    tsSmoState.DLM[,1] + (tsSmoState.DLM[,13] * X[,1]) + (tsSmoState.DLM[,14] * X[,2]), 
    "Figure 7.2",
    c("log UK driver KSI", "sto.level+beta*log(PETROL PRICE)+lambda*(SEATBELT LAW)")
  )
  sub.plotFigure_twoseries(
    NULL, 
    tsSmoState.DLM[,2], 
    "Figure 7.3",
    c("stochastic seasonal")
  )
  sub.plotFigure_residual(
    tsData - (tsSmoState.DLM[,1] + tsSmoState.DLM[,2] + (tsSmoState.DLM[,13] * X[,1]) + (tsSmoState.DLM[,14] * X[,2])), 
    "Figure 7.4",
    "irregular"
  )
}
# - - - - - - - - - - - - - - - - - 
exercise.Chapter7.3 <- function(){
  
  library(dlm)
  library(KFAS)
  
  # constants - - - - - -
  sFILENAME    <- "UKdriversKSI.txt"
  sXFILENAME   <- "logUKpetrolprice.txt"
  gCKLogLik    <- 0.9798650
  # - - - - - - - - - - - 
  
  # reading data
  dfData <- read.table(sFILENAME, skip=1)
  tsData <- ts(log(dfData[,1]), frequency=12, start=c(1969, 1))
  n      <- length(tsData)
  
  dfX1Data <- read.table(sXFILENAME, skip=1)
  x2 <- rep(0, n)
  x2[170:n] <- 1
  X <- ts(
    cbind(dfX1Data[,1], x2),
    frequency=12, start=c(1969, 1)
  )
  colnames(X) <- c("x1", "x2")
  
  # model specs
  nStateVar   <- 14
  nHyperParam <- 2
  
  # init value for variances (epsilon, xi)
  agInit     <- c(var(tsData), 0.001)
  
  cat("\n")
  
  cat("estimated by dlm package:\n")
  funModel <- function(parm){
    m1 <- dlmModPoly(order = 1, dV=exp(parm[1]), dW=exp(parm[2]))
    m2 <- dlmModSeas(12, dV=0, dW=rep(0,11))
    m3 <- dlmModReg(X, dV=0, addInt=FALSE)
    return( m1 + m2 + m3 )
  }
  ## print(funModel(log(c(123, 456, 789))))
  # Fitting
  oMLE.DLM <- dlmMLE(
    tsData, 
    parm  = log(agInit), 
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
  cat(sprintf("hat(sigma^2_xi)      = %e\n", agW[1]))
  cat(sprintf("hat(mu_1)            = %f\n", oSmoothed.DLM$s[2,1]))
  cat(sprintf("hat(gamma_1)         = %f\n", oSmoothed.DLM$s[2,2]))
  cat(sprintf("hat(beta_1)          = %f\n", oSmoothed.DLM$s[2,13]))
  cat(sprintf("hat(lambda_1)        = %f\n", oSmoothed.DLM$s[2,14]))
  cat("\n") 
  
  cat("estimated by KFAS package:\n")
  oModel.KFAS <- SSModel(
    tsData ~ SSMtrend(degree=1, Q=list(matrix(NA))) + 
             SSMseasonal(12, sea.type="dummy") + 
             SSMregression(~ -1 + X, Q=matrix(0, ncol=2, nrow=2)), 
    H=matrix(NA)
  )
  ## print(oModel.KFAS)   -> states are (x1, x2, level, seasondummies)
  ## print(oModel.KFAS$Q) -> Q is (3,3) matrix where diag is (0,0,NA)
  # Fitting
  oFitted.KFAS   <- fitSSM(
    oModel.KFAS, 
    inits  = log(agInit)
  )
  stopifnot(oFitted.KFAS$optim.out$convergence == 0)
  # Filtering & Smoothing
  # add disturbance smoothing to plot later
  oEstimated.KFAS <- KFS(oFitted.KFAS$model, smoothing=c("state", "mean", "disturbance"))
  agStdPredErr.KFAS  <- rstandard(oEstimated.KFAS, type="recursive")
  # LL
  gLLbyFun.KFAS    <- oEstimated.KFAS$logLik
  gLLbyErr.KFAS    <- sub.LogLik(oEstimated.KFAS$v[,1], t(oEstimated.KFAS$F), nStateVar)
  # report
  agQ            <- diag(oFitted.KFAS$model$Q[,,1])
  cat(sprintf("hat(sigma^2_epsilon) = %f\n", drop(oFitted.KFAS$model$H)))
  cat(sprintf("hat(sigma^2_xi)      = %e\n", agQ[3]))
  cat(sprintf("hat(mu_1)            = %f\n", coef(oEstimated.KFAS)[1,3]))
  cat(sprintf("hat(gamma_1)         = %f\n", coef(oEstimated.KFAS)[1,4]))
  cat(sprintf("hat(beta_1)          = %f\n", coef(oEstimated.KFAS)[1,1]))
  cat(sprintf("hat(lambda_1)        = %f\n", coef(oEstimated.KFAS)[1,2]))
  cat("\n") 
  
  cat("max loglikelihood and AIC:\n")
  sub.ShowLogLik (
    n, nStateVar, nHyperParam, 
    gLLbyFun.DLM, gLLbyErr.DLM, gLLbyFun.KFAS, gLLbyErr.KFAS, gCKLogLik*n
  )
  cat("\n") 
  
  # plot
  sub.plotFigure_ACF(
    tsData - (tsSmoState.DLM[,1] + tsSmoState.DLM[,2] + (tsSmoState.DLM[,13] * X[,1]) + (tsSmoState.DLM[,14] * X[,2])), 
    15, 
    "Figure 7.6", 
    "ACF-sto.level and det.seasonal model residuals"
  )
  
  # plot for section 8.5
  sub.plotFigure_twoseries(
    NULL, 
    window(agStdPredErr.DLM, start=c(1970,3)),
    "Figure 8.8",
    "standardized one-step prediction errors"
  )
  sub.plotFigure_ACF(
    window(agStdPredErr.DLM, start=c(1970,3)), 
    10, 
    "Figure 8.9", 
    "ACF-standardized one-step prediction errors"
  )
  cat("making Figure 8.10 ... \n")
  hist(
    agStdPredErr.DLM[-(1:nStateVar)], 
    breaks=17,
    main = "Figure 8.10"
  )
  
  sub.plotFigure_auxresidual(
    rstandard(oEstimated.KFAS, type="state"), 
    "Figure 8.12a", 
    "Structural level break t-tests"
  )
  sub.plotFigure_auxresidual(
    rstandard(oEstimated.KFAS, type="pearson"), 
    "Figure 8.12b", 
    "Outlier t-tests"
  )
}
# - - - - - - - - - - - - - - - - - 
exercise.Chapter7.4 <- function(){
  
  library(dlm)
  library(KFAS)
  
  # constants - - - - - -
  sFILENAME <- "UKinflation.txt"
  gCKLogLik <- 3.305023
  # - - - - - - - - - - - 
  
  # read data
  dfData <- read.table(sFILENAME, skip=1)
  tsData <- ts(dfData[,1], frequency=4, start=c(1950, 1))
  n      <- length(tsData)
  
  x1 <- ts(rep(0, n), frequency=4, start=c(1950, 1))
  x2 <- ts(rep(0, n), frequency=4, start=c(1950, 1))
  x1[which(time(x1) == 1975.25)] <- 1
  x2[which(time(x2) == 1979.50)] <- 1
  X <- cbind(x1, x2)
  
  # model specs
  nStateVar   <- 6
  nHyperParam <- 3
  
  # init value for variances (epsilon, xi, gamma)
  agInit     <- c(var(tsData), 0.00001, 0.00001)  
  
  cat("\n")
  
  cat("estimated by dlm package:\n")
  funModel <- function(parm){
    m1 <- dlmModPoly(order = 1, dV=exp(parm[1]), dW=exp(parm[2]))
    m2 <- dlmModSeas(4, dV=0, dW=c(exp(parm[3]), rep(0,2)))
    m3 <- dlmModReg(X, dV=0, addInt=FALSE)
    return( m1 + m2 + m3 )
  }
  # print(funModel(log(c(123,456,789))))
  oMLE.DLM <- dlmMLE(
    tsData, 
    parm  = log(agInit), 
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
  cat(sprintf("hat(sigma^2_epsilon) = %e\n", drop(oFitted.DLM$V)))
  cat(sprintf("hat(sigma^2_xi)      = %e\n", agW[1]))
  cat(sprintf("hat(sigma^2_gamma)   = %e\n", agW[2]))
  cat(sprintf("hat(mu_1)            = %f\n", oSmoothed.DLM$s[2,1]))
  cat(sprintf("hat(gamma_1)         = %f\n", oSmoothed.DLM$s[2,2]))
  cat(sprintf("hat(lambda1_1)       = %f\n", oSmoothed.DLM$s[2,5]))
  cat(sprintf("hat(lambda2_1)       = %f\n", oSmoothed.DLM$s[2,6]))
  cat("\n") 
  
  cat("estimated by KFAS package:\n")
  oModel.KFAS <- SSModel(
    tsData ~ SSMtrend(degree=1, Q=list(matrix(NA))) + 
             SSMseasonal(4, sea.type="dummy", Q=matrix(NA), n=n) + 
             SSMregression(~ -1 + X, Q=matrix(0, ncol=2, nrow=2)), 
    H=matrix(NA)
  )
  ## print(oModel.KFAS)   -> states are (x1, x2, level, seasondummies)
  ## print(oModel.KFAS$Q) -> Q is (3,3) matrix where diag is (0,0,NA,NA)
  # Fitting
  oFitted.KFAS   <- fitSSM(
    oModel.KFAS, 
    inits  = log(agInit)
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
  cat(sprintf("hat(sigma^2_epsilon) = %e\n", drop(oFitted.KFAS$model$H)))
  cat(sprintf("hat(sigma^2_xi)      = %e\n", agQ[3]))
  cat(sprintf("hat(sigma^2_omega)   = %e\n", agQ[4]))
  cat(sprintf("hat(mu_1)            = %f\n", coef(oEstimated.KFAS)[1,3]))
  cat(sprintf("hat(gamma_1)         = %f\n", coef(oEstimated.KFAS)[1,4]))
  cat(sprintf("hat(lambda1_1)       = %f\n", coef(oEstimated.KFAS)[1,1]))
  cat(sprintf("hat(lambda1_1)       = %f\n", coef(oEstimated.KFAS)[1,2]))
  cat("\n") 
  
  cat("max loglikelihood and AIC:\n")
  sub.ShowLogLik (
    n, nStateVar, nHyperParam, 
    gLLbyFun.DLM, gLLbyErr.DLM, gLLbyFun.KFAS, gLLbyErr.KFAS, gCKLogLik*n
  )
  cat("\n") 
  
  cat("diagnostics of standardized prediction error (estimated by dlm):\n")
  sub.ShowDiagnostics (
    agStdPredErr.DLM, nStateVar, nHyperParam, 10, c(1,4)
  )
  
  # plot
  sub.plotFigure_twoseries(
    tsData, 
    tsSmoState.DLM[,1] + tsSmoState.DLM[,5] * X[,1] + tsSmoState.DLM[,6] * X[,2], 
    "Figure 7.7a",
    c("quarterly price changes in U.K.", "sto.level+pulse intervention variables")
  )
  sub.plotFigure_twoseries(
    NULL, 
    tsSmoState.DLM[,2], 
    "Figure 7.7b",
    "stochastic seasonal"
  )
  sub.plotFigure_residual(
    tsData - (tsSmoState.DLM[,1] + tsSmoState.DLM[,2] + tsSmoState.DLM[,5] * X[,1] + tsSmoState.DLM[,6] * X[,2]), 
    "Figure 7.7c",
    "irregular"
  )
   
}
