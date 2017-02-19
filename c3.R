# - - - - - - - - - - - - - - - - - 
exercise.Chapter3.1 <- function(){
  
  library(dlm)
  library(KFAS)
  
  # constants - - - - - -
  sFILENAME    <- "UKdriversKSI.txt"
  gCKLogLik    <- 0.4140728
  # - - - - - - - - - - - 
  
  # read data
  dfData <- read.table(sFILENAME, skip=1)
  tsData <- ts(log(dfData[,1]), frequency=12, start=c(1969, 1))
  n      <- length(tsData)
  
  # model specs
  nStateVar   <- 2
  nHyperParam <- 1
  
  # init value for sigma^2_epsilon
  gInit     <- var(tsData)
  
  cat("\n")
  
  cat("estimated by dlm package:\n")
  funModel <- function(parm){
    dlmModPoly(order = 2, dV = exp(parm), dW = c(0,0))
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
  cat(sprintf("hat(nu_1)            = %f\n", oSmoothed.DLM$s[2,2]))
  cat("\n") 
  
  cat("estimated by KFAS package:\n")
  oModel.KFAS <- SSModel(
    tsData ~ SSMtrend(degree=2, Q=list(matrix(0), matrix(0))), 
    H=matrix(NA)
  )
  ## print(oModel.KFAS)  -> states are {level, slope}
  # fitting
  oFitted.KFAS   <- fitSSM(
    oModel.KFAS, 
    inits  = log(gInit), 
    method = "BFGS"
  )
  # filtering & smoothing
  oEstimated.KFAS <- KFS(oFitted.KFAS$model)
  agStdPredErr.KFAS  <- oEstimated.KFAS$v[,1] / sqrt(t(oEstimated.KFAS$F))
  # LL
  gLLbyFun.KFAS    <- oEstimated.KFAS$logLik
  gLLbyErr.KFAS    <- sub.LogLik(oEstimated.KFAS$v[,1], t(oEstimated.KFAS$F), nStateVar)
  # report
  cat(sprintf("hat(sigma^2_epsilon) = %f\n", drop(oFitted.KFAS$model$H)))
  cat(sprintf("hat(mu_1)            = %f\n", oEstimated.KFAS$alphahat[1,1]))
  cat(sprintf("hat(nu_1)            = %f\n", oEstimated.KFAS$alphahat[1,2]))
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
  cat("\n") 
  
}
# - - - - - - - - - - - - - - - - - 
exercise.Chapter3.2 <- function(){
  
  library(dlm)
  library(KFAS)
  
  # constants - - - - - -
  sFILENAME <- "UKdriversKSI.txt"
  gCKLogLik <- 0.6247935
  # - - - - - - - - - - - 
  
  # read data
  dfData <- read.table(sFILENAME, skip=1)
  tsData <- ts(log(dfData[,1]), frequency=12, start=c(1969, 1))
  n      <- length(tsData)
  
  # model specs
  nStateVar   <- 2
  nHyperParam <- 3
  
  # init value for variances (epsilon, xi, zeta)
  agInit     <- c(var(tsData),0.001,0.001)  
  
  cat("\n")
  
  cat("estimated by dlm package:\n")
  funModel <- function(parm){
    dlmModPoly(order = 2, dV = exp(parm[1]), dW = exp(parm[2:3]))
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
  oFiltered.DLM    <- dlmFilter(tsData, oFitted.DLM)
  agStdPredErr.DLM <- residuals(oFiltered.DLM, type = c("standardized"), sd=FALSE)
  lPredErr.DLM     <- residuals(oFiltered.DLM, type = c("raw"), sd=TRUE)  
  # Smoothing
  oSmoothed.DLM    <- dlmSmooth(tsData, oFitted.DLM)
  tsSmoState.DLM   <- dropFirst(oSmoothed.DLM$s)
  # LL
  gLLbyFun.DLM    <- -oMLE.DLM$value - 0.5*n*log(2*pi)
  gLLbyErr.DLM    <- sub.LogLik(lPredErr.DLM$res, lPredErr.DLM$sd^2, nStateVar)
  # report
  agW <- diag(oFitted.DLM$W)
  cat(sprintf("hat(sigma^2_epsilon) = %f\n", drop(oFitted.DLM$V)))
  cat(sprintf("hat(sigma^2_xi)      = %f\n", agW[1]))
  cat(sprintf("hat(sigma^2_zeta)    = %f\n", agW[2]))
  cat(sprintf("hat(mu_1)            = %f\n", oSmoothed.DLM$s[2,1]))
  cat(sprintf("hat(nu_1)            = %f\n", oSmoothed.DLM$s[2,2]))
  cat("\n") 
  
  cat("estimated by KFAS package:\n")
  oModel.KFAS <- SSModel(
    tsData ~ SSMtrend(degree=2, Q=list(matrix(NA), matrix(NA))), 
    H=matrix(NA)
  )
  ## print(oModel.KFAS) -> states are {level, slope}
  # Fitting
  oFitted.KFAS   <- fitSSM(
    oModel.KFAS, 
    inits = log(agInit)
  )
  # Filtering & Smoothing
  oEstimated.KFAS   <- KFS(oFitted.KFAS$model)
  agStdPredErr.KFAS <- oEstimated.KFAS$v[,1] / sqrt(t(oEstimated.KFAS$F))
  # LL
  gLLbyFun.KFAS    <- oEstimated.KFAS$logLik
  gLLbyErr.KFAS    <- sub.LogLik(oEstimated.KFAS$v[,1], t(oEstimated.KFAS$F), nStateVar)
  # report
  agQ <- diag(oFitted.KFAS$model$Q[,,1])
  cat(sprintf("hat(sigma^2_epsilon) = %f\n", drop(oFitted.KFAS$model$H)))
  cat(sprintf("hat(sigma^2_xi)      = %f\n", agQ[1]))
  cat(sprintf("hat(sigma^2_zeta)    = %f\n", agQ[2]))
  cat(sprintf("hat(mu_1)            = %f\n", oEstimated.KFAS$alphahat[1,1]))
  cat(sprintf("hat(nu_1)            = %f\n", oEstimated.KFAS$alphahat[1,2]))
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
  cat("\n") 
  
  # plot
  sub.plotFigure_twoseries(
    tsData, 
    tsSmoState.DLM[,1],
    "Figure 3.1",
    c("log UK drivers KSI","stochastic level and slope")
  )
  sub.plotFigure_twoseries(
    NULL, 
    dropFirst(tsSmoState.DLM[,2]),
    "Figure 3.2",
    "stochastic slope"
  )
  sub.plotFigure_residual(
    tsData - tsSmoState.DLM[,1], 
    "Figure 3.3",
    "irregular"
  )
  
}
# - - - - - - - - - - - - - - - - - 
exercise.Chapter3.3 <- function(){
  
  library(dlm)
  library(KFAS)
  
  # constants - - - - - -
  sFILENAME <- "UKdriversKSI.txt"
  gCKLogLik <- 0.6247935
  # - - - - - - - - - - - 
  
  # read data
  dfData <- read.table(sFILENAME, skip=1)
  tsData <- ts(log(dfData[,1]), frequency=12, start=c(1969, 1))
  n      <- length(tsData)
  
  # model specs
  nStateVar   <- 2
  nHyperParam <- 2
  
  # init value for variances (epsilon, xi)
  agInit     <- c(var(tsData),0.001)  
  
  cat("\n")
  
  cat("estimated by dlm package:\n")
  funModel <- function(parm){
    dlmModPoly(order = 2, dV = exp(parm[1]), dW = c(exp(parm[2]), 0))
  }
  # Fitting
  oMLE.DLM <- dlmMLE(
    tsData, 
    parm  = log(agInit), 
    build = funModel
  )
  stopifnot(oMLE.DLM$convergence == 0)
  oFitted.DLM <- funModel(oMLE.DLM$par)
  # Filtering
  oFiltered.DLM    <- dlmFilter(tsData, oFitted.DLM)
  agStdPredErr.DLM <- residuals(oFiltered.DLM, type = c("standardized"), sd=FALSE)
  lPredErr.DLM     <- residuals(oFiltered.DLM, type = c("raw"), sd=TRUE)  
  # Smoothing
  oSmoothed.DLM    <- dlmSmooth(tsData, oFitted.DLM)
  tsSmoState.DLM   <- dropFirst(oSmoothed.DLM$s)
  # LL
  gLLbyFun.DLM     <- -oMLE.DLM$value - 0.5*n*log(2*pi)
  gLLbyErr.DLM     <- sub.LogLik(lPredErr.DLM$res, lPredErr.DLM$sd^2, nStateVar)
  # report
  agW <- diag(oFitted.DLM$W)
  cat(sprintf("hat(sigma^2_epsilon) = %f\n", drop(oFitted.DLM$V)))
  cat(sprintf("hat(sigma^2_xi)      = %f\n", agW[1]))
  cat(sprintf("hat(mu_1)            = %f\n", oSmoothed.DLM$s[2,1]))
  cat(sprintf("hat(nu_1)            = %f\n", oSmoothed.DLM$s[2,2]))
  cat("\n") 
  
  cat("estimated by KFAS package:\n")
  oModel.KFAS <- SSModel(
    tsData ~ SSMtrend(degree=2, Q=list(matrix(NA), matrix(0))), 
    H=matrix(NA)
  )
  ## print(oModel.KFAS) -> states are (level, slope)
  # Fitting
  oFitted.KFAS   <- fitSSM(
    oModel.KFAS, 
    inits = log(agInit)
  )
  # Filtering & Smoothing
  oEstimated.KFAS <- KFS(oFitted.KFAS$model)
  agStdPredErr.KFAS  <- oEstimated.KFAS$v[,1] / sqrt(t(oEstimated.KFAS$F))
  # LL
  gLLbyFun.KFAS    <- oEstimated.KFAS$logLik
  gLLbyErr.KFAS    <- sub.LogLik(oEstimated.KFAS$v[,1], t(oEstimated.KFAS$F), nStateVar)
  # report
  agQ            <- diag(oFitted.KFAS$model$Q[,,1])
  cat(sprintf("hat(sigma^2_epsilon) = %f\n", drop(oFitted.KFAS$model$H)))
  cat(sprintf("hat(sigma^2_xi)      = %f\n", agQ[1]))
  cat(sprintf("hat(mu_1)            = %f\n", oEstimated.KFAS$alphahat[1,1]))
  cat(sprintf("hat(nu_1)            = %f\n", oEstimated.KFAS$alphahat[1,2]))
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
  cat("\n") 
  
  # plot
  sub.plotFigure_twoseries(
    tsData, 
    tsSmoState.DLM[,1], 
    "Figure 3.4",
    c("log UK drivers KSI", "stochastic level and deterministic slope")
  )
 
}
# - - - - - - - - - - - - - - - - - 
exercise.Chapter3.4a <- function(){
  # the first paragraph of section 3.4
  # local linear trend model with stochastic level & slope
  
  library(dlm)
  library(KFAS)
  
  # constants - - - - - -
  sFILENAME <- "NorwayFinland.txt"
  gCKLogLik <- 0.7864746
  # - - - - - - - - - - - 
  
  # read data
  dfData <- read.table(sFILENAME, skip=1)
  anTime <- as.integer(dfData[,1])
  stopifnot(anTime == seq(from=anTime[1], to=anTime[1]+length(anTime)-1))
  tsData <- ts(log(dfData[,3]), start=anTime[1])
  n      <- length(tsData)
  
  # model specs
  nStateVar   <- 2
  nHyperParam <- 3
  
  # init value for variances (epsilon, xi, zeta)
  agInit     <- c(var(tsData),0.001,0.001)  
  
  cat("\n")
  
  cat("estimated by dlm package: \n")
  funModel <- function(parm){
    dlmModPoly(order = 2, dV = exp(parm[1]), dW = exp(parm[2:3]))
  }
  # fitting
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
  # LL
  gLLbyFun.DLM    <- -oMLE.DLM$value - 0.5*n*log(2*pi)
  gLLbyErr.DLM    <- sub.LogLik(lPredErr.DLM$res, lPredErr.DLM$sd^2, nStateVar)
  # report
  agW <- diag(oFitted.DLM$W)
  cat(sprintf("hat(sigma^2_epsilon) = %f\n", drop(oFitted.DLM$V)))
  cat(sprintf("hat(sigma^2_xi)      = %e\n", agW[1]))
  cat(sprintf("hat(sigma^2_zeta)    = %f\n", agW[2]))
  cat(sprintf("hat(mu_1)            = %f\n", oSmoothed.DLM$s[2,1]))
  cat(sprintf("hat(nu_1)            = %f\n", oSmoothed.DLM$s[2,2]))
  cat("\n") 
   
  cat("estimated by KFAS package: \n")
  oModel.KFAS <- SSModel(
    tsData ~ SSMtrend(degree=2, Q=list(matrix(NA), matrix(NA))), 
    H=matrix(NA)
  )
  ## print(oModel.KFAS) -> states are {level, slopt}
  # fitting
  oFitted.KFAS <- fitSSM(
    oModel.KFAS, 
    inits = log(agInit)
  )
  stopifnot(oFitted.KFAS$optim.out$convergence == 0)
  # Filtering & Smoothing
  oEstimated.KFAS <- KFS(oFitted.KFAS$model)
  agStdPredErr.KFAS  <- oEstimated.KFAS$v[,1] / sqrt(t(oEstimated.KFAS$F))
  # LL
  gLLbyFun.KFAS    <- oEstimated.KFAS$logLik
  gLLbyErr.KFAS    <- sub.LogLik(oEstimated.KFAS$v[,1], t(oEstimated.KFAS$F), nStateVar)
  # report
  agQ <- diag(oFitted.KFAS$model$Q[,,1])
  cat(sprintf("hat(sigma^2_epsilon) = %f\n", drop(oFitted.KFAS$model$H)))
  cat(sprintf("hat(sigma^2_xi)      = %e\n", agQ[1]))
  cat(sprintf("hat(sigma^2_zeta)    = %f\n", agQ[2]))
  cat(sprintf("hat(mu_1)            = %f\n", oEstimated.KFAS$alphahat[1,1]))
  cat(sprintf("hat(nu_1)            = %f\n", oEstimated.KFAS$alphahat[1,2]))
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
  cat("\n") 
  
}
# - - - - - - - - - - - - - - - - - 
exercise.Chapter3.4b <- function(){
  # the rest of section 3.4
  # local linear trend model with deterministic level & stochastic slope
  
  library(dlm)
  library(KFAS)
  
  # constants - - - - - -
  sFILENAME <- "NorwayFinland.txt"
  gCKLogLik <- 0.7864746
  # - - - - - - - - - - - 
  
  # read data
  dfData <- read.table(sFILENAME, skip=1)
  anTime <- as.integer(dfData[,1])
  stopifnot(anTime == seq(from=anTime[1], to=anTime[1]+length(anTime)-1))
  tsData <- ts(log(dfData[,3]), start=anTime[1])
  n      <- length(tsData)
  
  # model specs
  nStateVar   <- 2
  nHyperParam <- 2
  
  # init value for variances (epsilon, zeta)
  agInit     <- c(var(tsData),0.001)  
  
  cat("\n")
  
  cat("estimated by dlm package: \n")
  funModel <- function(parm){
    dlmModPoly(order = 2, dV = exp(parm[1]), dW = c(0, exp(parm[2])))
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
  oFiltered.DLM     <- dlmFilter(tsData, oFitted.DLM)
  agStdPredErr.DLM  <- residuals(oFiltered.DLM, type = c("standardized"), sd=FALSE)
  lPredErr.DLM      <- residuals(oFiltered.DLM, type = c("raw"), sd=TRUE)  
  # Smoothing
  oSmoothed.DLM <- dlmSmooth(tsData, oFitted.DLM)
  tsSmoState.DLM   <- dropFirst(oSmoothed.DLM$s)
  # LL
  gLLbyFun.DLM    <- -oMLE.DLM$value - 0.5*n*log(2*pi)
  gLLbyErr.DLM    <- sub.LogLik(lPredErr.DLM$res, lPredErr.DLM$sd^2, nStateVar)
  # report
  agW <- diag(oFitted.DLM$W)
  cat(sprintf("hat(sigma^2_epsilon) = %f\n", drop(oFitted.DLM$V)))
  cat(sprintf("hat(sigma^2_zeta)    = %f\n", agW[2]))
  cat(sprintf("hat(mu_1)            = %f\n", oSmoothed.DLM$s[2,1]))
  cat(sprintf("hat(nu_1)            = %f\n", oSmoothed.DLM$s[2,2]))
  cat("\n") 
  
  cat("estimated by KFAS package: \n")
  oModel.KFAS <- SSModel(
    tsData ~ SSMtrend(degree=2, Q=list(matrix(0), matrix(NA))), 
    H=matrix(NA)
  )
  ## print(oModel.KFAS) -> states are (level, slope)
  # Fitting
  oFitted.KFAS <- fitSSM(
    oModel.KFAS, 
    inits = log(agInit)
  )
  stopifnot(oFitted.KFAS$optim.out$convergence == 0)
  # Filtering & Smoothing
  oEstimated.KFAS <- KFS(oFitted.KFAS$model)
  agStdPredErr.KFAS  <- oEstimated.KFAS$v[,1] / sqrt(t(oEstimated.KFAS$F))
  # LL
  gLLbyFun.KFAS    <- oEstimated.KFAS$logLik
  gLLbyErr.KFAS    <- sub.LogLik(oEstimated.KFAS$v[,1], t(oEstimated.KFAS$F), nStateVar)
  # report
  agQ <- diag(oFitted.KFAS$model$Q[,,1])
  cat(sprintf("hat(sigma^2_epsilon) = %f\n", drop(oFitted.KFAS$model$H)))
  cat(sprintf("hat(sigma^2_zeta)    = %f\n", agQ[2]))
  cat(sprintf("hat(mu_1)            = %f\n", oEstimated.KFAS$alphahat[1,1]))
  cat(sprintf("hat(nu_1)            = %f\n", oEstimated.KFAS$alphahat[1,2]))
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
  cat("\n") 
  
  # plot
  sub.plotFigure_twoseries(
    tsData, 
    tsSmoState.DLM[,1], 
    "Figure 3.5a",
    c("log fatalities Finland", "deterministic level, stochastic slope")
  )
  sub.plotFigure_twoseries(
    NULL, 
    tsSmoState.DLM[,2],
    "Figure 3.5b",
    "stochastic slope"
  )
  sub.plotFigure_residual(
    tsData - tsSmoState.DLM[,1], 
    "Figure 3.6",
    "irregular"
  )
  
  # plot for section 8.5, Figure 8.14 - - - - - 
  # observed
  tsData.Obs <- ts(c(tsData, rep(NA, 5)), start=anTime[1])
  # variance in filter
  lFltVar           <- dlmSvd2var(oFiltered.DLM$U.R, oFiltered.DLM$D.R)
  agLevelVar        <- sapply(lFltVar, function(x) x[1,1])
  # forecast
  oForecast         <- dlmForecast(oFiltered.DLM, nAhead = 5)
  agLevelVar.Future <- unlist(oForecast$Q)
  
  # filtered level
  tsData.Level    <- ts(
    c(oFiltered.DLM$f, oForecast$f),
    start=anTime[1]
  )
  # confidence interval
  tsData.Width    <- ts(
    qnorm(0.05, lower = FALSE) * sqrt(c(agLevelVar, agLevelVar.Future)),
    start=anTime[1]
  )
  
  sub.plotFigure_interval(
    tsData.Obs, 
    window(tsData.Level, start=1972), 
    window(tsData.Width, start=1972), 
    "Figure 8.14", 
    c("log fatalities in Finland", "filtered trend and forecasts")
  )
  
}
