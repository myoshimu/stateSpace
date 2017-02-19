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
# - - - - - - - - - - - - - - - - - 
exercise.Chapter4.1 <- function(){
  
  library(dlm)
  library(KFAS)
  
  # constants - - - - - -
  sFILENAME <- "UKdriversKSI.txt"
  gCKLogLik <- 0.4174873
  # - - - - - - - - - - - 
  
  # read data
  dfData <- read.table(sFILENAME, skip=1)
  tsData <- ts(log(dfData[,1]), frequency=12, start=c(1969, 1))
  n      <- length(tsData)
  
  # model specs
  nStateVar   <- 12
  nHyperParam <- 1
  
  # init value for variances (epsilon)
  gInit     <- c(var(tsData))  
  
  cat("\n")
  
  cat("estimated by dlm package:\n")
  funModel <- function(parm){
    m1 <- dlmModPoly(order = 1, dV=exp(parm), dW=0)
    m2 <- dlmModSeas(12, dV = 0, dW=rep(0,11)) 
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
  cat(sprintf("hat(gamma_1})        = %f\n", oSmoothed.DLM$s[2,2]))
  cat("\n") 
  
  cat("estimated by KFAS package:\n")
  oModel.KFAS <- SSModel(
    tsData ~ SSMtrend(degree=1, Q=list(matrix(0))) + 
             SSMseasonal(12, sea.type="dummy"), 
    H=matrix(NA)
  )
  ## print(oModel.KFAS)   -> states are (level, seasondummy1, seasondummy2, ...)
  ## print(oModel.KFAS$Q) -> Q is (1,1) matrix of 0
  # Fitting
  oFitted.KFAS   <- fitSSM(
    oModel.KFAS, 
    inits  = log(gInit), 
    method = "BFGS"
  )
  stopifnot(oFitted.KFAS$optim.out$convergence == 0)
  # Filtering & Smoothing
  oEstimated.KFAS    <- KFS(oFitted.KFAS$model)
  agStdPredErr.KFAS  <- rstandard(oEstimated.KFAS, type="recursive")
  # LL
  gLLbyFun.KFAS    <- oEstimated.KFAS$logLik
  gLLbyErr.KFAS    <- sub.LogLik(oEstimated.KFAS$v[,1], t(oEstimated.KFAS$F), nStateVar)
  # report
  cat(sprintf("hat(sigma^2_epsilon) = %f\n", drop(oFitted.KFAS$model$H)))
  cat(sprintf("hat(mu_1)            = %f\n", coef(oEstimated.KFAS)[1,1]))
  cat(sprintf("hat(gamma_1)         = %f\n", coef(oEstimated.KFAS)[1,2]))
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
    tsSmoState.DLM[,1]+tsSmoState.DLM[,2], 
    "Figure 4.2",
    c("log UK drivers KSI", "deterministic level")
  )
  sub.plotFigure_twoseries(
    tsData, 
    tsSmoState.DLM[,1], 
    "Figure 4.3",
    c("log UK drivers KSI", "deterministic level")
  )
  sub.plotFigure_twoseries(
    NULL, 
    tsSmoState.DLM[,2], 
    "Figure 4.4",
    "deterministic seasonal"
  )
  sub.plotFigure_residual(
    tsData - (tsSmoState.DLM[,1]+tsSmoState.DLM[,2]), 
    "Figure 4.5",
    "irregular"
  )

}
# - - - - - - - - - - - - - - - - - 
exercise.Chapter4.2 <- function(){
  
  library(dlm)
  library(KFAS)
  
  # constants - - - - - -
  sFILENAME <- "UKdriversKSI.txt"
  gCKLogLik <- 0.9369063
  # - - - - - - - - - - - 
  
  # read data
  dfData <- read.table(sFILENAME, skip=1)
  tsData <- ts(log(dfData[,1]), frequency=12, start=c(1969, 1))
  n      <- length(tsData)
  
  # model specs
  nStateVar   <- 12
  nHyperParam <- 3
  
  # init value for variances (epsilon, xi, omega)
  agInit <- c(var(tsData), 0.001, 0.001)  
  
  cat("\n")
  
  cat("estimated by dlm package:\n")
  funModel <- function(parm){
    m1 <- dlmModPoly(order = 1, dV=exp(parm[1]), dW=exp(parm[2]))
    m2 <- dlmModSeas(12, dV=0, dW=c(exp(parm[3]), rep(0,10)))
    return( m1 + m2 )
  }
  ## print(funModel(log(c(123,456,789))))
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
  cat("\n") 
  
  cat("estimated by KFAS package:\n")
  oModel.KFAS <- SSModel(
    tsData ~ SSMtrend(degree=1, Q=list(matrix(NA))) + 
             SSMseasonal(12, sea.type="dummy", Q=matrix(NA), n=n), 
    H=matrix(NA)
  )
  ## print(oModel.KFAS) -> states are (level, seasondummy1, seasondummy2, ...)
  ## print(oModel.KFAS$Q) -> Q is (2,2) matrix where diag is (NA,NA)
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
  cat(sprintf("hat(sigma^2_epsilon) = %f\n", drop(oFitted.KFAS$model$H)))
  cat(sprintf("hat(sigma^2_xi)      = %e\n", agQ[1]))
  cat(sprintf("hat(sigma^2_omega)   = %e\n", agQ[2]))
  cat(sprintf("hat(mu_1)            = %f\n", coef(oEstimated.KFAS)[1,1]))
  cat(sprintf("hat(gamma_1)         = %f\n", coef(oEstimated.KFAS)[1,2]))
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
    "Figure 4.6",
    c("log UK drivers KSI", "stochastic level")
  )
  sub.plotFigure_twoseries(
    NULL, 
    tsSmoState.DLM[,2], 
    "Figure 4.7",
    "stochastic seasonal"
  )
  sub.plotFigure_twoseries(
    NULL, 
    ts(window(tsSmoState.DLM[,2], start=c(1969,1), end=c(1969,12)), frequency = 1),
    "Figure 4.8",
    "seasonal 1969"
  )
  sub.plotFigure_residual(
    tsData - (tsSmoState.DLM[,1]+tsSmoState.DLM[,2]), 
    "Figure 4.9",
    "irregular"
  )
  
}
# - - - - - - - - - - - - - - - - - 
exercise.Chapter4.3 <- function(){
  
  library(dlm)
  library(KFAS)
  
  # constants - - - - - -
  sFILENAME <- "UKdriversKSI.txt"
  gCKLogLik <- 0.9363361
  # - - - - - - - - - - - 
  
  # read data
  dfData <- read.table(sFILENAME, skip=1)
  tsData <- ts(log(dfData[,1]), frequency=12, start=c(1969, 1))
  n      <- length(tsData)
  
  # model specs
  nStateVar   <- 12
  nHyperParam <- 2
  
  # init value for variances (epsilon, xi)
  agInit     <- c(var(tsData), 0.001)  
  
  cat("\n")
  
  cat("estimated by dlm package:\n")
  funModel <- function(parm){
    m1 <- dlmModPoly(order = 1, dV=exp(parm[1]), dW=exp(parm[2]))
    m2 <- dlmModSeas(12, dV=0, dW=rep(0,11))
    return( m1 + m2 )
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
  lPredErr.DLM     <- residuals(oFiltered.DLM, type = c("raw"), sd=TRUE)  
  agStdPredErr.DLM <- residuals(oFiltered.DLM, type = c("standardized"), sd=FALSE)
  # Smoothing
  oSmoothed.DLM    <- dlmSmooth(tsData, oFitted.DLM)
  tsSmoState.DLM   <- dropFirst(oSmoothed.DLM$s)
  # LL
  gLLbyFun.DLM    <- -oMLE.DLM$value - 0.5*n*log(2*pi)
  gLLbyErr.DLM    <- sub.LogLik(lPredErr.DLM$res, lPredErr.DLM$sd^2, nStateVar)
  # report
  agW           <- diag(oFitted.DLM$W)
  cat(sprintf("hat(sigma^2_epsilon) = %f\n", drop(oFitted.DLM$V)))
  cat(sprintf("hat(sigma^2_xi)      = %f\n", agW[1]))
  cat(sprintf("hat(mu_1)            = %f\n", oSmoothed.DLM$s[2,1]))
  cat(sprintf("hat(gamma_1)         = %f\n", oSmoothed.DLM$s[2,2]))
  cat("\n") 
  
  cat("estimated by KFAS package:\n")
  oModel.KFAS <- SSModel(
    tsData ~ SSMtrend(degree=1, Q=list(matrix(NA))) + 
             SSMseasonal(12, sea.type="dummy"), 
    H=matrix(NA)
  )
  ## print(oModel.KFAS) -> states are (level, seasondummy1, seasondummy2, ...)
  ## print(oModel.KFAS$Q)  -> Q is (1,1) matrix of NA
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
  cat(sprintf("hat(sigma^2_epsilon) = %f\n", drop(oFitted.KFAS$model$H)))
  cat(sprintf("hat(sigma^2_xi)      = %f\n", drop(oFitted.KFAS$model$Q[,,1])))
  cat(sprintf("hat(mu_1)            = %f\n", coef(oEstimated.KFAS)[1,1]))
  cat(sprintf("hat(gamma_1)         = %f\n", coef(oEstimated.KFAS)[1,2]))
  cat("\n") 
  
  cat("max loglikelihood and AIC:\n")
  sub.ShowLogLik (
    n, nStateVar, nHyperParam, 
    gLLbyFun.DLM, gLLbyErr.DLM, gLLbyFun.KFAS, gLLbyErr.KFAS, gCKLogLik*n
  )
  cat("\n") 
  
  # plot for section 8.3 - - - - - - - -
  # variance of smoothed states
  lEstVar       <- dlmSvd2var(oSmoothed.DLM$U.S, oSmoothed.DLM$D.S)[-1]
  agLevelVar    <- sapply(lEstVar, function(x) x[1,1])
  agSeasonVar   <- sapply(lEstVar, function(x) x[2,2])
  agTotalVar    <- agLevelVar + agSeasonVar

  # plot
  # level variance
  sub.plotFigure_twoseries(
    NULL, 
    ts(agLevelVar, frequency=12, start=c(1969, 1)), 
    "Figure 8.1", 
    "level estimation error variance"
  )
  # level interval
  sub.plotFigure_interval(
    tsData, 
    tsSmoState.DLM[,1], 
    ts(
      qnorm(0.05, lower = FALSE) * sqrt(agLevelVar), 
      frequency=12, start=c(1969, 1)
    ), 
    "Figure 8.2", 
    c("log UK drivers KSI", "stochatic level +/- 1.64SE")
  )
  # season
  tsSeasonWidth <- ts(
    qnorm(0.05, lower = FALSE) * sqrt(agSeasonVar), 
    frequency=12, start=c(1969, 1)
  )
  sub.plotFigure_interval(
    NULL, 
    window(tsSmoState.DLM[,2],  start=c(1981,1)),
    window(tsSeasonWidth, start=c(1981,1)), 
    "Figure 8.3", 
    c("deterministic seasonal +/- 1.64SE")
  )
  # signal
  tsTotalWidth <- ts(
    qnorm(0.05, lower = FALSE) * sqrt(agTotalVar), 
    frequency=12, start=c(1969, 1)
  )
  sub.plotFigure_interval(
    NULL, 
    window(tsSmoState.DLM[,1]+tsSmoState.DLM[,2], start=c(1981,1)),
    window(tsTotalWidth, start=c(1981,1)), 
    "Figure 8.4 (by dlm)", 
    c("signal +/- 1.64SE")
  )

  # ploting Figure 8.4 again by KFAS - - - - - - - - - -
  tsSmoState.KFAS <- oEstimated.KFAS$muhat
  agTotalVar.KFAS <- sapply(oEstimated.KFAS$V_mu, function(x) drop(x))

  plot(
    ts(
      cbind(agTotalVar, agTotalVar.KFAS), 
      frequency=12, start=c(1969, 1)
    ), 
    plot.type="single", 
    col = c("black", "red"),
    ylab = "variance", 
    main = "Variance of smoothed signal"
  )
  legend(
    "top", 
    legend = c(
      "dlmSvd2var(dlmSmooth(...)$U.S, dlmSmooth(...)$D.S)", 
      "KFS(...)$V_mu"
    ),
    lty = c(1, 1),
    col = c("black", "red"),
    cex = .8
  )
  ## I compared these two series with the output of ssfpack 
  ## and found ssfpack is similar with KFAS, not with ldm
  ## print(agTotalVar[144:192])
  ## print(agTotalVar.KFAS[144:192])
  
  tsTotalWidth.KFAS <- ts(
    qnorm(0.05, lower = FALSE) * sqrt(agTotalVar.KFAS), 
    frequency=12, start=c(1969, 1)
  )
  sub.plotFigure_interval(
    NULL, 
    window(tsSmoState.KFAS, start=c(1981,1)),
    window(tsTotalWidth.KFAS, start=c(1981,1)), 
    "Figure 8.4 (by KFAS)", 
    c("signal +/- 1.64SE")
  )

  # plot for section 8.5 - - - - - - - -
  sub.plotFigure_auxresidual(
    rstandard(oEstimated.KFAS, type="state"), 
    "Figure 8.11a", 
    "Structural level break t-tests"
  )
  sub.plotFigure_auxresidual(
    rstandard(oEstimated.KFAS, type="pearson"), 
    "Figure 8.11b", 
    "Outlier t-tests"
  )
}
# - - - - - - - - - - - - - - - - - 
exercise.Chapter4.4 <- function(){
  
  library(dlm)
  library(KFAS)
  
  # constants - - - - - -
  sFILENAME <- "UKinflation.txt"
  gCKLogLik <- 3.198464
  # - - - - - - - - - - - 
  
  # read data
  dfData <- read.table(sFILENAME, skip=1)
  tsData <- ts(dfData[,1], frequency=4, start=c(1950, 1))
  n      <- length(tsData)
  # plot(tsData)
  
  # model specs
  nStateVar   <- 4
  nHyperParam <- 3
  
  # init value for variances (epsilon, xi, gamma)
  agInit     <- c(var(tsData), 0.00001, 0.00001)  
  
  cat("\n")
  
  cat("estimated by dlm package:\n")
  funModel <- function(parm){
    m1 <- dlmModPoly(order = 1, dV=exp(parm[1]), dW=exp(parm[2]))
    m2 <- dlmModSeas(4, dV=0, dW=c(exp(parm[3]), rep(0,2)))
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
  cat(sprintf("hat(sigma^2_epsilon) = %e\n", drop(oFitted.DLM$V)))
  cat(sprintf("hat(sigma^2_xi)      = %e\n", agW[1]))
  cat(sprintf("hat(sigma^2_gamma)   = %e\n", agW[2]))
  cat(sprintf("hat(mu_1)            = %f\n", oSmoothed.DLM$s[2,1]))
  cat(sprintf("hat(gamma_1)         = %f\n", oSmoothed.DLM$s[2,2]))
  cat(sprintf("hat(mu_208)          = %f\n", oSmoothed.DLM$s[209,1]))
  cat("\n") 
  
  cat("estimated by KFAS package:\n")
  oModel.KFAS <- SSModel(
    tsData ~ SSMtrend(degree=1, Q=list(matrix(NA))) + 
             SSMseasonal(4, sea.type="dummy", Q=matrix(NA), n=n), 
    H=matrix(NA)
  )
  ## print(oModel.KFAS) -> states are (level, seasondummy1, seasondummy2, seasondummy3)
  ## print(oModel.KFAS$Q) -> Q is (2,2) matrix where diag is (NA,NA)
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
  agQ            <- diag(oFitted.KFAS$model$Q[,,1])
  cat(sprintf("hat(sigma^2_epsilon) = %e\n", drop(oFitted.KFAS$model$H)))
  cat(sprintf("hat(sigma^2_xi)      = %e\n", agQ[1]))
  cat(sprintf("hat(sigma^2_gamma)   = %e\n", agQ[2]))
  cat(sprintf("hat(mu_1)            = %f\n", coef(oEstimated.KFAS)[1,1]))
  cat(sprintf("hat(gamma_1)         = %f\n", coef(oEstimated.KFAS)[1,2]))
  cat(sprintf("hat(gamma_208)       = %f\n", coef(oEstimated.KFAS)[208,1]))
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
    "Figure 4.10a",
    c("quarterly price changes in U.K.", "stochastic level")
  )
  sub.plotFigure_twoseries(
    NULL, 
    tsSmoState.DLM[,2], 
    "Figure 4.10b",
    "stchastic seasonal"
  )
  sub.plotFigure_residual(
    tsData - (tsSmoState.DLM[,1]+tsSmoState.DLM[,2]), 
    "Figure 4.10c",
    "irregular"
  )
  
}
