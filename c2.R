# - - - - - - - - - - - - - - - - - 
exercise.Chapter2.1 <- function(){
  
  library(dlm)
  library(KFAS)
  
  # constants - - - - - -
  sDATAFILENAME   <- "UKdriversKSI.txt"
  gCKLogLik       <- 0.3297597  # AIC in C&K book
  # - - - - - - - - - - - 

  # read data
  dfData <- read.table(sDATAFILENAME, skip=1)
  tsData <- ts(log(dfData[,1]), frequency=12, start=c(1969, 1))
  n      <- length(tsData)
  
  # model specs
  nStateVar   <- 1  # number of state variables
  nHyperParam <- 1  # number of hyper parameters
  
  # init value for (sigma^2_epsilon)
  gInit <- var(tsData)
  
  cat("\n") 
  
  cat("estimated by dlm package: \n")
  funModel <- function(parm){
    dlmModPoly(order = 1, dV = exp(parm), dW = 0)
  }
  # Fitting
  oMLE.DLM <- dlmMLE(
    tsData, 
    parm  = log(gInit), 
    build = funModel
  )
  stopifnot(oMLE.DLM$convergence == 0)
  oFitted.DLM      <- funModel(oMLE.DLM$par)
  # Filtering
  oFiltered.DLM    <- dlmFilter(tsData, oFitted.DLM)
  lPredErr.DLM     <- residuals(oFiltered.DLM, type = c("raw"), sd=TRUE)  
  agStdPredErr.DLM <- residuals(oFiltered.DLM, type = c("standardized"), sd=FALSE)
  # Smoothing
  oSmoothed.DLM    <- dlmSmooth(tsData, oFitted.DLM)
  tsSmoState.DLM   <- dropFirst(oSmoothed.DLM$s)
  # LL
  gLLbyFun.DLM     <- -oMLE.DLM$value - 0.5*n*log(2*pi)
  gLLbyErr.DLM     <- sub.LogLik(lPredErr.DLM$res, lPredErr.DLM$sd^2, nStateVar)
  # report
  cat(sprintf("hat(sigma^2_epsilon) = %f\n", drop(oFitted.DLM$V)))
  cat(sprintf("hat(mu_1)            = %f\n", oSmoothed.DLM$s[2]))
  cat("\n") 
  
  cat("estimated by KFAS package: \n")
  oModel.KFAS <- SSModel(
    tsData ~ SSMtrend(degree=1, Q=list(matrix(0))), 
    H = matrix(NA)
  )
  # Fitting
  oFitted.KFAS <- fitSSM(
    oModel.KFAS, 
    inits  = log(gInit), 
    method = "BFGS"
  )
  # Filtering & Smoothing
  oEstimated.KFAS   <- KFS(oFitted.KFAS$model)
  agStdPredErr.KFAS <- oEstimated.KFAS$v[,1] / t(sqrt(oEstimated.KFAS$F))
  # LL
  gLLbyFun.KFAS     <- oEstimated.KFAS$logLik
  gLLbyErr.KFAS     <- sub.LogLik(oEstimated.KFAS$v[,1], t(oEstimated.KFAS$F), nStateVar)
  # report
  cat(sprintf("hat(sigma^2_epsilon) = %f\n", drop(oFitted.KFAS$model$H)))
  cat(sprintf("hat(mu_1)            = %f\n", oEstimated.KFAS$alphahat[1]))
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
    tsSmoState.DLM, 
    "Figure 2.1", 
    c("log UK drivers KSI", "deterministric level")
  )
  sub.plotFigure_residual(
    tsData - tsSmoState.DLM, 
    "Figure 2.2",
    "irregular"
  )
}
# - - - - - - - - - - - - - - - - - 
exercise.Chapter2.2 <- function(){
  
  library(dlm)
  library(KFAS)
  
  # constants - - - - - -
  sFILENAME <- "UKdriversKSI.txt"
  gCKLogLik    <- 0.6451960  # AIC in C&K book
  # - - - - - - - - - - - 
  
  # read data
  dfData <- read.table(sFILENAME, skip=1)
  tsData <- ts(log(dfData[,1]), frequency=12, start=c(1969, 1))
  n      <- length(tsData)
  
  # model specs
  nStateVar   <- 1
  nHyperParam <- 2
  
  # init value for (sigma^2_epsilon, sigma^2_xi)
  agInit     <- c(var(tsData), 0.001) 
  
  cat("\n")
  
  cat("estimated by dlm package:\n")
  funModel <- function(parm){
    dlmModPoly(order = 1, dV = exp(parm[1]), dW = exp(parm[2]))
  }
  # Fitting
  oMLE.DLM <- dlmMLE(
    tsData, 
    parm  = log(agInit), 
    build = funModel
  )
  stopifnot(oMLE.DLM$convergence == 0)
  oFitted.DLM      <- funModel(oMLE.DLM$par)
  # Filtering
  oFiltered.DLM    <- dlmFilter(tsData, oFitted.DLM)
  lPredErr.DLM     <- residuals(oFiltered.DLM, type = c("raw"), sd=TRUE)  
  agStdPredErr.DLM <- residuals(oFiltered.DLM, type = c("standardized"), sd=FALSE)
  # Smoothing
  oSmoothed.DLM    <- dlmSmooth(tsData, oFitted.DLM)
  tsSmoState.DLM   <- dropFirst(oSmoothed.DLM$s)
  # LL
  gLLbyFun.DLM     <- -oMLE.DLM$value - 0.5*n*log(2*pi)
  gLLbyErr.DLM     <- sub.LogLik(lPredErr.DLM$res, lPredErr.DLM$sd^2, nStateVar)
  # report
  cat(sprintf("hat(sigma^2_epsilon) = %f\n", drop(oFitted.DLM$V)))
  cat(sprintf("hat(sigma^2_xi)      = %f\n", drop(oFitted.DLM$W)))
  cat(sprintf("hat(mu_1)            = %f\n", oSmoothed.DLM$s[2]))
  cat("\n") 
  
  cat("estimated by KFAS package:\n")
  oModel.KFAS <- SSModel(
    tsData ~ SSMtrend(degree=1, Q=list(matrix(NA))), 
    H=matrix(NA)
  )
  # Fitting
  oFitted.KFAS   <- fitSSM(
    oModel.KFAS, 
    inits = log(agInit)
  )
  # Filtering & Smoothing
  oEstimated.KFAS   <- KFS(oFitted.KFAS$model)
  agStdPredErr.KFAS <- oEstimated.KFAS$v[,1] / sqrt(t(oEstimated.KFAS$F))
  # LL
  gLLbyFun.KFAS     <- oEstimated.KFAS$logLik
  gLLbyErr.KFAS     <- sub.LogLik(oEstimated.KFAS$v[,1], t(oEstimated.KFAS$F), nStateVar)
  # report
  cat(sprintf("hat(sigma^2_epsilon) = %f\n", drop(oFitted.KFAS$model$H)))
  cat(sprintf("hat(sigma^2_xi)      = %f\n", drop(oFitted.KFAS$model$Q)))
  cat(sprintf("hat(mu_1)            = %f\n", oEstimated.KFAS$alphahat[1]))
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
    tsSmoState.DLM, 
    "Figure 2.3", 
    c("log UK drivers KSI", "stochastic level")
  )
  sub.plotFigure_residual(
    tsData - tsSmoState.DLM, 
    "Figure 2.4",
    "irregular"
  )
}
# - - - - - - - - - - - - - - - - - 
exercise.Chapter2.3 <- function(){
  
  library(dlm)
  library(KFAS)
  
  # constants - - - - - -
  sFILENAME <- "NorwayFinland.txt"
  gCKLogLik <- 0.8468622
  # - - - - - - - - - - - 
  
  # read data
  dfData <- read.table(sFILENAME, skip=1)
  anTime <- as.integer(dfData[,1])
  stopifnot(anTime == seq(from=anTime[1], to=anTime[1]+length(anTime)-1))
  tsData <- ts(log(dfData[,2]), start=anTime[1])
  n      <- length(tsData)

  # model specs
  nStateVar   <- 1
  nHyperParam <- 2
  
  # init value for (sigma^2_epsilon, sigma^2_xi)
  agInit <- c(var(tsData), 0.001)
  
  cat("\n")
  
  cat("estimated by dlm package:\n")
  funModel <- function(parm){
    dlmModPoly(order = 1, dV = exp(parm[1]), dW = exp(parm[2]))
  }
  # Fitting
  oMLE.DLM <- dlmMLE(
    tsData, 
    parm  = log(agInit), 
    build = funModel
  )
  stopifnot(oMLE.DLM$convergence == 0)
  oFitted.DLM       <- funModel(oMLE.DLM$par)
  # Filtering
  oFiltered.DLM    <- dlmFilter(tsData, oFitted.DLM)
  lPredErr.DLM     <- residuals(oFiltered.DLM, type = c("raw"), sd=TRUE)  
  agStdPredErr.DLM <- residuals(oFiltered.DLM, type = c("standardized"), sd=FALSE)
  # Smoothing
  oSmoothed.DLM    <- dlmSmooth(tsData, oFitted.DLM)
  tsSmoState.DLM   <- dropFirst(oSmoothed.DLM$s)
  # LL
  gLLbyFun.DLM       <- -oMLE.DLM$value - 0.5*n*log(2*pi)
  gLLbyErr.DLM       <- sub.LogLik(lPredErr.DLM$res, lPredErr.DLM$sd^2, nStateVar)
  # report
  cat(sprintf("hat(sigma^2_epsilon) = %f\n", drop(oFitted.DLM$V)))
  cat(sprintf("hat(sigma^2_xi)      = %f\n", drop(oFitted.DLM$W)))
  cat(sprintf("hat(mu_1)            = %f\n", oSmoothed.DLM$s[2]))
  cat("\n") 
  
  cat("estimated by KFAS package:\n")
  oModel.KFAS <- SSModel(
    tsData ~ SSMtrend(degree=1, Q=list(matrix(NA))), 
    H=matrix(NA)
  )
  # Fitting
  oFitted.KFAS   <- fitSSM(
    oModel.KFAS, 
    inits = log(agInit)
  )
  # Filtering & Smoothing
  oEstimated.KFAS   <- KFS(oFitted.KFAS$model, filtering=c("state", "mean"))
  agStdPredErr.KFAS <- oEstimated.KFAS$v[,1] / sqrt(t(oEstimated.KFAS$F))
  # LL
  gLLbyFun.KFAS    <- oEstimated.KFAS$logLik
  gLLbyErr.KFAS    <- sub.LogLik(oEstimated.KFAS$v[,1], t(oEstimated.KFAS$F), nStateVar)
  # report
  cat(sprintf("hat(sigma^2_epsilon) = %f\n", drop(oFitted.KFAS$model$H)))
  cat(sprintf("hat(sigma^2_xi)      = %f\n", drop(oFitted.KFAS$model$Q)))
  cat(sprintf("hat(mu_1)            = %f\n", oEstimated.KFAS$alphahat[1]))
  cat("\n") 
  
  cat("max loglikelihood and AIC:\n")
  sub.ShowLogLik (
    n, nStateVar, nHyperParam, 
    gLLbyFun.DLM, gLLbyErr.DLM, gLLbyFun.KFAS, gLLbyErr.KFAS, gCKLogLik*n
  )
  cat("\n") 
  
  cat("diagnostics of standardized prediction error (estimated by dlm)\n")
  sub.ShowDiagnostics (
    agStdPredErr.DLM, nStateVar, nHyperParam, 10, c(1,4)
  )
  cat("\n") 
  
  # plot
  sub.plotFigure_twoseries(
    tsData, 
    tsSmoState.DLM, 
    "Figure 2.5", 
    c("log fatalities in Norway", "stochastic level")
  )
  sub.plotFigure_residual(
    tsData - tsSmoState.DLM, 
    "Figure 2.6",
    "irregular"
  )
  
  # plot for section 8.5, Fig8.5 & Fig8.7 - - - - - 
  sub.plotFigure_twoseries(
    dropFirst(oSmoothed.DLM$s), 
    dropFirst(oFiltered.DLM$f), 
    "Figure 8.5", 
    c("smoothed level", "filtered level")
  )
  sub.plotFigure_twoseries(
    NULL, 
    dropFirst(lPredErr.DLM$res), 
    "Figure 8.7a", 
    "prediction errors"
  )
  sub.plotFigure_twoseries(
    NULL, 
    dropFirst(lPredErr.DLM$sd)^2, 
    "Figure 8.7b", 
    "prediction error variance"
  )
  
  # plot for section 8.6, Fig8.13 - - - - - 
  # observed 
  tsActual <- ts(c(tsData, rep(NA, 5)), start=anTime[1])
  # past
  agLevel.Past <- as.numeric(oFiltered.DLM$a)
  agLevelVar.Past   <- sapply(
    dlmSvd2var(oFiltered.DLM$U.R, oFiltered.DLM$D.R), 
    function(x) x[1,1]
  )
  # future
  oForecast <- dlmForecast(oFiltered.DLM, nAhead = 5)
  agLevel.Future    <- as.numeric(oForecast$f)
  agLevelVar.Future <- unlist(oForecast$Q)

  ## cf. 4 possible sources of past data
  ##
  ## 1) filtered level by DLM
  ## print(oFiltered.DLM$a)  
  ## print(agLevelVar.Past)  
  ##
  ## 2) forecast by DLM
  ## print(oFiltered.DLM$f) # equal to 1) 
  ## not sure how to get variances, but it may be agLevelVar.Past
  ##
  ## 3) one step predictions of states by KFAS
  ## print(oEstimated.KFAS$a) 
  ## print(sapply(oEstimated.KFAS$P, function(x) drop(x)))
  ##
  ## 4) filtered estimate of signal by KFAS
  ## print(oEstimated.KFAS$m) 
  ## print(sapply(oEstimated.KFAS$P_mu, function(x) drop(x)))

  # filtered level
  tsLevel    <- ts(
    c(agLevel.Past, agLevel.Future), 
    start=anTime[1]
  )
  tsWidth    <- ts(
    qnorm(0.05, lower = FALSE) * sqrt(c(agLevelVar.Past, agLevelVar.Future)),
    start=anTime[1]
  )
  sub.plotFigure_interval(
    tsActual, 
    dropFirst(tsLevel), 
    dropFirst(tsWidth), 
    "Figure 8.13", 
    c("log fatalities in Norway", "filtered level and forecasts")
  )

}
