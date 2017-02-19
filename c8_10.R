# - - - - - - - - - - - - - - - - - 
exercise.Chapter8.6a <- function(){
  # log UK driver KSI before SEAT-BELT law
  # stochastic level, stochastic season, log UK petrol price
  
  library(dlm)
  library(KFAS)
  
  # constants - - - - - -
  sFILENAME    <- "UKdriversKSI.txt"
  sXFILENAME   <- "logUKpetrolprice.txt"
  gCKLogLik    <- 0.9556575
  # - - - - - - - - - - - 
  
  # reading data
  dfData <- read.table(sFILENAME, skip=1)
  # complete data
  ## tsData <- ts(log(dfData[1:169,1]), frequency=12, start=c(1969, 1))
  ## n      <- length(tsData)
  # short data
  tsShort <- ts(log(dfData[1:169,1]), frequency=12, start=c(1969, 1))
  nShort  <- length(tsShort)
  # explanatory data
  dfX1Data <- read.table(sXFILENAME, skip=1)
  x <- ts(dfX1Data[,1], frequency=12, start=c(1969, 1))
  
  # model specs
  nStateVar   <- 14
  nHyperParam <- 3
  
  # init value for variances (epsilon, xi, omega)
  agInit     <- c(var(tsShort), 0.0001, 0.0001)
  
  cat("\n")
  
  cat("estimated by dlm package:\n")
  funModel <- function(parm){
    m1 <- dlmModPoly(order = 1, dV=exp(parm[1]), dW=exp(parm[2]))
    m2 <- dlmModSeas(12, dV=0, dW=c(exp(parm[3]), rep(0,10)))
    m3 <- dlmModReg(x[1:nShort], dV=0, addInt=FALSE)
    return( m1 + m2 + m3 )
  }
  ## print(funModel(log(c(123, 456, 789))))
  # Fitting
  oMLE.DLM <- dlmMLE(
    tsShort, 
    parm  = log(agInit), 
    build = funModel
  )
  stopifnot(oMLE.DLM$convergence == 0)
  oFitted.DLM    <- funModel(oMLE.DLM$par)
  # Filtering
  oFiltered.DLM <- dlmFilter(tsShort, oFitted.DLM)
  agStdPredErr.DLM  <- residuals(oFiltered.DLM, type = c("standardized"), sd=FALSE)
  lPredErr.DLM      <- residuals(oFiltered.DLM, type = c("raw"), sd=TRUE)  
  # Smoothing
  oSmoothed.DLM <- dlmSmooth(tsShort, oFitted.DLM)
  # LL
  gLLbyFun.DLM    <- -oMLE.DLM$value - 0.5*nShort*log(2*pi)
  gLLbyErr.DLM    <- sub.LogLik(lPredErr.DLM$res, lPredErr.DLM$sd^2, nStateVar)
  # report
  agW           <- diag(oFitted.DLM$W)
  cat(sprintf("hat(sigma^2_epsilon) = %f\n", drop(oFitted.DLM$V)))
  cat(sprintf("hat(sigma^2_xi)      = %e\n", agW[1]))
  cat(sprintf("hat(sigma^2_omega)   = %e\n", agW[2]))
  cat(sprintf("hat(mu_1)            = %f\n", oSmoothed.DLM$s[2,1]))
  cat(sprintf("hat(gamma_1)         = %f\n", oSmoothed.DLM$s[2,2]))
  cat(sprintf("hat(beta_1)          = %f\n", oSmoothed.DLM$s[1,13]))
  cat("\n") 
  
  cat("estimated by KFAS package:\n")
  oModel.KFAS <- SSModel(
    tsShort ~ SSMtrend(degree=1, Q=list(matrix(NA))) + 
      SSMseasonal(12, sea.type="dummy", Q=matrix(NA), n=n) + 
      SSMregression(~ -1 + x[1:nShort], Q=matrix(0)), 
    H=matrix(NA)
  )
  ## print(oModel.KFAS)   # -> states are (x, level, seasondummies)
  ## print(oModel.KFAS$Q) # -> Q is (3,3) matrix where diag is (0,NA,NA)
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
  cat(sprintf("hat(sigma^2_xi)      = %e\n", agQ[2]))
  cat(sprintf("hat(sigma^2_omega)   = %e\n", agQ[3]))
  cat(sprintf("hat(mu_1)            = %f\n", coef(oEstimated.KFAS)[1,2]))
  cat(sprintf("hat(gamma_1)         = %f\n", coef(oEstimated.KFAS)[1,3]))
  cat(sprintf("hat(beta_1)          = %f\n", coef(oEstimated.KFAS)[1,1]))
  cat("\n") 
  
  cat("max loglikelihood and AIC:\n")
  sub.ShowLogLik (
    nShort, nStateVar, nHyperParam, 
    gLLbyFun.DLM, gLLbyErr.DLM, gLLbyFun.KFAS, gLLbyErr.KFAS, gCKLogLik*nShort
  )
  cat("\n") 
  
}
# - - - - - - - - - - - - - - - - - 
exercise.Chapter8.6b <- function(){
  # log UK driver KSI before SEAT-BELT law
  # stochastic level, deterministic season, log UK petrol price
  
  library(dlm)
  library(KFAS)
  
  # constants - - - - - -
  sFILENAME    <- "UKdriversKSI.txt"
  sXFILENAME   <- "logUKpetrolprice.txt"
  gCKLogLik    <- 0.9555823
  # - - - - - - - - - - - 
  
  # reading data
  dfData <- read.table(sFILENAME, skip=1)
  # complete data
  ## tsData <- ts(log(dfData[,1]), frequency=12, start=c(1969, 1))
  ## n      <- length(tsData)
  # short data
  tsShort    <- ts(
    log(dfData[1:169,1]), 
    frequency=12, start=c(1969, 1)
  )
  nShort     <- length(tsShort)
  # explanatory data
  dfX1Data <- read.table(sXFILENAME, skip=1)
  x <- ts(dfX1Data[,1], frequency=12, start=c(1969, 1))
  
  # model specs
  nStateVar   <- 13
  nHyperParam <- 2
  
  # init value for variances (epsilon, xi)
  agInit     <- c(var(tsShort), 0.0001)
  
  cat("\n")
  
  cat("estimated by dlm package:\n")
  funModel <- function(parm){
    m1 <- dlmModPoly(order = 1, dV=exp(parm[1]), dW=exp(parm[2]))
    m2 <- dlmModSeas(12, dV=0, dW=rep(0,11))
    m3 <- dlmModReg(x[1:nShort], dV=0, addInt=FALSE)
    return( m1 + m2 + m3 )
  }
  # Fitting
  oMLE.DLM <- dlmMLE(
    tsShort, 
    parm  = log(agInit), 
    build = funModel
  )
  stopifnot(oMLE.DLM$convergence == 0)
  oFitted.DLM    <- funModel(oMLE.DLM$par)
  # Filtering
  oFiltered.DLM <- dlmFilter(tsShort, oFitted.DLM)
  agStdPredErr.DLM  <- residuals(oFiltered.DLM, type = c("standardized"), sd=FALSE)
  lPredErr.DLM      <- residuals(oFiltered.DLM, type = c("raw"), sd=TRUE)  
  # Smoothing
  oSmoothed.DLM <- dlmSmooth(tsShort, oFitted.DLM)
  # LL
  gLLbyFun.DLM    <- -oMLE.DLM$value - 0.5*nShort*log(2*pi)
  gLLbyErr.DLM    <- sub.LogLik(lPredErr.DLM$res, lPredErr.DLM$sd^2, nStateVar)
  cat("LL by dlm (to compare with exercise.Chapter8.6c()):", oMLE.DLM$value, "\n")
  # report
  agW           <- diag(oFitted.DLM$W)
  cat(sprintf("hat(sigma^2_epsilon) = %f\n", drop(oFitted.DLM$V)))
  cat(sprintf("hat(sigma^2_xi)      = %e\n", agW[1]))
  cat(sprintf("hat(mu_1)            = %f\n", oSmoothed.DLM$s[2,1]))
  cat(sprintf("hat(gamma_1)         = %f\n", oSmoothed.DLM$s[2,2]))
  cat(sprintf("hat(beta_1)          = %f\n", oSmoothed.DLM$s[1,13]))
  cat("\n") 
  
  cat("estimated by KFAS package:\n")
  oModel.KFAS <- SSModel(
    tsShort ~ SSMtrend(degree=1, Q=list(matrix(NA))) + 
      SSMseasonal(12, sea.type="dummy", Q=matrix(0), n=n) + 
      SSMregression(~ -1 + x[1:nShort], Q=matrix(0)), 
    H=matrix(NA)
  )
  ## print(oModel.KFAS)   # -> states are (x, level, seasondummies)
  ## print(oModel.KFAS$Q) # -> Q is (3,3) matrix where diag is (0,NA,0)
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
  cat("LL by KFAS (to compare with exercise.Chapter8.6c()):", oEstimated.KFAS$logLik, "\n")
  # report
  agQ            <- diag(oFitted.KFAS$model$Q[,,1])
  cat(sprintf("hat(sigma^2_epsilon) = %f\n", drop(oFitted.KFAS$model$H)))
  cat(sprintf("hat(sigma^2_xi)      = %e\n", agQ[2]))
  cat(sprintf("hat(mu_1)            = %f\n", coef(oEstimated.KFAS)[1,2]))
  cat(sprintf("hat(gamma_1)         = %f\n", coef(oEstimated.KFAS)[1,3]))
  cat(sprintf("hat(beta_1)          = %f\n", coef(oEstimated.KFAS)[1,1]))
  cat("\n") 
  
  cat("max loglikelihood and AIC:\n")
  sub.ShowLogLik (
    nShort, nStateVar, nHyperParam, 
    gLLbyFun.DLM, gLLbyErr.DLM, gLLbyFun.KFAS, gLLbyErr.KFAS, gCKLogLik*nShort
  )
  cat("\n") 
  
  # plot
  oModelNew.KFAS <- SSModel(
    rep(NA, length(x[-(1:nShort)])) ~ SSMtrend(degree=1, Q=list(matrix(0))) + 
      SSMseasonal(12, sea.type="dummy", Q=matrix(0)) + 
      SSMregression(~ -1 + x[-(1:nShort)], Q=matrix(0))
  )
  oModelNew.KFAS$H <- oFitted.KFAS$model$H
  oModelNew.KFAS$Q <- oFitted.KFAS$model$Q
  tsPredict.KFAS <- predict(
    oFitted.KFAS$model, 
    oModelNew.KFAS, 
    interval = "confidence", 
    level = 0.90
  )  
  sub.plotFigure_interval(
    NULL, 
    ts(tsPredict.KFAS[,1], frequency = 12, start=c(1983, 1)),
    ts(tsPredict.KFAS[,1]-tsPredict.KFAS[,2], frequency = 12, start=c(1983, 1)),
    "Figure 8.15 (by predict() of KFAS)", 
    c("forecasts +/- 1.64SE")
  )

}
# - - - - - - - - - - - - - - - - - 
exercise.Chapter8.6c <- function(){
  
  library(dlm)
  library(KFAS)
  
  # constants - - - - - -
  sFILENAME    <- "UKdriversKSI.txt"
  sXFILENAME   <- "logUKpetrolprice.txt"
  # - - - - - - - - - - - 
  
  # reading data
  dfData <- read.table(sFILENAME, skip=1)
  # complete data
  ## tsData <- ts(log(dfData[,1]), frequency=12, start=c(1969, 1))
  ## n      <- length(tsData)
  # data with missing
  tsWithMiss <- ts(
    c( log(dfData[1:169,1]), rep(NA,nrow(dfData)-169) ), 
    frequency=12, start=c(1969, 1)
  )
  nWithMiss  <- length(tsWithMiss)
  # explanatory data
  dfX1Data <- read.table(sXFILENAME, skip=1)
  x <- ts(dfX1Data[,1], frequency=12, start=c(1969, 1))
  
  # model specs
  nStateVar   <- 13
  nHyperParam <- 2
  
  # init value for variances (epsilon, xi)
  agInit     <- c(var(tsWithMiss, na.rm=TRUE), 0.0001)
  
  cat("\n")
  
  # dlm - - - - - 
  funModel <- function(parm){
    m1 <- dlmModPoly(order = 1, dV=exp(parm[1]), dW=exp(parm[2]))
    m2 <- dlmModSeas(12, dV=0, dW=rep(0,11))
    m3 <- dlmModReg(x, dV=0, addInt=FALSE)
    return( m1 + m2 + m3 )
  }
  # fitting
  oMLE.DLM <- dlmMLE(
    tsWithMiss, 
    parm  = log(agInit), 
    build = funModel
  )
  stopifnot(oMLE.DLM$convergence == 0)
  oFitted.DLM  <- funModel(oMLE.DLM$par)
  cat("LL by dlm (to compare with exercise.Chapter8.6b()):", oMLE.DLM$value, "\n")
  # filtering
  oFiltered.DLM <- dlmFilter(tsWithMiss, oFitted.DLM)
  # get variances
  lFltVar       <- dlmSvd2var(oFiltered.DLM$U.R, oFiltered.DLM$D.R)
  agLevelVar    <- sapply(lFltVar, function(x) x[1,1])
  agSeasonVar   <- sapply(lFltVar, function(x) x[2,2])
  agRegVar      <- sapply(lFltVar, function(x) x[13,13])
  agSignalVar.DLM <- agLevelVar + agSeasonVar + x * agRegVar
  
  # KFAS - - - - - 
  oModel.KFAS <- SSModel(
    tsWithMiss ~ SSMtrend(degree=1, Q=list(matrix(NA))) + 
      SSMseasonal(12, sea.type="dummy", Q=matrix(0), n=n) + 
      SSMregression(~ -1 + x, Q=matrix(0)), 
    H=matrix(NA)
  )
  # Fitting
  oFitted.KFAS   <- fitSSM(
    oModel.KFAS, 
    inits  = log(agInit)
  )
  stopifnot(oFitted.KFAS$optim.out$convergence == 0)
  # Filtering 
  oEstimated.KFAS <- KFS(oFitted.KFAS$model, filtering=c("mean", "state"), smoothing="none")
  cat("LL by KFAS (to compare with exercise.Chapter8.6b()):", oEstimated.KFAS$logLik, "\n")
  # get variances
  agSignalVar.KFAS <- sapply(oEstimated.KFAS$P_mu, function(x) drop(x))
  
  tsSignalVar <- ts(
    cbind(agSignalVar.DLM, agSignalVar.KFAS), 
    frequency = 12, start=c(1969, 1)
  )
  plot(
    window(tsSignalVar, start=c(1975, 1)), 
    plot.type="single", 
    main = "variance of filtered signal", 
    ylab = "variance",
    col = c("black", "red")
  )
  legend(
    "top", 
    legend = c("dlm (maybe misspecified)", "KFAS"),
    lty = c(1, 1),
    col = c("black", "red"),
    cex = .8
  )
  
  tsTotal.Signal <- oEstimated.KFAS$m
  tsTotal.Width <- ts(
    qnorm(0.05, lower = FALSE) * sqrt(agSignalVar.KFAS),
    frequency = 12, start=c(1969, 1)
  )
  sub.plotFigure_interval(
    window(tsWithMiss, start=c(1981, 1)),
    window(tsTotal.Signal, start=c(1981,1)), 
    window(tsTotal.Width, start=c(1981,1)), 
    "Figure 8.15 (filtering by KFAS)", 
    c("log UK drivers KSI", "forecasts +/- 1.64SE")
  )
  
}
# - - - - - - - - - - - - - - - - - 
exercise.Chapter8.7 <- function(){
  
  library(dlm)
  library(KFAS)
  
  # constants - - - - - -
  sFILENAME <- "UKdriversKSI.txt"
  gCKLogLik <- 0.9363361
  # - - - - - - - - - - - 
  
  # read data
  dfData <- read.table(sFILENAME, skip=1)
  agData <- log(dfData[,1])
  agData[c(48:62, 120:140)] <- NA
  tsData <- ts(agData, frequency=12, start=c(1969, 1))
  n      <- length(tsData)
  
  # model specs
  nStateVar   <- 12
  nHyperParam <- 2
  
  # init value for variances (epsilon, xi)
  agInit     <- c(var(tsData, na.rm=TRUE), 0.001)  
  
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
  # Smoothing
  oSmoothed.DLM    <- dlmSmooth(tsData, oFitted.DLM)
  tsSmoState.DLM   <- dropFirst(oSmoothed.DLM$s)
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
  # Fitting  
  oFitted.KFAS   <- fitSSM(
    oModel.KFAS, 
    inits  = log(agInit) 
  )
  stopifnot(oFitted.KFAS$optim.out$convergence == 0)
  # Filtering & Smoothing
  oEstimated.KFAS <- KFS(oFitted.KFAS$model)
  # report
  cat(sprintf("hat(sigma^2_epsilon) = %f\n", drop(oFitted.KFAS$model$H)))
  cat(sprintf("hat(sigma^2_xi)      = %f\n", drop(oFitted.KFAS$model$Q[,,1])))
  cat(sprintf("hat(mu_1)            = %f\n", coef(oEstimated.KFAS)[1,1]))
  cat(sprintf("hat(gamma_1)         = %f\n", coef(oEstimated.KFAS)[1,2]))
  cat("\n") 
  
  # variance of smoothed states
  lEstVar       <- dlmSvd2var(oSmoothed.DLM$U.S, oSmoothed.DLM$D.S)[-1]
  agLevelVar    <- sapply(lEstVar, function(x) x[1,1])
  agSeasonVar   <- sapply(lEstVar, function(x) x[2,2])
  
  tsState.KFAS <- oEstimated.KFAS$alphahat
  agLevelVar.KFAS  <- oEstimated.KFAS$V[1,1,]
  agSeasonVar.KFAS <- oEstimated.KFAS$V[2,2,]
  agSignalVar.KFAS <- sapply(oEstimated.KFAS$V_mu, function(x) drop(x))
  
  # plot
  # level variance
  sub.plotFigure_twoseries(
    NULL, 
    ts(agLevelVar, frequency=12, start=c(1969, 1)), 
    "Figure 8.17 (by dlm)", 
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
    "Figure 8.18 (by dlm)", 
    c("log UK drivers KSI", "stochatic level +/- 1.64SE")
  )
  # season
  sub.plotFigure_twoseries(
    NULL, 
    ts(
      agSeasonVar, 
      frequency=12, start=c(1969, 1)
    ), 
    "Figure 8.19 (var(seasonal) by dlm)", 
    "seasonal estimation error variance"
  )
  sub.plotFigure_twoseries(
    NULL, 
    ts(
      agSeasonVar.KFAS, 
      frequency=12, start=c(1969, 1)
    ),
    "Figure 8.19 (var(seasonal) by KFAS)", 
    "seasonal estimation error variance"
  )
  sub.plotFigure_twoseries(
    NULL, 
    ts(
      agSignalVar.KFAS - agLevelVar.KFAS, 
      frequency=12, start=c(1969, 1)
    ),
    "Figure 8.19 (var(signal)-var(level) by KFAS)", 
    "seasonal estimation error variance"
  )
  
  # season interval
  sub.plotFigure_interval(
    NULL, 
    window(tsSmoState.DLM[,2], start=c(1971,1), end=c(1975,1)),
    window(
      ts(
        qnorm(0.05, lower = FALSE) * sqrt(agSeasonVar), 
        frequency=12, start=c(1969, 1)
      ),
      end=c(1975,1)
    ),
    "Figure 8.20 (by dlm)", 
    "deterministic seasonal +/- 1.64SE"
  )
  sub.plotFigure_residual(
    tsData - (tsSmoState.DLM[,1]+tsSmoState.DLM[,2]), 
    "Figure 8.21 (by dlm)",
    "irregular"
  )
}
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
