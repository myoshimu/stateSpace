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
