# R programs for 
#    Commandeur, J.J.F. & Koopman, S.J. (2007)
#    "An Introduction to State Space Time Series Analysis,"
#    Oxford University Press.
#
# notes:
#    prefixes: 
#      n  : iNteger
#      g  : numeric (g in "sinGle")
#      b  : binary  
#      a  : vector  (not sure where it comes from)
#      m  : matrix
#      l  : list
#      ts : ts object (Time Series)
#      o  : object of any class
#
# # # # # # # # # # # # # # # # # # #
#       some helper routines        #
# # # # # # # # # # # # # # # # # # #
# - - - - - - - - - - - - -
sub.LogLik <- function(tsNu, tsF, nStateVar){
  # purpose: get log likelihood.
  #          See Section 8.4
  # args:    tsNu:      (ts object)      time series of predicion error
  #          tsF:       (ts object)      time series of prediction error variance
  #          nStateVar: (numeric scalar) number of state variables
  # return:  (numeric scalar) log likelihood
  
  # trap: the length should be same
  stopifnot(length(tsF) == length(tsNu))
  
  # remove the first (nStateVar) timepoints
  agNu <- tsNu[-(1:nStateVar)]
  agF  <- tsF[-(1:nStateVar)]
  
  - 0.5 * length(tsNu) * log(2*pi) - 0.5 * sum( log(agF) + (agNu^2 / agF) )
}
# - - - - - - - - - - - - -- -
sub.AIC <- function(gLogLik, nStateVar, nHyperParam){
  # purpose: get AIC
  #          See Section 2.1.
  # args:    gLogLik:     (numeric scalar) log likelihood
  #          nStateVar:   (integer scalar) number of state variables
  #          nHyperParam: (integer scalar) number of hyper params
  # return:  (numeric scalar) AIC
  
  - 2 * gLogLik + 2 * ( nStateVar + nHyperParam )
  
}
# - - - - - - - - - - - - -- -
sub.BoxLjung <- function(agStdRes, nStateVar, nHyperParam, nLag=15){
  # purpose: Independence test for standardized residuals
  #          See Section 8.5.
  # args:    agStdRes:    (numeric vector) standardized residuals
  #          nStateVar:   (integer scalar) number of state variables
  #          nHyperParam: (integer scalar) number of hyper params
  #          nLag:        (integer scalar) max lag
  # return:  (list)
  #            stat:     (numeric scalar) Q(nLag) 
  #            critical: (numeric scalar) 5% critical value
  # notes:   Though the C&K book suggests all residuals can be used, 
  #          their ssfpack program removes the first (nStateVar) 
  #          residuals.
  
  # call Box.test
  gStat <- Box.test(
    agStdRes[-(1:nStateVar)], # remove the first (nStateVar) timepoints
    lag     = nLag, 
    type    = "Ljung-Box"
  )$statistic
  
  # get critical value
  gCritical <- qchisq(0.95, nLag-nHyperParam+1)
  
  list(
    stat     = gStat, 
    critical = gCritical
  )
}
# - - - - - - - - - - - - -- -
sub.ResidualAcf <- function(agStdRes, nStateVar, nLag=15){
  # purpose: autocorrelation of standardized residuals
  #          See Section 8.5.
  # args:    agStdRes:    (numeric vector) standardized residuals
  #          nStateVar:   (integer scalar) number of state variables
  #          nLag:        (integer scalar) max lag
  # return:  (list)
  #            acf:      (numeric vector) acf (1 ... nLag)
  #            critical: (numeric scalar) 95% confidence limit
  # notes:   Though the C&K book suggests all residuals can be used, 
  #          their ssfpack program removes the first (nStateVar) 
  #          residuals.
  
  list(
    acf = acf(agStdRes[-(1:nStateVar)], plot=FALSE)$acf[-1],
    critical = 2 / sqrt(length(agStdRes))
  )
  
}
# - - - - - - - - - - - - -- -
sub.Homoscedasticity <- function(agStdRes, nStateVar){
  # purpose: Homoscedasticity test of standardized residuals
  #          See the C&K Book, Section 8.5.
  # args:    agStdRes:    (numeric vector) standardized residuals
  #          nStateVar:   (integer scalar) number of state variables
  # return:  (list)
  #            size:  (integer scalar) block size
  #            stat:  (numeric scalar) H(block size)
  #            upper: (numeric scalar) 5% critical value (upper)
  #            lower: (numeric scalar) 5% critical value (lower)
  
  # get leength
  n <- length(agStdRes)
  
  # defint blocksize
  nBlockSize <- round((n - nStateVar)/3)
  gValue     <- sum(agStdRes[(n-nBlockSize+1):n]^2) / 
                  sum(agStdRes[(nStateVar+1):(nStateVar+nBlockSize)]^2)
  
  list(
    size  = nBlockSize,
    stat  = gValue,
    upper = qf(0.975, nBlockSize, nBlockSize), 
    lower = qf(0.025, nBlockSize, nBlockSize)
  )
}
# - - - - - - - - - - - - -- -
sub.JarqueBera <- function(agStdRes, nStateVar){
  # purpose: Normality test of standardized residuals
  #          See Section 8.5.
  # args:    agStdRes:    (numeric vector) standardized residuals
  #          nStateVar:   (integer scalar) number of state variables
  # return:  (list)
  #            stat:     (numeric vector) statistic
  #            critical: (numeric vector) 5% critical value (=5.99)
  # notes:   Though the C&K book suggests all residuals can be used, 
  #          their ssfpack program removes the first (nStateVar) 
  #          residuals.
  
  require(tseries)
  
  list(
    stat     = jarque.bera.test(agStdRes[-(1:nStateVar)])$statistic, 
    critical = qchisq(0.95,2)
  )
}
# - - - - - - - - - - - - -- -
sub.ShowLogLik <- function(
  n, nStateVar, nHyperParam, 
  gLLbyFun.DLM, gLLbyErr.DLM, gLLbyFun.KFAS, gLLbyErr.KFAS, gLL.CK
){
  # purpose: show tables of LogLikelihood(LL) and AIC
  # args:    n:           (numeric scalar) length of time series
  #          nSteteVar:   (integer scalar) number of state variables
  #          nHyperParam: (integer scalar) number of hyper params
  #          gLLbyFun.DLM:  (numeric scalar) LL computed by dlm
  #          gLLbyErr.DLM:  (numeric scalar) LL computed with predition error by dlm
  #          gLLbyFun.KFAS: (numeric scalar) LL computed by KFAS
  #          gLLbyErr.KFAS: (numeric scalar) LL computed with predition error by KFAS
  #          gLL.CK:      (numeric scalar) LL in the CK Book 
  # return:  NULL
  
  asTemplate <- c(
    "-----------------------------------------------------------------",
    "                   LogLik   [diff]     (/n)      AIC    (/n)",
    "-----------------------------------------------------------------",
    "dlm  (by pred.err) %8.4f            (%6.3f) %7.2f (%6.3f)",
    "     (by output)   %8.4f [%+8.3f] (%6.3f) %7.2f (%6.3f)", 
    "KFAS (by pred.err) %8.4f [%+8.3f] (%6.3f) %7.2f (%6.3f)", 
    "     (by output)   %8.4f [%+8.3f] (%6.3f) %7.2f (%6.3f)",
    "The CK Book        %8.4f [%+8.3f] (%6.3f) %7.2f (%6.3f)",
    "-----------------------------------------------------------------"
  )
  cat(sprintf(
    paste(asTemplate, collapse="\n"), 
    # DLM residual
    gLLbyErr.DLM, 
    gLLbyErr.DLM/n, 
    sub.AIC(gLLbyErr.DLM, nStateVar, nHyperParam),
    sub.AIC(gLLbyErr.DLM, nStateVar, nHyperParam)/n,
    # DLM function
    gLLbyFun.DLM, 
    gLLbyFun.DLM - gLLbyErr.DLM,
    gLLbyFun.DLM/n, 
    sub.AIC(gLLbyFun.DLM, nStateVar, nHyperParam),
    sub.AIC(gLLbyFun.DLM, nStateVar, nHyperParam)/n,
    # KFAS residual
    gLLbyErr.KFAS, 
    gLLbyErr.KFAS - gLLbyErr.DLM,
    gLLbyErr.KFAS/n, 
    sub.AIC(gLLbyErr.KFAS, nStateVar, nHyperParam),
    sub.AIC(gLLbyErr.KFAS, nStateVar, nHyperParam)/n,
    # KFAS function
    gLLbyFun.KFAS, 
    gLLbyFun.KFAS - gLLbyErr.DLM,
    gLLbyFun.KFAS/n, 
    sub.AIC(gLLbyFun.KFAS, nStateVar, nHyperParam),
    sub.AIC(gLLbyFun.KFAS, nStateVar, nHyperParam)/n,
    # C&K book
    gLL.CK, 
    gLL.CK - gLLbyErr.DLM,
    gLL.CK/n, 
    sub.AIC(gLL.CK, nStateVar, nHyperParam),
    sub.AIC(gLL.CK, nStateVar, nHyperParam)/n
  ))
  cat("\n")
}
# - - - - - - - - - - - - -- -
sub.ShowDiagnostics <- function(
  agStdRes, nStateVar, nHyperParam, nMaxLag, anACLag
){
  # purpose: Show Diagnostics Table
  # args:    agStdRes:    (numeric vector) standardized residuals
  #          nStateVar:   (integer scalar) number of state variables
  #          nHyperParam: (integer scalar) number of hyper params
  #          nMaxLag:     (integer scalar) max lag (for BoxLjung)
  #          anACLag:     (integer vector) two number of lag (for ACF)
  # return:  NULL
  
  # get various statistics
  lBoxLjung <- sub.BoxLjung(agStdRes, nStateVar, nHyperParam, nMaxLag)
  lACF      <- sub.ResidualAcf(agStdRes, nStateVar, nMaxLag)
  lHomo     <- sub.Homoscedasticity(agStdRes, nStateVar)
  lJB       <- sub.JarqueBera(agStdRes, nStateVar)
  
  # gHomoStat: H statistics, or its reciprocal
  gHomoStat <- ifelse(lHomo$stat > 1, lHomo$stat, 1/lHomo$stat)
  
  asTemplate <- c(
    "------------------------------------------------------",
    "                   stat     value  critical satisfied",
    "------------------------------------------------------",
    "independence      Q  (%2d)  %7.3f   %5.2f     %1s",  # BoxLJung, 4 args
    "                  r  (%2d)  %7.3f  +-%4.2f     %1s", # ACF,      4 args
    "                  r  (%2d)  %7.3f  +-%4.2f     %1s", # ACF,      4 args
    "homoscedasticity  %-3s(%2d)  %7.3f   %5.2f     %1s", # Homo,     5 args
    "normality         N        %7.3f   %5.2f     %1s",   # N,        3 args
    "------------------------------------------------------"
  )
  
  cat(
    sprintf(
      paste(asTemplate, collapse="\n"), 
      # BoxLjung, 4 args
      nMaxLag,    
      lBoxLjung$stat,
      lBoxLjung$critical,
      ifelse(lBoxLjung$stat < lBoxLjung$critical, "+", "-"),
      # ACF, 4 args
      anACLag[1], 
      lACF$acf[anACLag[1]],
      lACF$critical,
      ifelse(abs(lACF$acf[anACLag[1]]) < lACF$critical, "+", "-"),
      # ACF, 4 args
      anACLag[2],
      lACF$acf[anACLag[2]],
      lACF$critical,
      ifelse(abs(lACF$acf[anACLag[2]]) < lACF$critical, "+", "-"),
      # Homo, 5 args
      ifelse(lHomo$stat > 1, "H", "1/H"),
      lHomo$size, 
      gHomoStat, 
      lHomo$upper, 
      ifelse(gHomoStat < lHomo$upper, "+", "-"),
      # N, 3 args
      lJB$stat,  
      lJB$critical,
      ifelse(lJB$stat < lJB$critical, "+", "-")
    )
  )
  cat("\n")
    
}
# - - - - - - - - - - - - -- -
sub.plotFigure_scatterreg <- function(tsData, formula, sTitle, asLegend){
  # purpose: plot scatter plot + regression line (e.g. Figure 1.1)
  # args:    tsData:   (ts object)      time series 
  #          formula:  (formula object) formula for lm()
  #          sTitle:   (string scalar)  title
  #          asLegend: (string vector)  legend of the points and the line
  # return:  NULL
  
  # trap
  stopifnot(is.ts(tsData))
  
  cat("plotting", sTitle, "...\n")
  plot(
    tsData, 
    type = "p",
    pch  = 3,
    xlab = NA, 
    ylab = NA, 
    main = sTitle
  )
  
  oOut.lm <- lm(formula)
  abline(oOut.lm, col="red")
  
  legend(
    "top", 
    legend = asLegend,
    lty = c(-1, 1),
    pch = c(3, -1), 
    col = c("black", "red"),
    cex = .8
  )
}
# - - - - - - - - - - - - -- -
sub.plotFigure_xyscatterreg <- function(tsX, tsY, formula, sTitle, asLegend){
  # purpose: plot x-y scatter plot + regression line (e.g. Figure 5.2)
  # args:    tsY:      (ts object)      time series, dependent variable
  #          tsX:      (ts object)      time series, independent variable
  #          formula:  (formula object) formula for lm()
  #          sTitle:   (string scalar)  title
  #          asLegend: (string vector)  legend of the points and the line
  # return:  NULL
  
  # trap
  stopifnot(is.ts(tsY))
  stopifnot(is.ts(tsX))
  stopifnot(length(tsX) == length(tsY))
  
  cat("plotting", sTitle, "...\n")
  plot(
    tsY, 
    tsX,
    type = "p",
    pch  = 3,
    xlab = NA, 
    ylab = NA, 
    main = sTitle
  )
  
  oOut.lm <- lm(formula)
  abline(oOut.lm, col="red")
  
  legend(
    "top", 
    legend = asLegend,
    lty    = c(-1, 1),
    pch    = c(3, -1), 
    col    = c("black", "red"),
    cex    = .8
  )
}
# - - - - - - - - - - - - -- -
sub.plotFigure_twoseries <- function(tsData1, tsData2, sTitle, asLegend){
  # purpose: plot time series + smoothed series (e.g. Figure 2.1)
  # args:    tsData1:  (ts object) observed time series (black line)
  #          tsData2:  (ts object) smoothed time series (red line)
  #          sTitle:   (string scalar) title
  #          asLegend: (string vector) legends of the lines
  # return:  NULL
  # notes:   Either tsData1 or tsData2 can be NULL.
  #          - When both tsData1 and tsData2 are specified, 
  #            asLegend should be c(label of tsData1, label of tsData2)
  #          - When only tsData1 are specified, 
  #            asLegend should be the label of tsData1
  #          - When only tsData2 are specified, 
  #            asLegend should be the label of tsData2
  #          time(tsData1) can be different with time(tsData2)
  
  # trap
  stopifnot(!is.null(tsData1) | !is.ts(tsData1))
  stopifnot(!is.null(tsData2) | !is.ts(tsData2))
  abExist <- c(
    !is.null(tsData1), 
    !is.null(tsData2)
  )
  stopifnot(any(abExist))
  stopifnot(sum(abExist) == length(asLegend))
  
  cat("plotting", sTitle, "...\n")
  # set frame
  plot(
    if (abExist[1]) { tsData1 } else { tsData2 }, 
    type = "n" ,
    xlab = NA, 
    ylab = NA,
    main = sTitle
  )
  # plot 1
  if (abExist[1]) {
    lines(tsData1, col = "black")
  }
  # plot 2
  if (abExist[2]) {
    lines(tsData2, col = "red")
  }
  
  legend(
    "top", 
    legend = asLegend,
    lty    = 1, 
    col    = c("black", "red"),
    cex    = .8
  )
  
}
# - - - - - - - - - - - - -- -
sub.plotFigure_interval <- function(tsData1, tsData2, tsWidth, sTitle, asLegend){
  # purpose: plot time series + smoothed series + interval (e.g. Figure 8.2)
  # args:    tsData1:  (ts object) observed time series (black line)
  #          tsData2:  (ts object) smoothed time series (red line)
  #          tsWidth:  (ts object) interval width
  #          sTitle:   (string scalar) title
  #          asLegend: (string vector) legends of the lines
  # return:  NULL
  # notes:   tsData1 can be NULL.
  #          - When both tsData1 and tsData2 are specified, 
  #            asLegend should be c(label of tsData1, label of tsData2)
  #          - When only tsData2 are specified, 
  #            asLegend should be the label of tsData2
  #          time(tsData2) should be equal to time(tsWidth)
  
  # trap
  if (!is.null(tsData1)){
    stopifnot(is.ts(tsData1))
    stopifnot(length(asLegend) == 2)
  } else {
    stopifnot(length(asLegend) == 1)
  }
  stopifnot(is.ts(tsData2))
  stopifnot(is.ts(tsWidth))
  stopifnot(time(tsData2) == time(tsWidth))
  
  cat("plotting", sTitle, "...\n")
  # set frame
  if (!is.null(tsData1)){
    plot(
      cbind(tsData1, tsData2+tsWidth, tsData2-tsWidth),
      plot.type = "single",
      type = "n",
      xlab = NA, 
      ylab = NA,
      main = sTitle
    )
  } else {
    plot(
      cbind(tsData2+tsWidth, tsData2-tsWidth),
      plot.type = "single",
      type = "n",
      xlab = NA, 
      ylab = NA,
      main = sTitle
    )
  }  

  # plot 1
  if (!is.null(tsData1)){
    lines(tsData1, col = "black", lty=1)
  }
  
  # plot 2
  lines(tsData2, col = "red", lty=3)
  lines(tsData2+tsWidth, col = "red")
  lines(tsData2-tsWidth, col = "red")

  if (!is.null(tsData1)){
    legend(
      "top", 
      legend = asLegend,
      lty    = c(1,3), 
      col    = c("black", "red"), 
      cex    = .8
    )
  } else {
    legend(
      "top", 
      legend = asLegend,
      lty    = 1, 
      col    = "red", 
      cex    = .8
    )
  }
}
# - - - - - - - - - - - - -- -
sub.plotFigure_residual <- function(tsResidual, sTitle, sLegend){
  # purpose: plot residuals (e.g. Figure 1.3)
  # args:    tsResidual: (ts object)     residual time series
  #          sTitle:     (string scalar) title
  #          sLegend:    (string scalar) legend of the line
  # return:  NULL
  
  stopifnot(is.ts(tsResidual))
  
  cat("plotting", sTitle, "...\n")
  plot(
    tsResidual,
    type  = "l",
    lty   = 2,
    col   = "green",
    xlab  = NA, 
    ylab  = NA,
    main  = sTitle
  )
  
  abline(h=0, col="gray")
  
  legend(
    "top", 
    lty    = 2, 
    legend = sLegend,
    col    = "green",
    cex    = .8
  )

}
# - - - - - - - - - - - - -- -
sub.plotFigure_auxresidual <- function(tsAuxResidual, sTitle, sLegend){
  # purpose: plot auxiliary residuals (e.g. Figure 8.11)
  # args:    tsAuxResidual: (ts object)     auxiliary residual time series
  #          sTitle:        (string scalar) title
  #          sLegend:       (string scalar) legend of the line
  # return:  NULL
  
  stopifnot(is.ts(tsAuxResidual))
  
  cat("plotting", sTitle, "...\n")
  plot(
    tsAuxResidual,
    type  = "l",
    lty   = 2,
    col   = "green",
    xlab  = NA, 
    ylab  = NA,
    main  = sTitle
  )
  
  abline(h=0, col="gray")
  abline(h=-1.96, col="red")
  abline(h=+1.96, col="red")
  
  legend(
    "top", 
    lty    = 2, 
    legend = sLegend,
    col    = "green",
    cex    = .8
  )
  
}
# - - - - - - - - - - - - -- -
sub.plotFigure_ACF <- function(tsData, nMaxLag, sTitle, sLegend){
  # purpose: plot a correlogram 
  # args:    tsData:     (ts object)      time series
  #          nMaxLag:    (integer scalar) max lag
  #          sTitle:     (string scalar)  title
  #          sLegend:    (string scalar)  legend of the line
  # return:  NULL
  
  stopifnot(is.ts(tsData))
  
  cat("plotting", sTitle, "...\n")
  acf(
    ts(tsData, frequency=1), 
    lag.max = nMaxLag,
    lwd     = 12, 
    col     = "green", 
    main    = sTitle,
    xlab    = NA, 
    ylab    = NA , 
    xlim    = c(1,nMaxLag)
  )
  
  legend(
    "top", 
    legend = sLegend, 
    lty    = 1, 
    col    = "green",
    cex    = .8
  ) 
  
}
