# # # # # # # # # # # # # # # # # # #
#          main routines            #
# # # # # # # # # # # # # # # # # # #
# - - - - - - - - - - - - -- -
exercise.Chapter1 <- function(){
  
  # constants - - - - - -
  sFILENAME <- "UKdriversKSI.txt"
  # - - - - - - - - - - - 
  
  # reading data
  dfData <- read.table(sFILENAME, skip=1)
  tsData <- ts(log(dfData[,1]))
  
  # linear regression
  cat("Linear regression: \n")
  oOut.lm <- lm(tsData ~ time(tsData))
  print(oOut.lm)
  
  # plot
  sub.plotFigure_scatterreg(
    tsData, 
    tsData ~ time(tsData),
    "Figure 1.1", 
    c("log UK drivers KSI against time (in months)", "regression line")
  )
  sub.plotFigure_twoseries(
    tsData, 
    NULL, 
    "Figure 1.2", 
    "log UK drivers KSI"
  )
  sub.plotFigure_residual(
    tsData - predict(oOut.lm), 
    "Figure 1.3",
    "residuals"
  )
  sub.plotFigure_ACF(
    as.ts(rnorm(length(tsData))), 
    15, 
    "Figure 1.4", 
    "ACF-random residuals"
  )
  sub.plotFigure_ACF(
    tsData - predict(oOut.lm), 
    15, 
    "Figure 1.5", 
    "ACF-regression residuals"
  )
  
}
