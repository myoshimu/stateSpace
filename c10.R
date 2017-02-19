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

