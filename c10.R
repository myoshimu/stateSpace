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
