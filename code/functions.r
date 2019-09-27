# return a p-value testing if a linear regression slope = 1? 
slope_eq_1 <- function(lm){
  
  slope = exp(coef(summary(lm))[2, 1])
  se_slope = exp(coef(summary(lm))[2, 2])
  df = summary(lm)$df[2]
  
  pt(-abs((slope - 1) / se_slope), df = df) * 2
  
}
