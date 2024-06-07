hist(soil$P, main = "", xlab = expression("soil P/ mg l"^-1), ylim = c(0,2500)) # From results, P can not be assumed to be from a normal distribution

hist(log(soil$P), main = "", xlab = expression("soil P/ log (mg l"^-1*")"))


hist(soil$pH, main = "", xlab = expression("soil pH"), ylim=c(0,1200)) # From results, P can not be assumed to be from a normal distribution

hist(pH_yj, main = "", xlab = expression("soil pH/ Yeo-Johnson transformed"*", " *lambda == 2 ))
