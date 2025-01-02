
## Not part of script###
library(devtools)
install_github("chawezimiti/BLA") # install the development version of BLA from Github
##########################################################################################


################################# BLA R Package ##########################################

# This is a script to reproduce the outputs in the paper "BLA: An R package for carryout
# boundary line analysis for biological response."

# i). Installation and Loading Necessary Packages-----------------------------------------

# This illustration will use the BLA R package and some other accompanying packages. These 
# need to be installed and loaded in R.


required_packages<-c("BLA","MASS","bestNormalize","aplpack","Matrix","HDInterval")#required packages


install_if_missing <- function(packages) {# Function to install packages not already installed
  missing_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(missing_packages) > 0) {
    install.packages(missing_packages)
  }
}


install_if_missing(required_packages) # Call the function with the list of required packages

# ii). Load the necessary packages--------------------------------------------------------

library(BLA)
library(MASS)
library(bestNormalize) # For data transformation i.e Yeo_Johnson method
library(aplpack)       # For determining bivariate outliers
library(HDInterval)    # For determining the High Density Interval

### 1. INTRODUCTION ----------------------------------------------------------------------

# The boundary line methodology has been recognized as means to derive standard values 
# from which fertilizer recommendations can be made and also to determine from potential 
# limiting factors, the most limiting factor to production. In this illustration, we apply 
# the boundary line methodology to a dataset collected by AgSpace Agriculture Ltd to extract 
# information on the possible limitations to increased productivity and provide target 
# values that maximize production in order to increase productivity for soils of England.
# 
# The data set used for this session is called "soil" which consists of wheat yield, 
# soil phosphorus (P) and soil pH which were collected from sites across the England 
# as part of a Survey by AgSpace Agriculture Ltd in 2016. Boundary line models will be 
# fitted to soil P and soil pH data for the purpose of yield gap analysis.


### 2. BOUNDARY LINE FITTING FOR SOIL P --------------------------------------------------

# The fitting of the boundary line involves 5 steps. These include (1) data distribution 
# checks, (2) outlier detection, (3) test for evidence of boundary structures and 
# (5) fitting of the boundary line model. In this section, We apply these steps to the 
# soil P data.

# 2.1. Normality check for the soil P-----------------------------------------------------

# When the censored bivariate normal model is used, the variables x and y should plausibly 
# be drawn from an underlying bivariate normal distribution albeit with some censoring of 
# the y variable. Here we check this assumption using the summastat() function.
#

# 2.1.1. Soil P

summastat(soil$P) # From results, P can not be assumed to be from a normal distribution and
                  # so we try a log transformation.

summastat(log(soil$P))# The log-transformed P can be assumed to be from a normal distribution

# 2.1.2 Wheat yield

summastat(soil$yield)# Wheat yield can be assumed to be from a normal distribution

# 2.2. Outlier detection using bagplot() function ----------------------------------------

# The boundary line model is highly sensitive to outliers and for this reason, bivariate 
# outliers are identified and removed from the data. This is achieved using the bagplot 
# method. The bag plot is a bivariate equivalent of the univariate boxplot. It is composed 
# of a depth median, bag and loop. The depth median describes the center of the data cloud 
# and is equivalent to the median value in the univariate boxplot, the bag contains 50% of 
# the data and the loop contains data outside the bag which are not outliers. Here we identify 
# and remove outliers from the dataset using the bagplot() function from the aplpack package.

dat<-data.frame(x=log(soil$P), y=soil$yield) #Input for the bagplot() is a dataframe of x and y. 

out<-bagplot(dat,show.whiskers = F)# Figure 2 in article

legend("bottomright", legend = c("Bag","Loop","Depth median", "outlier"), #adds legend
       pch = c(15,15,8,3), col = c(adjustcolor( "#7799ff", alpha.f = 0.7),
               "#aaccff","red","red"))

# Creating a new data set without bivariate outliers i.e points in the bag and loop only.

vals<-rbind(out$pxy.bag,out$pxy.outer) 

head(vals)

# 2.3. Testing evidence for presence of boundary in dataset---------------------------------

# Fitting boundary line models to data works on the assumption that data has boundary structure
# at the upper edges of the data. This assumption can be assessed using expl_boundary() function.
# The inputs are the x and y variables.

x<-vals[,1]
y<-vals[,2]

expl_boundary(x,y,shells=10,simulations=1000) 

# Results indicate that there is strong evidence of bounding structures in both the right 
# and left sections.

# 2.4. Fitting boundary line model----------------------------------------------------------

# Before fitting the boundary line model, the data is plotted to get a feel of the distribution 
# of points at the upper edges. This allows the selection of an appropriate boundary model. This 
# should be supported by biological/agronomic plausibility.

plot(x,y,pch=16, col="grey", xlab="soil P", ylab = "Yield (t/ha)") 

# A trapezium boundary model is appropriate for this data. The censored bivariate normal
# procedure (cbvn) will be used to fit the boundary model. To fit the cbvn, some arguments
# need to be estimated. These include the initial starting values of the trapezium model
# and value for the standard deviation of measurement error.

## (a) Determination of initial starting values 

# Here we come up with the initial starting values for the optimization process.
# Since the selected model is the trapezium, the starting values will be a vector of length 
# 10 comprising the y-intercept, slope, plateau value, second y-intercept, second slope and 
# data distribution properties i.e  means of x and y, standard deviation of x and y, and 
# their correlation.

# To determine start values for the selected trapezium model, the startValues() function 
# is used. We determine several start values to avoid sticking to a local optimum solution.
# Note that for the startValues() function to work effectively, the appearance zoom in the 
# global options of R and your computer display zoom should be equal i.e.both should be at 
# 100%.

startValues("trapezium") 

# click on the plot 4 points that make up the trapezium structure at the upper edge of
# the data cloud in a clockwise movement. We run this step multiple times and the values 
# are added to the list of starting values below.

start<-list(# Table 3 in article
  c(4.3,3.4,13.8,32.8,-4.9,mean(x),mean(y),sd(x),sd(y),cor(x,y)),
  c(2.5,4.1,13.42,32,-4.8,mean(x),mean(y),sd(x),sd(y),cor(x,y)),
  c(3.5,3.7,13.35,47.7,-8.4,mean(x),mean(y),sd(x),sd(y),cor(x,y)),
  c(2.83,4.11,13.7,32.,-4.6,mean(x),mean(y),sd(x),sd(y),cor(x,y)),
  c(4.1,3.4,13.6,29,-4.1,mean(x),mean(y),sd(x),sd(y),cor(x,y))
)

## (b) Determination of standard deviation of measurement error

# The standard deviation of measurement error (sigh) is a fixed parameter in the boundary
# line determination using cbvn. As we don't have a direct measure of sigh, it is determined by 
# log-likelihood profiling. In this process, the negative log-likelihoods of several suggested 
# values of sigh are determined given the data distribution and the suggested model. The sigh value 
# with the smallest negative log-likelihood is selected. This is done using the ble_profile() 
# function at the multiple start values. 

sigh<-c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1) # we suggest sigh values ranging from 0.3 to 1 t/ha

# Note: This function will take a few minutes to run as it is evaluating each sigh value
# for each set of start values. We use the lapply function to loop over each sigh and start 
# value set.

profile <- lapply(start, function(start) {#estimates the log-likelihood ~ sigh for each set of start values
  
  lapply(sigh, function(sigh) {
    
    ble_profile(data = vals, start = start, sigh = sigh, model = "trapezium", plot = F)
  })
})

# Extract log-likelihood and sigh (Merror) values from profile list

log_likelihood <- unlist(lapply(profile, function(sublist) {
  sapply(sublist, function(x) x[["log-likelihood"]])
}))

Merror<- unlist(lapply(profile, function(sublist) {
  sapply(sublist, function(x) x[["Merror"]])
}))

# The sigh with the smallest negative log-likelihood 

Merror[which(log_likelihood==max(log_likelihood, na.rm = T))]# profile maximized  at 0.4 t/ha


# Plotting the log-likelihood profile for visualization

me_profile<-data.frame(x=Merror,y=log_likelihood)
me_profile<-me_profile[order(Merror),] # ordering data from smallest to largest Merror

vec<-vector()

for(i in unique(me_profile$x)){ # getting largest log-likelihood at each Merror value
  maxy<-max(me_profile[which(me_profile$x==i),]$y, na.rm = T)
  vec<-c(vec,maxy)
}

# ploting the largest value at each Merror

par(mar=c(5,5,4,4))

plot(unique(me_profile$x),vec,pch=16, # Figure 3 in article
     xlab=expression(bold(sigma[e]*"/t ha"^-1)),
     ylab=expression(bold(log-Likelihood)),
     cex.lab=1.8, cex.axis=1.8)
lines(unique(me_profile$x),vec, lty=5, lwd=1.5) 
abline(v=Merror[which(log_likelihood==max(log_likelihood, na.rm = T))],
       lty=5, col="red", lwd=1.3) # sigh with largest log-likelihood

## (c) Fitting the boundary line using censored bivariate normal model 

# All the arguments for the censored bivariate normal model function are set. To fit the 
# boundary line, the function cbvn() is used. We do this using multiple start values to 
# avoid sticking at a local optimum solution. The solution with the smallest negative 
# log-likelihood value is selected. This is equivalent to selecting the model with the 
# smallest Aikaike Information Criterion (AIC) value, which is an output of the cbvn() 
# function. The lapply() function is used here to allows us to loop the cbvn() over the 
# different start values. The tryCatch() allows the loop function to continue to the next 
# input when the function fails to converge at a given input combination.

models <- lapply(start, function(start) {
  
  tryCatch(
    
    cbvn(data=vals, start = start, sigh=0.4, model = "trapezium",
         optim.method = "Nelder-Mead",
         xlab=expression("Phosphorus/ mg L"^-1), 
         ylab=expression("Yield/ t ha"^-1), pch=16, col="grey"),
    
    error=function(e) NA) 
  })

## Selects the solution with the smallest value of AIC

model_1 <- models[[which.min(unlist(lapply(X=models,FUN = function(a){
  
  b<-tryCatch(a$AIC[2,1],error=function(e) NA)
  return(b)})))]]

model_1


# 2.5. Predicting boundary yield for each data point in our dataset-------------------------

# The largest expected yield for each value of soil P in our data set is estimated
# from the boundary line model parameters. This is done using the function predictBL(). 
# The inputs are the model object and a vector of soil P values.

P_data <-log(soil$P) # extracting soil P from the data set

P<-predictBL(model_1,P_data) # boundary yield for soil P


# 2.6. Post-hoc analysis for soil P boundary line model-------------------------------------

## (a) Determining the Critical soil P concentration 

# This is calculated as the quotient of the difference between the plateau yield and the y-intercept,
# and the slope.

crit_P<-(model_1$Parameters[3,1]-model_1$Parameters[1,1])/model_1$Parameters[2,1]

crit_P #Critical soil P on log-scale

exp(crit_P) #Critical soil P on the original scale

# (b) Determining uncertainty around critical P as 95% Confidence Interval.

# This is done by generating several model parameters using the mvrnorm() with means equal
# to the determined model parameters and sigma equal to the variance-covariance of the 
# parameters. The variance-covariance is estimated from the Hessian matrix (output of the 
# cbvn function). The critical value is estimated for each set of parameters and the 95% CI 
# is determined from these values.

params <- mvrnorm(n = 100000,mu = model_1$Parameters[1:3,1], Sigma = solve(model_1$Hessian[1:3,1:3]))

crit_points<-(params[,3]-params[,1])/params[,2] # calculate critical points for simulated parameters

CI_95P<-hdi(crit_points, credMass = 0.95) # 95% High Density Interval

CI_95P # 95% Confidence interval of critical soil P concentration on log scale
exp(CI_95P) # 95% Confidence interval of Critical soil P concentration on the original scale

# (c) Probability that soil P is limiting at each point

# This is done by defining a cumulative probability function of the variable being less 
# than critical value given the confidence interval around.

P_limProb<-1-pnorm(P_data, mean = mean(crit_points), sd = sd(crit_points))

head(P_limProb)

# (d) Plot the boundary yields as function of soil P on the log scale
par(mar=c(5,5,4,4))
plot(vals[,1],vals[,2], xlab=expression(bold("Phosphorus/ ln(mg L"^-1*")")), # plot the data 
     ylab=expression(bold("Yield/ t ha"^-1)), pch=16, col="16", cex.axis=1.8,cex.lab=1.8)

# Generating soil P concentration values to predict boundary yield in the range of our data

xa<-seq(min(vals[,1]),max(vals[,1]),length.out =1000)

ya<-numeric() # empty vector to contain predicted boundary yield values

for(i in 1: length(xa)){ # predicting boundary yield using the determined model parameters
  ya[i]<-min(model_1$Parameters[1,1]+model_1$Parameters[2,1]*xa[i]
             ,model_1$Parameters[3,1],
             model_1$Parameters[4,1]+model_1$Parameters[5,1]*xa[i])
  
}

lines(xa,ya, col="red", lwd=2) # adding the boundary yield for each soil P value
abline(v=crit_P, col="black", lty=5, lwd=2) # adding critical soil P to plot
abline(v=CI_95P, col="blue", lty=1, lwd=2) # adding confidence interval around soil P 


# (e) Plot data on the original scale of soil P concentration in mg/l (Figure 4a in article)

plot(exp(vals[,1]),vals[,2], xlab=expression(bold("Phosphorus/ mg L"^-1)),
     ylab=expression(bold("Yield/ t ha"^-1)), col="grey", pch=16,
     cex.axis=1.8,cex.lab=1.8)
lines(exp(xa),ya, col="red", lwd=2) # adding the boundary yield for each soil P value
abline(v=exp(crit_P), col="black", lty=5, lwd=2) # adding critical soil P to plot

CI_95P_original <- exp(CI_95P)

polygon(c(CI_95P_original[1],CI_95P_original[2],CI_95P_original[2],CI_95P_original[1],
          CI_95P_original[1]), c(0,0,18,18,0), col=adjustcolor( "red", alpha.f = 0.2),
        border=NA) # add the confidence interval around critical soil P

legend("bottomright", legend = c("Boundary line","Critical Soil P", "95% CI of critical Soil P"),
       lty = c(1,5,NA), pch = c(NA,NA,15), lwd = c(1.5,2, NA),
       col = c("red","black",col=adjustcolor( "red", alpha.f = 0.2) ), cex=1.5)



### 3. BOUNDARY LINE FITTING FOR SOIL pH --------------------------------------------------

# The boundary line model is fitted to soil pH following the same procedure outlined above.

# 3.1. Normality check for pH and wheat yield ----------------------------------------------

# As was done with soil P, we test the plausibly that soil pH is drawn from an underlying 
# bivariate normal distribution

# (a) Original data

summastat(soil$pH) 
# pH can not be assumed to be from normal distribution. We try transforming the data to log.

# (b) Log transformation

summastat(log(soil$pH)) 
# log transformed pH can not be assumed to be from normal distribution. We try the box-cox transformation

# (c) box transformation

## Fit a linear model (using a dummy variable)

model_bc <- lm(soil$pH ~ 1)

# Use box cox function directly on the numeric vector
bcox<- MASS::boxcox(model_bc, plotit = TRUE)
maxlik<-max(bcox$y)
lambda<-bcox$x[which(bcox$y==maxlik)]
print(lambda)

## Define the Box-Cox transformation function
box_cox_transform <- function(data, lambda) {
  if (lambda == 0) {
    log(data)
  } else {
    (data^lambda - 1) / lambda
  }
}

## Transform the data

pH_bc <- box_cox_transform(soil$pH, lambda)

summastat(pH_bc) 
# Box cox transformed pH can not be assumed to be from normal distribution. We try the 
# Yeo-Johnson transformation.

## (d) Yeo-Johnson Transformation

yj <- yeojohnson(soil$pH)

pH_yj <- yj$x.t

summastat(pH_yj)
# The Yeo-Johnson transformed pH can be assumed to be from normal distribution as O.skewness 
# is within [-0.2,0.2]


# 3.2. Outlier detection using bagplot() function ------------------------------------------

# The boundary line model is highly sensitive to outliers. For this reason, bivariate outliers
# are identified and removed from the data. This is achieved using the function bagplot()
# from the aplpack package. The bag plot is a bivariate equivalent of the univariate boxplot.
# it is composed of a bag and loop. The bag contains 50% of the data and the loop contains
# data outside the bag which are not outliers.

dat2<-data.frame(pH_yj,soil$yield)# Input for bagplot() is a dataframe of x and y.

out2<-bagplot(dat2, show.whiskers = F)
legend("bottomright", legend = c("Bag","Loop", "depth median","outlier"), pch = c(15,15,8,3),
       col = c(adjustcolor( "blue", alpha.f = 0.7),
               "lightblue","red","red" ))

# Creating a new data set without bivariate outliers i.e points in the bag and loop only.

vals2<-rbind(out2$pxy.bag,out2$pxy.outer) 

# 3.3. Fitting the censored boundary model to data -----------------------------------------

# The data can be plotted to get a feel of the distribution of points at the upper edges.
# This allows to selected an appropriate boundary model. This should be supported by 
# biological/agronomic plausibility.

x2<-vals2[,1]
y2<-vals2[,2]

plot(x2,y2, pch=16, col="grey", xlab="soil pH", ylab = "Yield (t/ha)") 

# A linear-plateau boundary model is appropriate for this data. The censored bivariate normal
# procedure (cbvn) will be used to fit the boundary model. To fit the cbvn, some arguments
# need to be determined. These include the initial starting values of the linear-plateau model
# and the measurement error standard deviation (already determined using soil P).

## (a) Determination of starting values for the linear-plateau model

# Since the selected model is the linear-plateau, the starting values will be a vector of 
# length 8 comprising the y-intercept, slope, plateau value and data distribution properties 
# i.e  means of x and y, standard deviation of x and y, and their correlation.

# To determine start values for the selected boundary model, the startValues() is used.
# We determine several start values to avoid selection of a local optimum solution.
# Note that for this functions to work effectively, the Appearance zoom in the global 
# options of R and your computer display zoom should be equal i.e.both should be at 
# 100%.


startValues("lp") 

# click on the plot 2 points that make up the linear-plateau structure at the upper edge of
# the data cloud in a clockwise movement. We do this Multiple times and the values are added
# to the list of starting values below.


start_pH<-list( # Table 3 in article
  c(20,4.3,13.57,mean(x2),mean(y2),sd(x2),sd(y2),cor(x2,y2)),
  c(21.3,5.4,13.44,mean(x2),mean(y2),sd(x2),sd(y2),cor(x2,y2)),
  c(20.9,5.2,13.57,mean(x2),mean(y2),sd(x2),sd(y2),cor(x2,y2)),
  c(18.47,3.54,13.87,mean(x2),mean(y2),sd(x2),sd(y2),cor(x2,y2)),
  c(20.7,4.5,13.67,mean(x2),mean(y2),sd(x2),sd(y2),cor(x2,y2))
)

## (b) Fitting the boundary line using censored bivariate normal model 

# All the arguments for the censored bivariate normal model function are set. To fit the 
# boundary line, the function cbvn() is used. We do this using multiple start values to 
# avoid local optima solution. The solution with the smallest negative log-likelihood value is 
# selected. This is equivalent to selecting the model with the smallest Aikaike Information 
# Criterion (AIC) value, which is an output of the cbvn() function. # The lapply() function 
# is used here to allows us to loop the cbvn() over the different start values. The tryCatch() 
# allows the loop function to continue to the next input when the function fails to converge 
# at a given input combination.


models2 <- lapply(start_pH, function(start_pH) {
  
  tryCatch(
    
    cbvn(vals2, start = start_pH, sigh=0.4, model = "lp",
         optim.method = "Nelder-Mead",
                     xlab=expression("pH/tranformed YeoJohnson"), 
                     ylab=expression("Yield/ t ha"^-1), 
                     pch=16, col="grey"), 
    error=function(e) NA) 
  })

## Select the solution with the smallest AIC

model_2 <- models2[[which.min(unlist(lapply(X=models2,FUN = function(a){
  b<-tryCatch(a$AIC[2,1],error=function(e) NA)
  return(b)})))]]

model_2

# 3.4. Predicting boundary yield for each data point ---------------------------------------

# The largest expected yield for each value of soil pH in our data set is estimated
# using the boundary line model parameters. This is done using the function predictBL(). 
# The inputs are the model object and a vector of soil P values.

pH_data<- yeojohnson(soil$pH)$x.t # extracting soil pH from data and transforming it

pH<-predictBL(model_2,pH_data) # predicting boundary yield


# 3.5. Post-hoc analysis for soil pH boundary line model------------------------------------


# (a) Determining critical soil pH value

crit_ph<-(model_2$Parameters[3,1]-model_2$Parameters[1,1])/model_2$Parameters[2,1]

crit_ph # critical soil pH on log-scale
predict(yj, newdata = crit_ph, inverse = TRUE) #critical soil pH on the original scale


# (b) Determining uncertainty around critical soil pH as 95% Confidence Interval.

# This is done by generation several model parameters using the mvrnorm() with means equal
# to the determined model parameters and sigma equal to the variance-covariance of the 
# parameters. The variance-covariance is estimated from the Hessian matrix (output of the 
# cbvn function). The critical value is estimated for each set of parameters and the 95% CI 
# is determined from these values.

params3 <- mvrnorm(n = 100000,mu = model_2$Parameters[1:3,1], Sigma = solve(model_2$Hessian[1:3,1:3]))

crit_points3<-(params3[,3]-params3[,1])/params3[,2] # generated critical points
CI_95ph<-quantile(crit_points3, c(0.025,0.975)) # 95% confidence interval

CI_95ph # confidence interval on Yeo-Johnson transformed scale
predict(yj, newdata = CI_95ph, inverse = TRUE) # confidence interval on the original scale


# (c) Probability that pH is limiting at each point.

# This is done by defining a cumulative probability function of the variable being less 
# than critical value given the confidence interval around it.

pH_limProb<-1-pnorm(pH_data, mean = mean(crit_points3), sd = sd(crit_points3))

head(pH_limProb)

# (d) Plot the boundary yields as function of soil pH on the Yeo-Johnson scale

plot(vals2[,1],vals2[,2], xlab="pH", ylab=expression("Yield/ t ha"^-1),# plot the data 
     col="grey", pch=16) 

# Generating soil pH values to predict boundary yield in the range of our data

xc<-seq(min(vals2[,1]), max(vals2[,1]),length.out =1000)

yc<-numeric() # empty vector for adding predicted boundary yield values from xc

for(i in 1: length(xc)){  # predicting boundary yield using the model parameters
  yc[i]<-min(model_2$Parameters[1,1]+model_2$Parameters[2,1]*xc[i]
             ,model_2$Parameters[3,1])
  
}

lines(xc,yc, col="red", lwd=2) # adding the boundary line
abline(v=crit_ph, col="black", lty=5, lwd=2) # adding critical soil pH to plot
abline(v=CI_95ph, col="blue", lty=1, lwd=2)# adding confidence interval around soil pH 


# (e) Plot data on the original scale of soil pH (Figure 4b in article)

x2_original<-predict(yj, newdata = x2, inverse = TRUE)

plot(x2_original,y2, xlab=expression(bold("pH")),
     ylab=expression(bold("Yield/ t ha"^-1)), col="grey", 
     pch=16, cex.axis=1.8,cex.lab=1.8)

xc2 <- predict(yj, newdata = xc, inverse = TRUE)# back-transforming the generated soil pH data
lines(xc2,yc, col="red", lwd=2) # adding the boundary yield for each soil pH value
abline(v=predict(yj, newdata = crit_ph, inverse = TRUE), col="black", lty=5, lwd=2) # adding critical soil pH to plot

CI_95_original <- predict(yj, newdata = quantile(crit_points3, c(0.025,0.975)), inverse = TRUE)# back-transforming CI to original scale

polygon(c(CI_95_original[1],CI_95_original[2],CI_95_original[2],CI_95_original[1],
          CI_95_original[1]), c(0,0,18,18,0), col=adjustcolor( "red", alpha.f = 0.2),
        border=NA) # add the confidence interval around critical soil P

legend("bottomright", legend = c("Boundary line","Critical Soil pH", "95% CI of critical pH"),
       lty = c(1,5,NA), pch = c(NA,NA,15), lwd = c(2,2, NA), cex = 1.8,
       col = c("red","black",col=adjustcolor( "red", alpha.f = 0.2) ))


### 4. DETERMING MOST LIMITING FACTOR-----------------------------------------------------

# The most limiting factor at each point in our dataset can be determined using the boundary line
# models of soil P and soil pH. The factor that predicts the lowest yield a given point is 
# considered as the most limiting factor. However, if both soil P and soil pH are greater
# than the critical value, the limiting factor is considered unknown. This is achieved 
# using the limfactor() function. The inputs are the predicted values of yield for each factor.

# 4.1. Extracting the most limiting factor

yieldlim0<-limfactor(P,pH)

# the output is a list of length 2 that includes (1) a dataframe of the identified most limiting 
# factor and the largest yield predicted by the most limiting factor at that point, and 
# (2) largest predicted yield from the dataset (attainable yield)

yieldlim<-yieldlim0[[1]] # extracts the dataframe of the identified most limiting factor and the largest predicted yield

yieldlim$P_limProb<-P_limProb # adds the probability of P being limiting at each point
yieldlim$pH_limProb<-pH_limProb # adds the probability of pH being limiting at each point

head(yieldlim)

# 4.2 Visualize the limiting factor proportions (Figure 5 in article)

props<- prop.table(table(yieldlim$Lim_factor))

par(mar=c(5,5,4,2))

barplot(props*100, ylab="Proportion identified (%)", xlab = "Limiting factor",
        ylim=c(0,100), cex.lab=1.8, cex.axis=1.5, cex=1.5)
abline(h = 0, col = "black", lwd = 1)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ END ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~










