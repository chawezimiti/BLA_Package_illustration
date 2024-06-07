

## BOUNDARY LINE ANALYSIS USING CENSORED BiVARIATE NORMAL MODEL #############

# Step 1: Install and Load Necessary Packages

## List of packages required

packages <- c("BLA", "MASS", "bestNormalize", "aplpack", "Matrix")

## Function to check if a package is installed and install it if not

install_if_missing <- function(p) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p, dependencies = TRUE)
  }
}

## Load the packages

library(BLA)
library(MASS)
library(bestNormalize)
library(aplpack)

### Reading in data ----------------------------------------------------------------------

# The data set used for this session is called "soil" consisting of wheat yield, 
# soil phosphorus (P) and soil pH which were collected from sites across the UK 
# as part of a Survey. Boundary line methodology will be conducted on soil P and soil pH 
# data for yield gap analysis.


### Boundary line fitting for Soil P -----------------------------------------------------

# The fitting of the boundary line involves 5 steps. The include (1) normal distribution 
# checks, (2) outlier detection, (3) test for evidence of boundary presence and 
# (5) fitting of the boundary line.

# 1. Normality check for the P and yield data---------------------------------------------

summastat(soil$P) # From results, P can not be assumed to be from a normal distribution

summastat(log(soil$P)) # The log-transformed P can be assumed to be from a normal distribution

# 2. Outlier detection using bagplot() function ------------------------------------------

dat<-data.frame(x=log(soil$P), y=soil$yield) #Input for bagplot is a dataframe of x and y. 

bag<-bagplot(dat, show.whiskers = F)

vals<-rbind(bag$pxy.bag,bag$pxy.outer) # new data set without bivariate outliers

# 3. Testing evidence for presence of boundary in dataset---------------------------------

x<-vals[,1]
y<-vals[,2]

expl.boundary(x,y,10,1000) # there is evidence of bounding structures

# 4. Fitting boundary line model----------------------------------------------------------

plot(x,y) # a trapezium boundary model is appropriate

## (a) determination of starting values for the trapezium model

startValues("trapezium")

start<-c(3.8,3.6,13.5,35,-5,mean(x),mean(y),sd(x),sd(y),cor(x,y))

## (b) determination of standard deviation of measurement error

sigh<-c(0.45,0.55,0.6,0.65,0.7)

ble_profile(data=vals, start=start, sigh=sigh, model = "trapezium")

## (c) fitting the boundary line using censored bivariate normal model 

model_1 <-cbvn(data=vals, start = start, sigh=0.6, model = "trapezium",
               optim.method =  "BFGS",
               xlab=expression("Phosphorus/ mg L"^-1), 
               ylab=expression("Yield/ t ha"^-1), 
               pch=16, col="grey")

model_1


# 5. Predicting boundary yield for each data point --------------------------------------

P_data <-log(soil$P) # extracting soil P from the data

P_data[which(is.na(P_data)==T)]<-mean(P_data,na.rm=T) # replacing missing values with the mean value

P<-predictBL(model_1,P_data) # boundary yield for soil P

points(P_data, P, col="red", pch=16)



### Boundary line fitting for Soil PH -----------------------------------------------------


# 1. Normality check for pH and wheat yield ----------------------------------------------

# (a) Original data

summastat(soil$pH) # can not be assumed to be from normal distribution

# (b) Log transformation

summastat(log(soil$pH)) # can not be assumed to be from normal distribution

# (c) box transformation

## Fit a linear model (using a dummy variable)

model_bc <- lm(soil$pH ~ 1)

# Simplified approach: Use boxcox function directly on the numeric vector
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

summastat(pH_bc) # can not be assumed to be from normal distribution

## (d) Yeo-Johnson Transformation

yj <- yeojohnson(soil$pH)

pH_yj <- yj$x.t

summastat(pH_yj) # can be assumed to be from normal distribution


# 2. Outlier detection using bagplot() function ------------------------------------------

dat2<-data.frame(pH_yj,soil$yield)

bag<-bagplot(dat2, show.whiskers = F)

vals2<-rbind(bag$pxy.bag,bag$pxy.outer) # new excludes bivariate outliers

# 3. Exploring presence of boundary using the clustering methodology --------------------

# This has already been done using soil P


# 4. Fitting the censored boundary model to data -----------------------------------------

plot(vals2[,1],vals2[,2]) # a linear model is appropriate

## (a) determination of starting values for the trapezium model

startValues("lp") 

mean(vals2[,1])
mean(vals2[,2])
sd(vals2[,1])
sd(vals2[,2])
cor(vals2[,1],vals2[,2])

start<-c(20,4.6,13.3, 0.003, 9.27, 0.99,1.76, 0.118)

## (b) fitting the boundary line using censored bivariate normal model 

model_2<-cbvn(vals2, start = start, sigh=0.71, model = "lp",
              xlab=expression("pH/tranformed YeoJohnson"), 
              ylab=expression("Yield/ t ha"^-1), 
              pch=16, col="grey")


model_2

# 5. Predicting boundary yield for each data point ---------------------------------------

pH_data<- yeojohnson(soil$pH)$x.t # extracting soil pH from data and transforming it

pH_data[which(is.na(pH_data)==T)]<-mean(pH_data,na.rm=T)#replacing missing values with mean

pH<-predictBL(model_2,pH_data) # predicting boundary yield

points(pH_data, pH, col="red", pch=16)


# Post-hoc analysis ------------------------------------------------------------------

## 1. Critical values and uncertainty-----------------------------------------------------

### soil P --------------------------------------------------------------------------------

# (a) Plot of yield against soil P

plot(vals[,1],vals[,2], xlab=expression("Phosphorus/ln(mg kg"^-1*")"), 
     ylab="yield/ t ha", pch=16, col="16")

xa<-seq(min(vals[,1]),max(vals[,1]),length.out =1000) #values of soil P to predict yield

ya<-numeric() # empty vector to contain predicted boundary yield values

for(i in 1: length(xa)){ # predicting boundary yield values
  ya[i]<-min(model_1$Parameters[1,1]+model_1$Parameters[2,1]*xa[i]
          ,model_1$Parameters[3,1],
          model_1$Parameters[4,1]+model_1$Parameters[5,1]*xa[i])
  
}

lines(xa,ya, col="red", lwd=1.5) # adding the boundary values for each soil P value


# (b) Critical soil P value 

crit_P<-(model_1$Parameters[3,1]-model_1$Parameters[1,1])/model_1$Parameters[2,1]

abline(v=crit_P, col="red", lty=5) # adding critical soil P to plot

# (c) calculating uncertainty around critical P as 95% Confidence Interval

params <- mvrnorm(n = 100000,mu = model_1$Parameters[1:3,1], Sigma = solve(model_1$Hessian[1:3,1:3]))
crit_points<-(params[,3]-params[,1])/params[,2]
CI_95 <- quantile(crit_points, c(0.025,0.975))

abline(v=CI_95, col="blue", lty=1, lwd=1.5)

# (d) Probability that soil P is not limiting at each point

P_limProb<-pnorm(P_data, mean = mean(crit_points), sd = sd(crit_points))

head(P_limProb)

# (e) Back transform to the original scale in mg/l

plot(exp(vals[,1]),vals[,2], xlab=expression(bold("Phosphorus/ mg kg"^-1)),
     ylab=expression(bold("Yield/ t ha"^-1)), col="grey", pch=16, xlim=c(0,120),
     font.axis = 2)
lines(exp(xa),ya, col="red", lwd=1.5)
abline(v=exp(crit_P), col="red")

CI_95 <- exp(quantile(crit_points, c(0.025,0.975)))

polygon(c(12.2,19.0,19.0,12.2,12.2), c(0,0,18,18,0),
        col=adjustcolor( "red", alpha.f = 0.2),border=NA) # the confidence interval 

legend("bottomright", legend = c("Boundary line", "95% CI"),
       lty = c(1,NA), pch = c(NA,15), lwd = c(1.5, NA),
       col = c("red",col=adjustcolor( "red", alpha.f = 0.2) ))

### Soil pH ------------------------------------------------------------------------------

# (a) plot of soil pH vs yield

plot(vals2[,1],vals2[,2], xlab="pH", ylab=expression("Yield/ t ha"^-1),
     col="grey", pch=16) # ploting soil pH vs yield

xc<-seq(min(vals2[,1]), max(vals2[,1]),length.out =1000)# values of soil pH to be predicted

yc<-numeric() # vector for predicted boundary yield values

for(i in 1: length(xc)){ # predicting boundary yield values
  yc[i]<-min(model_2$Parameters[1,1]+model_2$Parameters[2,1]*xc[i]
             ,model_2$Parameters[3,1])
  
}

lines(xc,yc, col="red", lwd=1.5) # adding the boundary yield values to plot


# (b) Determining critical soil pH value and uncertainty

crit_ph<-(model_2$Parameters[3,1]-model_2$Parameters[1,1])/model_2$Parameters[2,1]

abline(v=crit_ph, col="red", lty=5) # adding to plot

# (c) Calculating uncertainty around critical pH

params3 <- mvrnorm(n = 100000,mu = model_2$Parameters[1:3,1], Sigma = solve(model_2$Hessian[1:3,1:3]))

# Since 'Sigma' is not positive definite, we convert it to a near positive matrix. 
# Another approach is to use Cholesky decomposition with a small positive value added 
# to the diagonal until the matrix is positive definite. This method ensures that the 
# adjustments are minimal.

Hessian<- model_2$Hessian[1:3,1:3]

# Define a small initial value for epsilon

epsilon <- 1e-10

# Function to check if a matrix is positive definite

is_positive_definite <- function(m) {
  all(eigen(m)$values > 0)
}

# Increment epsilon until the matrix is positive definite

Hessian_pd <- Hessian
while (!is_positive_definite(Hessian_pd)) {
  Hessian_pd <- Hessian + diag(epsilon, nrow(Hessian))
  epsilon <- epsilon * 5
}

# Compute the covariance matrix as the inverse of the positive definite matrix

cov_matrix <- solve(Hessian_pd)

params3 <- mvrnorm(n = 100000,mu = model_2$Parameters[1:3,1], Sigma = cov_matrix)#

crit_points3<-(params3[,3]-params3[,1])/params3[,2] # critical points

CI_95ph<-quantile(crit_points3, c(0.025,0.975)) # 95% confidence interval

abline(v=CI_95ph, col="blue", lty=1, lwd=1.5)

polygon(c(CI_95ph[1],CI_95ph[2],CI_95ph[2],CI_95ph[1],CI_95ph[1]), c(0,0,18,18,0),
        col=adjustcolor( "red", alpha.f = 0.2),border=NA)

# (d) Probability that pH is not limiting at each point

pH_limProb<-pnorm(pH_data, mean = mean(crit_points3), sd = sd(crit_points3))

head(pH_limProb)

# (e) Back transform the plot of yield~pH to the original scale in pH

x2 <- predict(yj, newdata = vals2[,1], inverse = TRUE)# back-transforming the pH data

plot(x2,vals2[,2], xlab=expression(bold("pH")),
     ylab=expression(bold("Yield/ t ha"^-1)), col="grey", 
     pch=16)

xc2 <- predict(yj, newdata = xc, inverse = TRUE)# back-transforming the pH data
lines(xc2,yc, col="red", lwd=1.5) # boundary line

abline(v=predict(yj, newdata = crit_ph, inverse = TRUE), col="red") # critical pH

CI_95 <- predict(yj, newdata = quantile(crit_points3, c(0.025,0.975)), inverse = TRUE)# back-transforming CI

polygon(c(CI_95[1],CI_95[2],CI_95[2],CI_95[1],CI_95[1]), c(0,0,18,18,0),
        col=adjustcolor( "red", alpha.f = 0.2),border=NA) # the confidence interval on plot

legend("bottomright", legend = c("Boundary line", "95% CI"),
       lty = c(1,NA), pch = c(NA,15), lwd = c(1.5, NA),
       col = c("red",col=adjustcolor( "red", alpha.f = 0.2) ))

### 2. DETERMING MOST LIMITING FACTOR--------------------------------------------------------------------


yieldlim0<-limfactor(P,pH)

yieldlim<-yieldlim0[[1]]

yieldlim$P_limProb<-P_limProb
yieldlim$pH_limProb<-pH_limProb

head(yieldlim)

# (b) Visualize the limiting factor proportions

props<- prop.table(table(yieldlim$Lim_factor))
pie(props, main="",
    labels = paste(names(props), "\n",scales::percent(as.vector(unname(props)))), 
    cex = 0.8,col = gray(seq(0.4, 1.0, length.out = 6)))



