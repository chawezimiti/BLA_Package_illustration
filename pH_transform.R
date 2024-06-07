# Step 1: Install and Load Necessary Packages

library(BLA)
library(MASS)
library(bestNormalize)

# Step 2: Load the BLA Package and Access the Soil Dataset
data("soil", package = "BLA")

# Extract the pH data and ensure it is numeric
ph_data <- as.numeric(soil$pH)
print(ph_data)
summastat(ph_data)

# Step 3: Apply the Box-Cox Transformation
# Check for non-negative values
if (min(ph_data) <= 0) {
  constant <- abs(min(ph_data)) + 1
  ph_data_bc <- ph_data + constant
} else {
  constant <- 0
  ph_data_bc <- ph_data
}

# Ensure ph_data_bc is numeric
ph_data_bc <- as.numeric(ph_data_bc)
print(ph_data_bc)

# Fit a linear model (using a dummy variable)
model_bc <- lm(ph_data_bc ~ 1)

# Verify that the fitted model's response is numeric
response <- model_bc$model[[1]]
print(response)
print(is.numeric(response))

# Simplified approach: Use boxcox function directly on the numeric vector
lambda <- MASS::boxcox(model_bc, plotit = TRUE)$x[which.max(boxcox(model_bc, plotit = FALSE)$y)]
print(lambda)

# Define the Box-Cox transformation function
box_cox_transform <- function(data, lambda) {
  if (lambda == 0) {
    log(data)
  } else {
    (data^lambda - 1) / lambda
  }
}

# Transform the data
transformed_bc <- box_cox_transform(ph_data_bc, lambda)

# If a constant was added, subtract it after transformation
if (constant > 0) {
  transformed_bc <- transformed_bc - constant
}

# Ensure transformed_bc is numeric
transformed_bc <- as.numeric(transformed_bc)
hist(transformed_bc)

# Step 4: Apply the Yeo-Johnson Transformation
# Apply the Yeo-Johnson transformation using bestNormalize
library(bestNormalize)
yj <- yeojohnson(ph_data)
transformed_yj <- yj$x.t

# Ensure transformed_yj is numeric
transformed_yj <- as.numeric(transformed_yj)

# Print transformed data
hist(transformed_yj)


# Step 5: Visualize the Data
# Visualize original data
hist(ph_data, main = "Original Data", xlab = "pH", breaks = 10)

# Visualize Box-Cox transformed data
hist(transformed_bc, main = "Box-Cox Transformed Data", xlab = "Transformed pH", breaks = 10)

# Visualize Yeo-Johnson transformed data
hist(transformed_yj, main = "Yeo-Johnson Transformed Data", xlab = "Transformed pH", breaks = 10)


summastat(transformed_yj)

x<-transformed_yj
y<-soil$yield
plot(x,y)
dat<-data.frame(x,y)
bag<-bagplot(dat)

dat1<-rbind(bag$pxy.bag,bag$pxy.outer)
x<-dat1[,1]
y<-dat1[,2]
plot(x,y)

vals<-data.frame(x,y)
vals<-vals[order(vals$x),]
startValues("lp")
start<-c(18,4,13.5, 0.0025, 9.27, 0.99,1.75, 0.118)
model<-cbvn(vals, start = start, model="lp", sigh = 0.7)
plot(x,y, col="grey", pch=16)

y2<- predictBL(model,vals$x)
lines(vals$x,y2, col="red", lwd=1.5)

#backtransform of x

x2 <- predict(yj, newdata = vals$x, inverse = TRUE)
plot(x2, vals$y, col="grey", pch=16)




lines(x2,y2, col="red", lwd=1.5)



