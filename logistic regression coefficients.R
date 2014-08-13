# Code to perform a logistic regression from the raw occupancy simulations
# of a species




#Create some fake simulations
spMat <- matrix(NA,1000,6)
set.seed(23)
for (i in 1:1000){
  spMat[i,] <- plogis(2*rnorm(1) + 0.1*rnorm(1,1,.5)*(1:6))
}
  
  #Function to do this with one row of data (one simulation)
  logRegression <- function(row,n,years){
    #code to perform a logistic regression on each row of data and return coefficients
    # row is a vector with occupancy values for each year
    # n is the number of camera traps at a site
    # y is numeric vector enumerating years. Eg. 1:6 or 2007:2011
    y <- round(cbind(row*n,(1-row)*n))
    model<-glm(y~years,family="binomial")
    model$coeff[2]
  }

coeffs <- apply(spMat,1,logRegression,60,1:6) # might get some warnings. Its ok.

median(coeffs) # close to 0.1 the original value we simulated for the slope

#Transform to odds ratios

coeffs.OR <- plogis(coeffs)

#Report median and confidence limits for the odds ratio
# If odds ratio > 1 and confidence limits do not overlap 1, species is increasing
# If odds ratio < 1 nd confidence limits do not overlap 1, species is decreasing
# Otherwise species is stable

median <- median(coeffs.OR) # species seems to be decreasing
conflimits80 <- quantile(coeffs.OR,c(0.1,0.9))
conflimits95 <- quantile(coeffs.OR,c(0.025,0.975))

# Since conflimits do not overlap with 1 this species is decreasing
# With each year the probability of presence over no presence decreases by 54%
# median - value

# We want to report the median and the confidence limits for both cases (80 and 95)
# But decide on change when conflimits80 does not overlap 1