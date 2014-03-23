# Code to fit two simple occupancy models for rare species

require(R2jags)

model0<-"Model-const-phi-gam-p.txt"
model1<-"Model-Time-spec-phi-gam-p(.).txt"

# Start with with the trinary array for each species (spMat)
# This array is nsites x 15 x nyears

# Fitting the occupancy models to each species
# MCMC settings
ni <- 50000 # Number of iterations per chain in the Gibbs sampler
nt <- 3    # Thining rate every third observation is discarded
nb <- 49000 # Burn-in rate. First nb iterations are discarded
nc <- 3    # Number of chains


#Values for JAGS model

params <- c("psi","gamma","phi","p")

initial <- apply(spMat, c(1,3), max, na.rm = TRUE)
initial[initial=="-Inf"] < -NA # remove the -Inf's that result when the camera trap 
inits <- function(){ list(z = initial)}
jags.data <- list(y = spMat, nsite = dim(spMat)[1], nrep = dim(spMat)[2], nyear = dim(spMat)[3])


#run the model - jags() function calls in JAGS and runs the MCMC chain

model.fit <- jags(jags.data, inits, params, "Model-Time-spec-phi-gam-p(.).txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb))

#extracting psi values should proceed the same as in model with covariates