#Demo code to calculate WPI
#Jorge Ahumada - Conservation International - Copyright 2013

require(R2jags)
require(ggplot2)

#Aside from these requirements, need to have JAGS installed in the machine


#source file with some additional helper functions
source("helper functions WPI.R")

#read in covariates
covs<-read.csv("data/covariates.csv",h=T)
ELEVATION<-(covs[,3]-mean(covs[,3]))/sd(covs[,3])
CANOPY <-(covs[,5]-mean(covs[,5]))/sd(covs[,5])


#import species list

splist<-read.csv("data/spVB.txt",h=F)
splist<-as.character(splist[,1])

path<-"data/"
#create some code to import the data and put into a matrix
mat<-array(NA,dim=c(60,15,5,length(splist)))

for(s in 1:length(splist)){
  sp<-read.csv(paste(path,splist[s],sep=""))
  sp<-sp[,-1]

  nyears<-5
  
  #transfer the data into a three dimensional array
  sec<-seq(1,75,by=15)
  for(i in 1:nyears)
    mat[,,i,s]<-as.matrix(sp[,sec[i]:(sec[i]+14)]) # need to put as.matrix() here.. otherwise does not work
  
}
#First step is to calculate occupancies of each species

# Species 1 ---------------------------------------------------------------

#Species 1 - Cuniculus paca
# main data and parameters of JAGS model are passed on as a list
jags.data <- list(y = mat[,,,1], nsite = dim(mat)[1], nrep = dim(mat)[2], nyear = dim(mat)[3])

#Initial values for JAGS model
initial <- apply(mat[,,,1], c(1,3), max, na.rm = TRUE)
#clean up and put NAs
initial[initial=="-Inf"]<-NA # remove the -Inf's that result when the camera trap was stolen
# initial values need to be passed as a function to JAGS
inits <- function(){ list(z = initial)}



# MCMC settings
ni <- 8000 # Number of iterations per chain in the Gibbs sampler
nt <- 3    # Thining rate every third observation is discarded
nb <- 3000 # Burn-in rate. First nb iterations are discarded
nc <- 5    # Number of chains

# parameters to be monitored by JAGS. These are model parameters and derived values
params <- c("psi", "phi", "gamma", "p", "fit","fit.new") 
#run the model - jags() function calls in JAGS and runs the MCMC chain
#model is in a text file in the default path - working directory
outm0s1 <- jags(jags.data, inits, params, "models/Dynocc-jags-m0-WPI-PPC.txt", n.chains = nc, 
            n.thin = nt, n.iter = ni, n.burnin = nb)

#check results
print(outm0s1)
#calculate goodness of fit
f.ppc(model=outm0s1)
#graph occupancy dynamics
graph.psi(outm0s1$BUGSoutput$sims.list$psi,fun=f.mode,"Cuniculus paca")


# Species 2 -----------------------------------------
#Species 2 - Dasyprocta punctata - model 2



initial <- apply(mat[,,,2], c(1,3), max, na.rm = TRUE)
initial[initial=="-Inf"]<-NA # remove the -Inf's that result when the camera trap was stolen
#number of sites with camera traps in year 1 needs tobe passed in jags.data for this model
nsite1<-sum(!is.na(initial[,1]))
inits <- function(){ list(z = initial)}

# Notice now there are some covariates included
jags.data <- list(y = mat[,,,2], nsite = dim(mat)[1], nrep = dim(mat)[2], nyear = dim(mat)[3],nsite1=nsite1, elevation=ELEVATION,canopy=CANOPY)

params <- c("psi", "phi", "gamma", "p", "fit","fit.new","beta0","beta1","beta2")
#run the model
outm2.1s2 <- jags(jags.data, inits, params, "models/Dynocc-jags-m2.1-WPI-PPC.txt", n.chains = nc, 
                n.thin = nt, n.iter = ni, n.burnin = nb)
#check for goodness of fit

f.ppc(outm2.1s2)
#look at results
print(outm2.1s2)
#graph occupancy predicted by the model
graph.psi(outm2.1s2$BUGSoutput$sims.list$psi,initial,fun=f.mode,"Dasyprocta punctata")


# Species 3 ---------------------------------------------------------------
#Dasypus novemcinctus

initial <- apply(mat[,,,3], c(1,3), max, na.rm = TRUE)
initial[initial=="-Inf"]<-NA # remove the -Inf's that result when the camera trap was stolen
nsite1<-sum(!is.na(initial[,1]))
inits <- function(){ list(z = initial)}
jags.data <- list(y = mat[,,,3], nsite = dim(mat)[1], nrep = dim(mat)[2], nyear = dim(mat)[3],nsite1=nsite1, elevation=ELEVATION,canopy=CANOPY)
params <- c("psi", "phi", "gamma", "p", "fit","fit.new","beta0","beta1","beta2","beta3") 

outm2s3 <- jags(jags.data, inits, params, "models/Dynocc-jags-m2-WPI-PPC.txt", n.chains = nc, 
                n.thin = nt, n.iter = ni, n.burnin = nb)

print(outm2s3)
plot(outm2s3)
f.ppc(outm2s3)
graph.psi(outm2s3$BUGSoutput$sims.list$psi,initial,fun=f.mode,"Dasypus")



# Put together WPI --------------------------------------------------------


#First extract the psi distributions for each species

psi_sp1<-outm0s1$BUGSoutput$sims.list$psi
psi_sp2<-outm2.1s2$BUGSoutput$sims.list$psi
psi_sp3<-outm2s3$BUGSoutput$sims.list$psi

psi_all<-abind(psi_sp1,psi_sp2,psi_sp3,along=3)

#calculate the wpi posteriors
wpi_all<-f.WPI(psi_all)

graph.WPI(wpi_all,median,"WPI for three species")