# Code to estimate occupancy for any species using the following covariates:
# camera trap point covariates: edge, elevation, people,
# yearly covariates: rain6(6 months before), rain12(12 months before),
# maxtemp6, maxtemp12, mintemp6, mintemp12,
# year, edge vs. year, elevation vs. year,  people vs. year
require(R2jags)
require(plotrix)
require(ggplot2)
source("helper functions WPI.R")

#Read in a table with the names of covariates, associated parameters and indicators.
#This file will change from site to site based on the number of years of data available for that site
parTable<-read.csv("Variables-indicators table-test code.csv",h=T)
# JAGS model code
sink("fullmodel3.txt")
cat("
    model {
    
    # Priors for constant parameters
    for (i in 1:4){
    alpha[i] ~ dunif(-10,10)
    }
    
    # Priors for parameters associated with a covariate, but the parameter is time invariant
    # and priors for associated w's
    for(i in 1:12){
    beta[i] ~ dunif(-10,10)
    w[i] ~ dbern(0.5)
    }
    # Priors for parameters that are time variant and associated w's
    for (i in 1:4){
    for(j in 1:4){
    # for gammas
    betag[i,j] ~ dunif(-10,10)
    
    #for phi's
    betaph[i,j] ~ dunif(-10,10)
    #for p's
    # w's for gammas
    wg[i,j] ~ dbern(0.5)
    wph[i,j] ~ dbern(0.5)
    
    }
    }
    # Last betap and wp
    for(i in 1:nyear){
    for(j in 1:4) {
    betap[i,j] ~ dunif(-10,10)
    wp[i,j] ~ dbern(0.5)
    }#j
    }#i
    #Ecological Model-define gamma and phi
    
    for (i in 1:nsite){
    for (t in 2:nyear){
    logit(gamma[i,t-1])<- alpha[2]
    +w[4]*beta[4]*rainfall6[t]
    +w[5]*beta[5]*maxtemp6[t]
    +w[6]*beta[6]*mintemp6[t]
    +wg[t-1,1]*betag[t-1,1]*elevation[i,t]
    +wg[t-1,2]*betag[t-1,2]*edge[i,t]
    +wg[t-1,3]*betag[t-1,3]*people[i,t]
    +wg[t-1,4]*betag[t-1,4]*year[t]
    
    
    logit(phi[i,t-1])<- alpha[3]
    +w[7]*beta[7]*rainfall6[t]
    +w[8]*beta[8]*maxtemp6[t]
    +w[9]*beta[9]*mintemp6[t]
    +wph[t-1,1]*betaph[t-1,1]*elevation[i,t]
    +wph[t-1,2]*betaph[t-1,2]*edge[i,t]
    +wph[t-1,3]*betaph[t-1,3]*people[i,t]
    +wph[t-1,4]*betaph[t-1,4]*year[t]
    
    }
    }
    
    #Ecological Model - define psi1 and subsequent ones
    for (i in 1:nsite){
    logit(psi1[i])<-alpha[1]+w[1]*beta[1]*elevation[i,1]+w[2]*beta[2]*people[i,1]+w[3]*beta[3]*edge[i,1]
    z[i,1] ~ dbern(psi1[i])
    for (t in 2:nyear){
    muZ[i,t]<- z[i,t-1]*phi[i,t-1] + (1-z[i,t-1])*gamma[i,t-1]
    z[i,t] ~ dbern(muZ[i,t])
    } #t
    } #i
    
    #Observation process
    for(t in 1:nyear){
    for(i in 1:nsite){
    logit(p[i,t]) <- alpha[4]
    +w[10]*beta[10]*rainfall6[t]
    +w[11]*beta[11]*maxtemp6[t]
    +w[12]*beta[12]*mintemp6[t]
    +wp[t,1]*betap[t,1]*elevation[i,t]
    +wp[t,2]*betap[t,2]*edge[i,t]
    +wp[t,3]*betap[t,3]*people[i,t]
    +wp[t,4]*betap[t,4]*year[t] 
    } #i
    } #t
    
    for (i in 1:nsite){
    for (j in 1:nrep){
    for (t in 1:nyear){
    muy[i,j,t] <- z[i,t]*p[i,t]
    y[i,j,t] ~ dbern(muy[i,j,t])
    } #t
    } #j
    } #i
    
    # Derived quantities
    psi[1] <- sum(z[,1])/nsites[1]
    for(t in 2:nyear){
    for(i in 1:nsite){
    tpsi[i,t-1] <- z[i,t-1]*phi[i,t-1] + (1-z[i,t-1])*gamma[i,t-1]
    }#i
    psi[t]<-sum(tpsi[,t-1])/nsites[t]
    }#t
    }
    ",fill = TRUE)
sink()


#COVARIATES
elevationVec<-c(-0.8432954,-0.8834697,-0.7987446,-0.9081438,-0.8443755,-0.9012985,-0.8616432,-0.8212585,-0.7828937,-0.7468715,-0.7700867,-0.9082139,-0.7418497,-0.9106827,-0.9018315,-0.8951685,-0.8807063,-0.8247934,-0.8412895,-0.8307269,-0.6632405,-0.4463779,-0.4166399,-0.3009004,-0.2136083,-0.134424,-0.04164716,0.09453022,-0.3754837,-0.280673,-0.2325031,-0.6320297,-0.1641198,-0.5763973,-0.514691,-0.4916161,-0.4851354,-0.3649071,-0.4555658,-0.5518775,0.2508088,1.17089,1.25828,1.403,1.540384,1.644803,1.663656,1.851973,1.899792,2.375992,2.483877,0.1804196,2.60663,0.04800155,0.2460395,0.4137082,0.563408,0.6972568,0.8337989,1.011932)

elevation <- array(dim=c(60,5))
elevation[,1] <- elevationVec;
elevation[,2] <- elevationVec;
elevation[,3] <- elevationVec;
elevation[,4] <- elevationVec;
elevation[,5] <- elevationVec;

edgeVec <-c(-0.2976453,-1.725904,0.692566,-1.733318,-0.3659969,-1.59561,-0.317112,0.7639836,0.6112632,1.034628,0.5466121,-1.666433,-0.3876706,-0.9422138,-0.8748138,-1.675684,-0.7825684,-0.7298378,-1.247774,0.3784161,0.05678602,1.057345,1.446243,0.4794103,-0.138555,-1.043869,-0.878316,0.1654189,0.007821922,1.291765,-0.1050136,0.37498,1.215061,-0.3613317,-0.2829892,0.3299806,-0.9198529,-0.7981497,0.5859684,-0.04338885,2.591052,1.241968,0.8390619,0.7367593,1.188696,0.2600036,-0.3253718,-0.6390461,-1.502506,-0.5444218,-1.430388,2.137859,-0.6034562,0.7954237,-0.3710452,-0.3237595,0.2881794,0.8981492,1.187533,1.451107)
edge <- array(dim=c(60,5))
edge[,1] <- edgeVec;
edge[,2] <- edgeVec;
edge[,3] <- edgeVec;
edge[,4] <- edgeVec;
edge[,5] <- edgeVec;



people<-array(dim=c(60,5))
for(i in 1:5) {
  a<-rnorm(1,0,0.05)
  cat(a)
  people[,i]<-rbinom(n=60,1,0.4+a)
}

rainfall6<-c(730.2131,1843.735,845.0721,1094.583,1169.649)
#standardize
rainfall6<-(rainfall6-mean(rainfall6))/sd(rainfall6)

maxtemp6<-c(32.85997,34.54642,29.00966,34.22395,29.62746)
#standardize
maxtemp6<-(maxtemp6-mean(maxtemp6))/sd(maxtemp6)

mintemp6<-c(15.60596,10.00254,16.31606,17.79081,14.2639)
#standardize
mintemp6<-(mintemp6-mean(mintemp6))/sd(mintemp6)

year<-1:5


#Get some species data - these are R objects - attached with original email
# These species are from site 1
load("Daspun")
load("Cunpac")

#Fit Daspun species first
initial<-apply(Daspun,c(1,3),max,na.rm=T)
initial[initial=='-Inf']<-NA
nsites<-apply(initial,2,function(x) sum(!is.na(x)))

inits <- function(){ list(z = initial)}
jags.data <- list(y = Daspun, nsite = dim(Daspun)[1], nrep = dim(Daspun)[2], nyear = dim(Daspun)[3],elevation=elevation, edge=edge,people=people,rainfall6=rainfall6,maxtemp6=maxtemp6,mintemp6=mintemp6,year=year,nsites=nsites)

# MCMC settings

ni <- 30000 # Number of iterations per chain in the Gibbs sampler
nt <- 3    # Thining rate every third observation is discarded
nb <- 10000 # Burn-in rate. First nb iterations are discarded
nc <- 3    # Number of chains



params <- c("psi","beta","betag","betaph","betap","w","wg","wph","wp","alpha")

#run the model - jags() function calls in JAGS and runs the MCMC chain
#model is in a text file in the default path - working directory
# takes about 24 min to run in my machines
system.time(out <- jags(jags.data, inits, params, "fullmodel3.txt", n.chains = nc, 
                        n.thin = nt, n.iter = ni, n.burnin = nb))

print(out)
outDaspun<-out

# Needs more iterations probably good to increase to 50,000

#Fit Cunpac species

initial<-apply(Cunpac,c(1,3),max,na.rm=T)
initial[initial=='-Inf']<-NA
nsites<-apply(initial,2,function(x) sum(!is.na(x)))

inits <- function(){ list(z = initial)}
jags.data <- list(y = Cunpac, nsite = dim(Cunpac)[1], nrep = dim(Cunpac)[2], nyear = dim(Cunpac)[3],elevation=elevation, edge=edge,people=people,rainfall6=rainfall6,maxtemp6=maxtemp6,mintemp6=mintemp6,year=year,nsites=nsites)

# MCMC settings

ni <- 30000 # Number of iterations per chain in the Gibbs sampler
nt <- 3    # Thining rate every third observation is discarded
nb <- 10000 # Burn-in rate. First nb iterations are discarded
nc <- 3    # Number of chains


#This took 120 mins in TEAM server
system.time(outCunpac <- jags(jags.data, inits, params, "fullmodel3.txt", n.chains = nc, 
                              n.thin = nt, n.iter = ni, n.burnin = nb))/60

print(outCunpac)

#Get the psis for these two species

#get the psi's

psiMatDaspun<-outDaspun$BUGSoutput$sims.list$psi
psiMatCunpac<-outCunpac$BUGSoutput$sims.list$psi
psiMatDaspun[psiMatDaspun>1]<-1
psiMatCunpac[psiMatCunpac>1]<-1

#graph the occupancy of one of the species
graph.psi(psi=psiMatCunpac,initial=initial,fun=mean,title="Cuniculus paca")

#Calculate the WPI for site 1
#put them in the same object - note: since each species was run with different
# ni's the largest one needs to be subsetted.
psi_all1<-abind(psiMatDaspun,psiMatCunpac,along=3)

#calculate the wpi posteriors
wpiSite1<-f.WPI(psi_all1)
#graph it
graph.WPI(wpiSite1,median,"WPI for site 1")

# Now take two other species from Site 2 and fit occupancy models to those
# For demo purposes I will assume the covariate values for Site 2 are the same
# but in reality the values would be different between sites

load("Tapbai") # species # 1 for site 2
load("Pectaj") # species # 2 for site 2

#fit an occupancy model to species 1 Tapbai

initial<-apply(Tapbai,c(1,3),max,na.rm=T)
initial[initial=='-Inf']<-NA
nsites<-apply(initial,2,function(x) sum(!is.na(x)))

inits <- function(){ list(z = initial)}
jags.data <- list(y = Tapbai, nsite = dim(Tapbai)[1], nrep = dim(Tapbai)[2], nyear = dim(Tapbai)[3],elevation=elevation, edge=edge,people=people,rainfall6=rainfall6,maxtemp6=maxtemp6,mintemp6=mintemp6,year=year,nsites=nsites)

# MCMC settings

ni <- 30000 # Number of iterations per chain in the Gibbs sampler
nt <- 3    # Thining rate every third observation is discarded
nb <- 10000 # Burn-in rate. First nb iterations are discarded
nc <- 3    # Number of chains



params <- c("psi","beta","betag","betaph","betap","w","wg","wph","wp","alpha")

#run the model - jags() function calls in JAGS and runs the MCMC chain
#model is in a text file in the default path - working directory
# takes about 22 min to run in my machine
system.time(outTapbai <- jags(jags.data, inits, params, "fullmodel3.txt", n.chains = nc, 
                              n.thin = nt, n.iter = ni, n.burnin = nb))

print(outTapbai)

#fit an occupancy model to species 2 Pectaj
save(outDaspun,file="outDaspun_site1")


initial<-apply(Pectaj,c(1,3),max,na.rm=T)
initial[initial=='-Inf']<-NA
nsites<-apply(initial,2,function(x) sum(!is.na(x)))

inits <- function(){ list(z = initial)}
jags.data <- list(y = Pectaj, nsite = dim(Pectaj)[1], nrep = dim(Pectaj)[2], nyear = dim(Pectaj)[3],elevation=elevation, edge=edge,people=people,rainfall6=rainfall6,maxtemp6=maxtemp6,mintemp6=mintemp6,year=year,nsites=nsites)

# MCMC settings

ni <- 30000 # Number of iterations per chain in the Gibbs sampler
nt <- 3    # Thining rate every third observation is discarded
nb <- 10000 # Burn-in rate. First nb iterations are discarded
nc <- 3    # Number of chains



params <- c("psi","beta","betag","betaph","betap","w","wg","wph","wp","alpha")

#run the model - jags() function calls in JAGS and runs the MCMC chain
#model is in a text file in the default path - working directory
# takes about 24 min to run in my machines
system.time(outPectaj <- jags(jags.data, inits, params, "fullmodel3.txt", n.chains = nc, 
                              n.thin = nt, n.iter = ni, n.burnin = nb))

print(outPectaj)

psiMatTapbai<-outTapbai$BUGSoutput$sims.list$psi
psiMatPectaj<-outPectaj$BUGSoutput$sims.list$psi
psiMatTapbai[psiMatTapbai>1]<-1
psiMatPectaj[psiMatPectaj>1]<-1
#graph the occupancy of one of the species
graph.psi(psi=psiMatTapbai,initial=initial,fun=f.mode,title="Tapbai")

#Calculate the WPI for site 2
#put them in the same object - note: since each species was run with different
# ni's the largest one needs to be subsetted.
psi_all2<-abind(psiMatTapbai,psiMatPectaj,along=3)

#calculate the wpi posteriors
wpiSite2<-f.WPI(psi_all2)
#graph it
graph.WPI(wpiSite2,median,"WPI for site 2")

#Consolidate the two WPIs for Site 1 and Site 2 into 1

wpi_all<-rbind(wpiSite1,wpiSite2)
graph.WPI(wpi_all,median,"WPI for both sites")

#determine which parameters are important for Daspun species
#parmatrix<-out$BUGSoutput$sims.matrix

# Extracting parameters ---------------------------------------------------
#Read in a table with the names of covariates, associated parameters and indicators.
#This file will change from site to site based on the number of years of data available for that site
parTable<-read.csv("Variables-indicators table-test code.csv",h=T)
#determine which parameters are important for Daspun species

# First: extract the matrix of indicators and parameters
parmatDaspun<-outDaspun$BUGSoutput$sims.matrix
# figure out which columns have the indicator variables
indx<-grep("w",colnames(parmatDaspun))
#Put the indicators (w's) in a different matrix
ImatDaspun<-parmatDaspun[,indx]
#Calculate weigths
WDaspun<-apply(ImatDaspun,2,sum)/dim(ImatDaspun)[1]

SigWDaspun<-sort(WDaspun[WDaspun>0.5], dec=T) #The parameters associated with these indicators are the variables that are contributing to the changes in occupancy of this species 

#Extract these w's from the larger ImatDaspun matrix
cols_w<-which(colnames(ImatDaspun) %in% names(SigWDaspun))
SigWDaspun_mat<-ImatDaspun[,cols_w]

#Subset the parTable to get the parameters that matter for Daspun species
pars_Daspun<-subset(parTable,Indicator.parameter %in% names(SigWDaspun))

#Get the values of these parameters from the parameter matrix
#First get the columns in the matrix corresponding to those parameters
cols <- which(colnames(parmatDaspun) %in% as.character(pars_Daspun$Model.Parameter))

#Now, get the values
sub_parmatDaspun<-parmatDaspun[,cols]

#sub_parmatDaspun has all the values for the betas, but we need to ignore the values where the w for that particular parameter was 0

#First remove the 0's and replace by NAs
SigWDaspun_mat<-ifelse(SigWDaspun_mat<1,NA,SigWDaspun_mat)
#Now do an element by element multiplication of the indicator and parameter
#value matrices
sub_parmatDaspun<-sub_parmatDaspun*SigWDaspun_mat

#Calculate the mean(),median, st dev. and 95% confidence intervals for the beta's
mean.betas<-apply(sub_parmatDaspun,2,mean,na.rm=T)
median.betas<-apply(sub_parmatDaspun,2,median,na.rm=T)
sd.betas<-apply(sub_parmatDaspun,2,sd,na.rm=T)
conf.betas<-apply(sub_parmatDaspun,2,quantile,na.rm=T,c(0.025,0.975))

#add to the beta table
pars_Daspun<-data.frame(pars_Daspun,mean.betas=mean.betas, median.betas=median.betas,sd.betas=sd.betas,lo95=conf.betas[1,],hi95=conf.betas[2,])

#some graphs (I did not do any pie charts here..so please explore)
#one complication is that some beta's are positive and others are negative
#how to convey this information in a pie chart?

p<-ggplot(data=pars_Daspun, aes(x=covariate))
p+geom_bar(aes(y=mean.betas,fill=variable),stat="identity")+facet_grid(variable~variable.family)
