# Code to estimate occupancy for common species (>10 detections per year) using the following covariates:
# camera trap point covariates: edge, elevation, people,
# yearly covariates: rain6 maxtemp6, mintemp6

require(R2jags)
require(plotrix)
require(ggplot2)

# JAGS model code
sink("fullmodel7.txt")
cat("
     model {
    
    # Priors for constant parameters
    
    alpha ~ dnorm(0,0.01)
    alphap ~ dnorm(0,0.01)
    alphag ~ dnorm(0,0.01)
    
    # Priors for parameters associated with a covariate, but the parameter is time invariant
    for(i in 1:11){
    beta[i] ~ dnorm(0,0.01)
    w[i] ~ dbern(0.5)
    }
    
    
    #Ecological Model-define gamma and phi
    
    for (i in 1:nsite){
      for (t in 2:nyear){
        logit(gamma[i,t-1])<- alphag + w[4]*beta[4]*rainfall6[t]
        + w[5]*beta[5]*maxtemp6[t]
        + w[6]*beta[6]*mintemp6[t]
        + w[10]*beta[10]*people[i,t]
    
    
    
        logit(phi[i,t-1])<-alphap + w[7]*beta[7]*rainfall6[t]
        + w[8]*beta[8]*maxtemp6[t]
        + w[9]*beta[9]*mintemp6[t]
        + w[11]*beta[11]*people[i,t]
      }
    }
    
    #Ecological Model - define psi1 and subsequent ones
    for (i in 1:nsite){
      logit(psi1[i])<-alpha + w[1]*beta[1]*elevation[i,1]
      + w[2]*beta[2]*people[i,1]
      + w[3]*beta[3]*edge[i,1]
    
      z[i,1] ~ dbern(psi1[i])
    
        for (t in 2:nyear){
          muZ[i,t]<- z[i,t-1]*phi[i,t-1] + (1-z[i,t-1])*gamma[i,t-1]
          z[i,t] ~ dbern(muZ[i,t])
        }#t
      } #i
    
    #Observation process - no covariates
    for(t in 1:nyear){
      p[t] ~ dunif(0,1)    
    } #t
    
    for (i in 1:nsite){
      for (j in 1:nrep){
        for (t in 1:nyear){
          muy[i,j,t] <- z[i,t]*p[t]
          y[i,j,t] ~ dbern(muy[i,j,t])
        } #t
      } #j
    } #i
    
    # Derived quantities (finite sample estimates)
    
    for(t in 1:nyear){
      psi[t]<-sum(z[1:nsite,t])/nsite
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
load("Cunpac")

#Fit Cunpac species

initial<-apply(Cunpac,c(1,3),max,na.rm=T)
initial[initial=='-Inf']<-NA
nsites<-apply(initial,2,function(x) sum(!is.na(x)))
params <- c("psi","beta","w","alpha","alphag","alphap")

inits <- function(){ list(z = initial)}
jags.data <- list(y = Cunpac, nsite = dim(Cunpac)[1], nrep = dim(Cunpac)[2], nyear = dim(Cunpac)[3],elevation=elevation, edge=edge,people=people,rainfall6=rainfall6,maxtemp6=maxtemp6,mintemp6=mintemp6)

# MCMC settings

ni <- 30000 # Number of iterations per chain in the Gibbs sampler
nt <- 3    # Thining rate every third observation is discarded
nb <- 10000 # Burn-in rate. First nb iterations are discarded
nc <- 3    # Number of chains

#Run JAGS
outCunpac <- jags(jags.data, inits, params, "fullmodel7.txt", n.chains = nc,n.thin = nt, n.iter = ni, n.burnin = nb)
print(outCunpac)

#get the psi's
psiMatCunpac<-outCunpac$BUGSoutput$sims.list$psi

# Extracting parameters ---------------------------------------------------

#Read in a table with the names of covariates, associated parameters and indicators.
parTable<-read.csv("Variables-indicators table-test code-simplified model.csv",h=T)
#determine which parameters are important for Daspun species
#Extract matrix of results from JAGS object
parmatCunpac<-outCunpac$BUGSoutput$sims.matrix
# figure out which columns have the indicator variables
indx<-grep("w",colnames(parmatCunpac))
#Put the indicators (w's) in a different matrix
ImatCunpac<-parmatCunpac[,indx]
#Calculate weigths
WCunpac<-apply(ImatCunpac,2,mean)

SigWCunpac<-sort(WCunpac[WCunpac>0.5], dec=T) #The parameters associated with these indicators are the variables that are contributing to the changes in occupancy of this species 

#Extract these w's from the larger ImatCunpac matrix
cols_w <- which(colnames(ImatCunpac) %in% names(SigWCunpac))
SigWCunpac_mat<-ImatCunpac[,cols_w]

#Subset the parTable to get the parameters that matter for Cunpac species
pars_Cunpac<-subset(parTable,Indicator.parameter %in% names(SigWCunpac))

#Get the values of these parameters from the parameter matrix
#First get the columns in the matrix corresponding to those parameters
cols <- which(colnames(parmatCunpac) %in% as.character(pars_Cunpac$Model.Parameter))

#Now, get the values
sub_parmatCunpac<-parmatCunpac[,cols]

#sub_parmatCunpac has all the values for the betas, but we need to ignore the values where the w for that particular parameter was 0

#First remove the 0's and replace by NAs
SigWCunpac_mat<-ifelse(SigWCunpac_mat<1,NA,SigWCunpac_mat)
#Now do an element by element multiplication of the indicator and parameter
#value matrices
sub_parmatCunpac<-sub_parmatCunpac*SigWCunpac_mat

#Calculate the mean(),median, st dev. and 95% confidence intervals for the beta's
mean.betas<-apply(as.matrix(sub_parmatCunpac),2,mean,na.rm=T)
median.betas<-apply(as.matrix(sub_parmatCunpac),2,median,na.rm=T)
sd.betas<-apply(as.matrix(sub_parmatCunpac),2,sd,na.rm=T)
conf.betas<-apply(as.matrix(sub_parmatCunpac),2,quantile,na.rm=T,c(0.025,0.975))

#add to the beta table
pars_Cunpac<-data.frame(pars_Cunpac,mean.betas=mean.betas, median.betas=median.betas,sd.betas=sd.betas,lo95=conf.betas[1,],hi95=conf.betas[2,])