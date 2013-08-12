#HELPER functions in WPI analysis and rare species analyses
#create some fake data for a species
#that is rare

#function to calculate the mode of a distribution
f.mode<-function(data){
	qwe<-density(data)
	qwe$x[which(qwe$y==max(qwe$y))]
	
}

#function to generate the WPI from the output simulations in JAGS
# psi is a three dimensional matrix with the psi of each species in each year
f.WPI <-function(psi){
  nsim <- dim(psi)[1]
  nyears <- dim(psi)[2]
  nsp <- dim(psi)[3]
  rel_psi<-numeric()
  wpi<-matrix(NA,nr=nsim,nc=nyears)
  for(i in 1:nsim){
    for(t in 1:nyears){
      for(s in 1:nsp){
        rel_psi[s] <- log(psi[i,t,s]/psi[i,1,s])
      }
      wpi[i,t]<-exp(1/nsp*sum(rel_psi))
    }
  }
  wpi
  
}

#graph the WPI through time with 95% confidence limits
#WPI is a matrix of n x t values (n = number of runs, t=number of years)
#calculate with mean, mode or median
graph.WPI <- function(wpi,fun=mean,title){
  
  require(ggplot2)
  FUN<-match.fun(fun)
  ct<-apply(wpi,2,FUN)
  lo<-apply(wpi,2,quantile,0.025)
  hi<-apply(wpi,2,quantile,0.975)
  res<-data.frame(year=2008:2012,ct=ct,lo=lo,hi=hi)
  #res<-melt(res,id.vars=c('year'))
  
  p<-ggplot(data=res, aes(x=year))
  p<-p+geom_line(aes(y=ct),size=2)
  p<-p+geom_ribbon(aes(ymin=lo,ymax=hi),alpha=0.2)+xlab("Year")+ylab("Wildlife Picture Index")+labs(title="")+geom_hline(yintercept=1,size=0.5,linetype=2) + ylim(0,2)+labs(title=title)
  p
  
  #ggsave("SpeciesRichness.pdf",p,width=15,height=8,units="cm")
}

graph.psi <- function(psi,initial,fun=mean,title){
  
  require(ggplot2)
  FUN<-match.fun(fun)
  ct<-apply(psi,2,FUN)
  lo<-apply(psi,2,quantile,0.025)
  hi<-apply(psi,2,quantile,0.975)
  naive<-apply(initial,2,function(x) sum(x,na.rm=T)/sum(!is.na(x)))
  res<-data.frame(year=2008:2012,ct=ct,lo=lo,hi=hi)
  #res<-melt(res,id.vars=c('year'))
  
  p<-ggplot(data=res, aes(x=year))
  p<-p+geom_line(aes(y=ct),size=2)+geom_point(y=naive,size=3)
  p<-p+geom_ribbon(aes(ymin=lo,ymax=hi),alpha=0.2)+xlab("Year")+ylab("Occupancy")+labs(title="")+ ylim(0,1)+labs(title=title)
  p
  
  #ggsave(paste("Occ_",title,".pdf",sep=""),p,width=15,height=8,units="cm")
}

#function to check for posterior predictive checks

f.ppc<-function(model){
  fit<-model$BUGSoutput$sims.list$fit
  fit.new<-model$BUGSoutput$sims.list$fit.new
  plot(fit,fit.new)
  abline(0,1)
  return(mean(fit.new>fit))
  
}
