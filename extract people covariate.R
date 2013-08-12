# Code to extract human presence at a TEAM site

#download data first

require(lubridate)
source("camera trap analysis code-04-26-13.R")
#load data from CT csv downloaded from teamnetwork.org
data<-f.readin.fix.data()
nyears<-length(year<-unique(data$Sampling.Period))

mat<-list()
for(i in 1:nyears)
  mat[[i]]<-f.matrix.creator2(data,year[i])  


shmat<-list()
for(i in 1:nyears)
  shmat[[i]]<-lapply(mat[[i]],FUN=f.shrink.matrix.to15)  


spnames<-names(mat[[1]])
indx<-which(spnames=="Homo sapiens")
#extract species 8 ready for analysis Homo sapiens
Homsap<-f.multyear.sp(shmat,indx)
people<-apply(Homsap,c(1,3),max,na.rm=T)
(people[people=='-Inf']<-NA)
save(people,file="covariates/people")