######################################process data#######################################
source("./load.R")
gen<-read.csv(file="./07.HY-MISS20P.Map-Genotype.csv", header=T);
stem_16<-read.table("./zhugao_16.txt",header=F,na.strings = c("NA"))
taproot_16<-read.table("./zhugenchang_16_new.txt",header=F,na.strings = c("NA"))
taproot_16<-taproot_16[taproot_16[,1]%in%stem_16[,1],]
dat1<- load.data(datt=gen,pheno=matrix(as.numeric(unlist(stem_16)),ncol=17,nrow=321,byrow=FALSE)[,-2],nstart=1,nlen=8305, times=c(3,5,7,18,20,22,24,26,28,31,34,38,47,54,62))
dat2<- load.data(datt=gen,pheno=matrix(as.numeric(unlist(taproot_16)),ncol=17,nrow=321,byrow=FALSE)[,-2],nstart=1,nlen=8305, times=c(3,5,7,18,20,22,24,26,28,31,34,38,47,54,62))
#save(dat1,file="stem.Rdata")
#save(dat2,file="taproot.Rdata")



##################################################################################
library(deSolve)
library(ggplot2)
library(mvtnorm)
library(patchwork)

load("stem.Rdata")#dat1
load("taproot.Rdata")#dat2
load("N.s.Rdata")#dat3
load("L.s.Rdata")#dat4

source("2pairs.ode.R")
source("2pairs.covar.R")
source("2pairs.curve.R")
source("2pairs.optim.R")



init.par<-pairs.H0(dat1,dat2,parin2=c(0.5,11.8071,0.5,23.3910,0.5,182.1515,0.2040,0.1694,-0.0394,20.1145,0.1193,1.1515,2.3619))
ret1 <- pairs.est1(dat1,dat2,
                   parin2=c(0.5,11.8071,0.5,23.3910,0.5,182.1515,0.2040,0.1694,-0.0394,20.1145,0.1193,1.1515,2.3619),interval = c(1,4000))
#save(ret1,file="ret-stem-taproot.Rdata")


init.par<-pairs.H0(dat3,dat4,parin2=c(0.5,14.3665,0.5,19.7714,0.5,57.3800,0.1864,0.2190,0.0141,61.4581,0.2565,0.1781,0.0731))
ret <- pairs.est1(dat3,dat4,interval = c(1,2),
                      parin2=c(0.5,14.3665,0.5,19.7714,0.5,57.3800,0.1864,0.2190,0.0141,61.4581,0.2565,0.1781,0.0731))
#save(ret,file="ret-number-length.Rdata")







########################################network###########################################
source("2pairs.network.R")
sigqtl<-read.csv(file="./result/sigqtl2.csv", header=T)
sigqtl.index<-sigqtl[,1]
effects<-SD(dat3,dat4,sigqtl.index)

t<-dat3$sample_times
all_para1 <- t(apply(effects[,c(15:(2*length(t)))], 1,smooth.optim_ind,times=t,
                     para=rep(0.01,6)))
rownames(all_para1) <- c(1:dim(effects)[1])

smooth_data1 <- t(apply(all_para1, 1, Legendre.model,t=seq(1,max(t),0.5)))
rownames(smooth_data1) <- c(1:dim(effects)[1])

Net<-optim_inter(all_cluster_value=smooth_data1,ind_Lpar11=all_para1,norder=5,time=t)
save(Net,file="Net-4-L.Rdata")
#load("Net-1-stem.Rdata")
result<-out_Netdata_cytoscape(optim_cluster=Net,filename="network4.L.csv",clusterAA=smooth_data1)







