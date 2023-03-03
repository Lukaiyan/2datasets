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



########################################PLOT########################################
library(ggplot2)
library(dplyr)
library(patchwork)
library(deSolve)

source("2pairs.plot.R")

allviolin.plot1(dat1,dat2,color="#56B4E9",ylab1="Stem length (mm)",ylab2="Taproot length (mm)",filename="violinplot1.tif")
allviolin.plot1(dat3,dat4,color="#E69F00",ylab1="Root number (count)",ylab2="Root length (cm)",filename="violinplot2.tif")

cv.plot(dat1,dat2,color="#56B4E9",filename="cv1.tif")
cv.plot(dat3,dat4,color="#E69F00",filename="cv2.tif")



#pairs.mean(dat1,dat2,init.par=c(84.1011,0.1597,0.2956,0.0658,87.8234,0.3284,0.1169,0.2290))
pairs.mean(dat1,dat2,color="#56B4E9",init.par=c(182.1515,0.2040,0.1694,-0.0394,20.1145,0.1193,1.1515,2.3619))
pairs.mean1(dat3,dat4,color="#E69F00",init.par=c(57.3800,0.1864,0.2190,0.0141,61.4581,0.2565,0.1781,0.0731))



pairs.fit(dat1,dat2,color="#56B4E9",init.par=c(182.1515,0.2040,0.1694,-0.0394,20.1145,0.1193,1.1515,2.3619),
          lpar1=c(71.8931,0.0897,193.497),gpar1=c(409.2132,0.0201,7.9497),rpar1=c(113.8355,0.0314,0.5157,0.9008),
          lpar2=c(167.0514,0.0574,33.5669),gpar2=c(797.0535,0.0140,5.3180),rpar2=c(284.7532,0.0208,0.1581,0.9637))
pairs.fit1(dat3,dat4,color="#E69F00",init.par=c(57.3800,0.1864,0.2190,0.0141,61.4581,0.2565,0.1781,0.0731),
           lpar1=c(50.5633,0.0695,38.3047),gpar1=c(67.6613,0.0322,5.3221),rpar1=c(67.1708,0.0313,0.2922,0.9402),
           lpar2=c(68.6949,0.0731,53.3358),gpar2=c(95.4887,0.0325,5.9905),rpar2=c(113.8355,0.0333,0.245,0.9700))




####heatmap
library(pheatmap)
library(gplots)
library(RColorBrewer)
library(viridis)
library(factoextra)

source("2pairs.plot.R")
sigqtl<-read.csv(file="./result/sigqtl1.csv", header=T)
all_sigqtl<-sigqtl[,1]
H2.plot(dat1,all_sigqtl,filename="heatmap1.tif",cutree_num=3)
H2.plot11(dat2,all_sigqtl,filename="heatmap2.tif",cutree_num=3)


sigqtl<-read.csv(file="./result/sigqtl2.csv", header=T)
all_sigqtl<-sigqtl[,1]
H2.plot1(dat3,all_sigqtl,filename="heatmap3.tif",cutree_num=3)
H2.plot1(dat4,all_sigqtl,filename="heatmap4.tif",cutree_num=3)




sigqtl<-read.csv(file="./result/sigqtl1.csv", header=T)
all_sigqtl<-sigqtl[,1]
H2.plot1(dat1,dat2,all_sigqtl,snpname=c(4017,5033),rowname=c("stem","taproot"),filename="heatmap3505.1.tif")

sigqtl<-read.csv(file="./result/sigqtl2.csv", header=T)
all_sigqtl<-sigqtl[,1]
H2.plot1(dat3,dat4,all_sigqtl,snpname=c(2601,8300),rowname=c("root number","root length"),filename="heatmap3505.2.tif")


bar.plot(twotrait=c("Stem","Taproot"),
         net1=read.csv("./result/network1.stem.csv"),
         net2=read.csv("./result/network2.taproot.csv"))

bar.plot(twotrait=c("Root number","Root length"),
         net1=read.csv("./result/network3.N.csv"),
         net2=read.csv("./result/network4.L.csv"))







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







