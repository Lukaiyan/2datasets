allviolin.plot1<-function(dat1,dat2,color,ylab1="Stem length (mm)",ylab2="Taproot length (mm)",filename="violinplot1-2.tif"){
  
  data1<-c()
  for(i in 1:length(dat1$sample_times)){
    tmp<-cbind(dat1$pheno[,i],dat2$pheno[,i],rep(dat1$sample_times[i],dim(dat1$pheno)[1]),n=c(1:dim(dat1$pheno)[1]))
    data1<-rbind(data1,tmp)
  }
  data<-cbind(data1,rep(1,dim(data1)[1]),rep(2,dim(data1)[1]))
  colnames(data)<-c("phe1","phe2","time","n","factor1","factor2")
  data<-as.data.frame(data)
  
  max_height<-c(50,75,100,150,200,250)
  yunit<-c(10,25,25,50,50,50)
  ry11<-max(data[,1])
  ry1 <- max_height[min(which(max_height>=ry11))]
  runit1<-yunit[min(which(max_height>=ry11))]
  
  ry22<-max(data[,2])
  ry2 <- max_height[min(which(max_height>=ry22))]
  runit2<-yunit[min(which(max_height>=ry22))]
  
  
  g1 <- ggplot(data,aes(x=factor(time),y=phe1,fill=factor(factor1)))
  g1 <- g1 + geom_violin(colour=color)##56B4E9
  g1 <- g1 + geom_boxplot(aes(size = factor(factor2)),width=0.05,outlier.colour = NA)
  g1 <- g1 + scale_fill_manual(values = c(color))
  g1 <- g1 + scale_size_manual(values=c(0.2))
  g1 <- g1 + xlab("Time (day)")+ylab(ylab1)
  g1 <- g1 + scale_y_continuous(limits = c(0,ry1),breaks=seq(0,ry1,runit1))
  g1 <- g1 + theme_zg()
  g1 <- g1 + theme(axis.text.x = element_blank(),axis.title.x=element_blank(),
                   axis.title.y = element_text(vjust=-0.2,size=10))
  
  
  g2 <- ggplot(data,aes(x=factor(time),y=phe2,fill=factor(factor2)))
  g2 <- g2 + geom_violin(colour=color)#"#E69F00"
  g2 <- g2 + geom_boxplot(aes(size = factor(factor2)),width=0.05,outlier.colour = NA)
  g2 <- g2 + scale_fill_manual(values = c(color))
  g2 <- g2 + scale_size_manual(values=c(0.2))
  g2 <- g2 + xlab("Time (day)")+ylab(ylab2)
  g2 <- g2 + scale_y_continuous(limits = c(0,ry2),breaks=seq(0,ry2,runit2))
  g2 <- g2 + theme_zg()
  g2 <- g2 + theme(axis.title.y = element_text(vjust=-0.2,size=10))
  #g2 <- g2 + theme(axis.text.x = element_blank(),axis.title.x=element_blank())
  
  g <- g1 + g2 + plot_layout(ncol = 1)
  
  #tiff(filename,width=1800,height=800,res=300)
  tiff(filename,width=1250,height=900,res=300)
  print(g)
  dev.off()
}



cv.plot<-function(dat1,dat2,color,filename="cv.tif"){
 
  y1<-apply(dat1$pheno,2,sd)/colMeans(dat1$pheno)
  y2<-apply(dat2$pheno,2,sd)/colMeans(dat2$pheno)
  cv<-data.frame(time=dat1$sample_times,
                    y1=y1,y2=y2)
  
  p <- ggplot(data,aes(x=factor(time),y=phe))
  p <- p + xlab("Time (day)")+ylab("CV")
  p <- p + geom_line(data = cv,aes(x = time,y = y1),size=0.8,linetype=1,colour=color)
  p <- p + geom_line(data = cv,aes(x = time,y = y2),size=0.8,linetype=2,colour=color)
  p <- p + geom_point(data = cv,aes(x = time,y = y1),size=1,colour=color,fill="white",shape=21)
  p <- p + geom_point(data = cv,aes(x = time,y = y2),size=1,colour=color,fill="white",shape=22)
  p <- p + scale_x_continuous(limits = c(min(dat1$sample_times),max(dat1$sample_times)))
  p <- p + scale_y_continuous(limits = c(0,2),breaks=seq(0,2,0.5))
  p <- p + theme_zg()
  
  tiff(filename,width=700,height=400,res=300)
  print(p)
  dev.off()
}







pairs.mean <- function(dat1,dat2,color,init.par){
  
  X1 <- dat1$pheno
  Y1 <- dat2$pheno
  times<-dat1$sample_times
  
  phe_mean<-colMeans(X1,na.rm=TRUE)
  phe_mean1<-colMeans(Y1,na.rm=TRUE)
  a1 <- ind.get_mu(init.par,times,x1=phe_mean[1],x2=phe_mean1[1])
  a2<- ind.get_mu1(init.par[c(1,2,3,5,6,7)],times,x1=phe_mean[1],x2=phe_mean1[1])
  
  
  n1 <- dim(X1)[1]
  n<-length(dat1$sample_times)
  tn <- c()
  for(i in 1:n1){
    
    nn <- as.numeric(X1[i,])
    nn.1 <- cbind(rep(i,length(nn)),1:n,nn)
    tn <- rbind(tn,nn.1)
  }
  
  tl <- c()
  for(i in 1:n1){
    
    nl <- as.numeric(Y1[i,])
    nl.1 <- cbind(rep(i,length(nl)),1:n,nl)
    tl <- rbind(tl,nl.1)
  }
  colnames(tn) <- c("index","time","pheno")
  tn <- as.data.frame(tn)
  colnames(tl) <- c("index","time","pheno")
  tl <- as.data.frame(tl)
  
  fit.tn <- data.frame(index=rep(1,n),time=times,pheno=a1[1:n],tpheno=phe_mean)
  fit.tn1 <- data.frame(index=rep(1,n),time=times,pheno=a2[1:n])
  fit.tn2 <- data.frame(index=rep(1,n),time=times,pheno=a1[1:n]-a2[1:n])
  
  
  Maxti <- dat1$sample_times
  Minti <- rev(Maxti)
  p2.data1 <- c(a1[1:n]*2,rev(a1[1:n]/1.4))
  p2.data2 <- c(Maxti,Minti)
  p2.data <- data.frame(s=p2.data1,TT=p2.data2)
  
  g1 <- ggplot(tn)
  g1 <- g1 + geom_polygon(data=p2.data,aes(x=TT,y=s),fill=color,alpha=0.3)
  g1 <- g1 + geom_line(data=fit.tn,aes(x=time,y=pheno,group=index),colour="red",size=1)
  #g1 <- g1 + geom_point(data=fit.tn,aes(x=time,y=tpheno,group=index),colour="orange",size=3,pch=2)
  g1 <- g1 + geom_line(data=fit.tn1,aes(x=time,y=pheno,group=index),colour="blue",size=1,linetype=2)
  g1 <- g1 + geom_line(data=fit.tn2,aes(x=time,y=pheno,group=index),colour="blue",size=1,linetype=3)
  g1 <- g1 + scale_x_continuous(limits=c(3,63),breaks=seq(3,63,10),labels=seq(3,63,10))
  g1 <- g1 + geom_hline(yintercept=0,linetype=1,size=0.3)
  g1 <- g1 + scale_y_continuous(limits=c(-5,90),breaks=c(-5,seq(0,90,30)),labels=c(-5,seq(0,90,30)))
  g1 <- g1 + xlab("Time (day)")+ylab("Stem length (mm)") + theme_zg()
  #g1 <- g1 + ggtitle("A")
  #g1 <- g1 + theme(plot.title = element_text(vjust = 0.5,hjust = -0.08,face="bold",family="sans"))
  g1 <- g1 + theme(axis.text.x = element_blank(),axis.title.x=element_blank())
  
  
  fit.tl <- data.frame(index=rep(1,n),time=times,pheno=a1[(n+1):(2*n)],tpheno=phe_mean1,tpheno=phe_mean,lpheno=l2,gpheno=g2,rpheno=r2)
  fit.tl1 <- data.frame(index=rep(1,n),time=times,pheno=a2[(n+1):(2*n)])
  fit.tl2 <- data.frame(index=rep(1,n),time=times,pheno=a1[(n+1):(2*n)]-a2[(n+1):(2*n)])
  
  Maxti <- dat2$sample_times
  Minti <- rev(Maxti)
  p2.data1 <- c(a1[(n+1):(2*n)]*1.6,rev(a1[(n+1):(2*n)]/1.4))
  p2.data2 <- c(Maxti,Minti)
  p2.data <- data.frame(s=p2.data1,TT=p2.data2)
  
  
  g2 <- ggplot(tl)
  g2 <- g2 + geom_polygon(data=p2.data,aes(x=TT,y=s),fill=color,alpha=0.3)
  g2 <- g2 + geom_line(data=fit.tl,aes(x=time,y=pheno,group=index),colour="red",size=1)
  #g2 <- g2 + geom_point(data=fit.tl,aes(x=time,y=tpheno,group=index),colour="orange",size=3,pch=2)
  g2 <- g2 + geom_line(data=fit.tl1,aes(x=time,y=pheno,group=index),colour="blue",size=1,linetype=2)
  g2 <- g2 + geom_line(data=fit.tl2,aes(x=time,y=pheno,group=index),colour="blue",size=1,linetype=3)
  g2 <- g2 + scale_x_continuous(limits=c(0,63),breaks=seq(3,63,10),labels=seq(3,63,10))
  g2 <- g2 + geom_hline(yintercept=0,linetype=1,size=0.3)
  g2 <- g2 + scale_y_continuous(limits=c(0,150),breaks=seq(0,150,50),labels=seq(0,150,50))
  g2 <- g2 + xlab("Time (day)")+ylab("Taproot length (mm)") + theme_zg()
  #g2 <- g2 + ggtitle("B")
  g2 <- g2 + theme(plot.title = element_text(vjust = 0.5,hjust = -0.08,face="bold",family="sans"))
  
  g <- g1 + g2 + plot_layout(ncol = 1)
  
  tiff("stem-taproot.tif",height=1200,width=900,res=300)
  print(g)
  dev.off()
  
  
}

pairs.mean1 <- function(dat3,dat4,color,init.par){
  
  require(ggplot2)
  
  X1 <- dat3$pheno
  Y1 <- dat4$pheno
  times<-dat3$sample_times
  
  phe_mean<-colMeans(X1,na.rm=TRUE)
  phe_mean1<-colMeans(Y1,na.rm=TRUE)
  a1 <- ind.get_mu(init.par,times,x1=phe_mean[1],x2=phe_mean1[1])
  a2<- ind.get_mu1(init.par[c(1,2,3,5,6,7)],times,x1=phe_mean[1],x2=phe_mean1[1])
  
  n1 <- dim(X1)[1]
  n<-length(times)
  tn <- c()
  for(i in 1:n1){
    
    nn <- as.numeric(X1[i,])
    nn.1 <- cbind(rep(i,length(nn)),1:n,nn)
    tn <- rbind(tn,nn.1)
  }
  
  tl <- c()
  for(i in 1:n1){
    
    nl <- as.numeric(Y1[i,])
    nl.1 <- cbind(rep(i,length(nl)),1:n,nl)
    tl <- rbind(tl,nl.1)
  }
  colnames(tn) <- c("index","time","pheno")
  tn <- as.data.frame(tn)
  colnames(tl) <- c("index","time","pheno")
  tl <- as.data.frame(tl)
  
  fit.tn <- data.frame(index=rep(1,n),time=times,pheno=a1[1:n],tpheno=phe_mean)
  fit.tn1 <- data.frame(index=rep(1,n),time=times,pheno=a2[1:n])
  fit.tn2 <- data.frame(index=rep(1,n),time=times,pheno=a1[1:n]-a2[1:n])
  
  Maxti <- dat3$sample_times
  Minti <- rev(Maxti)
  p2.data1 <- c(a1[1:n]*2,rev(a1[1:n]/1.4))
  p2.data2 <- c(Maxti,Minti)
  p2.data <- data.frame(s=p2.data1,TT=p2.data2)
  
  g1 <- ggplot(tn)
  g1 <- g1 + geom_polygon(data=p2.data,aes(x=TT,y=s),fill=color,alpha=0.3)
  g1 <- g1 + geom_line(data=fit.tn,aes(x=time,y=pheno,group=index),colour="red",size=1)
  #g1 <- g1 + geom_point(data=fit.tn,aes(x=time,y=tpheno,group=index),colour="orange",size=3,pch=2)
  g1 <- g1 + geom_line(data=fit.tn1,aes(x=time,y=pheno,group=index),colour="blue",size=1,linetype=2)
  g1 <- g1 + geom_line(data=fit.tn2,aes(x=time,y=pheno,group=index),colour="blue",size=1,linetype=3)
  g1 <- g1 + scale_x_continuous(limits=c(13,78),breaks=seq(13,78,13),labels=seq(13,78,13))
  g1 <- g1 + geom_hline(yintercept=0,linetype=1,size=0.3)
  g1 <- g1 + scale_y_continuous(limits=c(-5,90),breaks=c(-5,seq(0,90,30)),labels=c(-5,seq(0,90,30)))
  g1 <- g1 + xlab("Time (day)")+ylab("Root number (count)") + theme_zg()
  #g1 <- g1 + ggtitle("A")
  #g1 <- g1 + theme(plot.title = element_text(vjust = 0.5,hjust = -0.08,face="bold",family="sans"))
  g1 <- g1 + theme(axis.text.x = element_blank(),axis.title.x=element_blank())
  
  
  fit.tl <- data.frame(index=rep(1,n),time=times,pheno=a1[(n+1):(2*n)],tpheno=phe_mean1)
  fit.tl1 <- data.frame(index=rep(1,n),time=times,pheno=a2[(n+1):(2*n)])
  fit.tl2 <- data.frame(index=rep(1,n),time=times,pheno=a1[(n+1):(2*n)]-a2[(n+1):(2*n)])
  
  Maxti <- dat4$sample_times
  Minti <- rev(Maxti)
  p2.data1 <- c(a1[(n+1):(2*n)]*1.6,rev(a1[(n+1):(2*n)]/1.4))
  p2.data2 <- c(Maxti,Minti)
  p2.data <- data.frame(s=p2.data1,TT=p2.data2)
  
  
  g2 <- ggplot(tl)
  g2 <- g2 + geom_polygon(data=p2.data,aes(x=TT,y=s),fill=color,alpha=0.3)
  g2 <- g2 + geom_line(data=fit.tl,aes(x=time,y=pheno,group=index),colour="red",size=1)
  #g2 <- g2 + geom_point(data=fit.tl,aes(x=time,y=tpheno,group=index),colour="orange",size=3,pch=2)
  g2 <- g2 + geom_line(data=fit.tl1,aes(x=time,y=pheno,group=index),colour="blue",size=1,linetype=2)
  g2 <- g2 + geom_line(data=fit.tl2,aes(x=time,y=pheno,group=index),colour="blue",size=1,linetype=3)
  g2 <- g2 +  scale_x_continuous(limits=c(13,78),breaks=seq(13,78,13),labels=seq(13,78,13))
  g2 <- g2 + geom_hline(yintercept=0,linetype=1,size=0.3)
  g2 <- g2 + scale_y_continuous(limits=c(0,100),breaks=seq(0,100,25),labels=seq(0,100,25))
  g2 <- g2 + xlab("Time (day)")+ylab("Root length (cm)") + theme_zg()
  #g2 <- g2 + ggtitle("B")
  #g2 <- g2 + theme(plot.title = element_text(vjust = 0.5,hjust = -0.08,face="bold",family="sans"))
  
  g <- g1 + g2 + plot_layout(ncol = 1)
  
  tiff("rootnumber-rootlength.tif",height=1200,width=900,res=300)
  print(g)
  dev.off()
  
  
}




pairs.fit <- function(dat1,dat2,color,init.par,lpar1,gpar1,rpar1,lpar2,gpar2,rpar2){
  
  X1 <- dat1$pheno
  Y1 <- dat2$pheno
  times<-dat1$sample_times
  times1<-seq(2,59,3)
  times2<-seq(2,59,3)+1
  times3<-seq(2,59,3)+2
  n<-length(times)
  n1<-length(times1)
  
  phe_mean<-colMeans(X1,na.rm=TRUE)
  phe_mean1<-colMeans(Y1,na.rm=TRUE)
  a1 <- ind.get_mu(init.par,times,x1=phe_mean[1],x2=phe_mean1[1])
  a2<- ind.get_mu1(init.par[c(1,2,3,5,6,7)],times,x1=phe_mean[1],x2=phe_mean1[1])
  
  l1<-lpar1[1]/(1+lpar1[3]*exp(-lpar1[2]*times1))
  g1<-gpar1[1]*exp(-gpar1[3]*exp(-gpar1[2]*times2))
  r1<-rpar1[1]*(1-rpar1[3]*exp(-rpar1[2]*times3))^(1/(1-rpar1[4]))
  l2<-lpar2[1]/(1+lpar2[3]*exp(-lpar2[2]*times1))
  g2<-gpar2[1]*exp(-gpar2[3]*exp(-gpar2[2]*times2))
  r2<-rpar2[1]*(1-rpar2[3]*exp(-rpar2[2]*times3))^(1/(1-rpar2[4]))

  fit.tn <- data.frame(index=rep(1,n),time=times,pheno=a1[1:n],tpheno=phe_mean)
  fit.tn1 <- data.frame(index=rep(1,n1),time1=times1,time2=times2,time3=times3,lpheno=l1,gpheno=g1,rpheno=r1)
  

  g1 <- ggplot(fit.tn)
  g1 <- g1 + geom_line(data=fit.tn1,aes(x=time1,y=lpheno,group=index),colour="#838B83",size=0.5)
  g1 <- g1 + geom_line(data=fit.tn1,aes(x=time2,y=gpheno,group=index),colour="#6959CD",size=0.5)
  g1 <- g1 + geom_line(data=fit.tn1,aes(x=time3,y=rpheno,group=index),colour="#EEEE00",size=0.5)
  g1 <- g1 + geom_line(aes(x=time,y=pheno,group=index),colour=color,size=0.5)
  g1 <- g1 + geom_point(data=fit.tn1,aes(x=time1,y=lpheno,group=index),colour="#838B83",fill="#838B83",size=1,shape=8)
  g1 <- g1 + geom_point(data=fit.tn1,aes(x=time2,y=gpheno,group=index),colour="#6959CD",fill="#6959CD",size=1,shape=22)
  g1 <- g1 + geom_point(data=fit.tn1,aes(x=time3,y=rpheno,group=index),colour="#EEEE00",fill="#EEEE00",size=1,shape=23)
  #g1 <- g1 + geom_point(aes(x=time,y=pheno,group=index),colour=color,size=1,shape=24)
  g1 <- g1 + geom_point(aes(x=time,y=tpheno,group=index),colour="black",size=1)
  g1 <- g1 + scale_x_continuous(limits=c(3,63),breaks=seq(3,63,10),labels=seq(3,63,10))
  g1 <- g1 + geom_hline(yintercept=0,linetype=1,size=0.3)
  g1 <- g1 + scale_y_continuous(limits=c(0,60),breaks=seq(0,60,20),labels=seq(0,60,20))
  g1 <- g1 + xlab("Time (day)")+ylab("Stem length (mm)") + theme_zg()
  #g1 <- g1 + ggtitle("A")
  #g1 <- g1 + theme(plot.title = element_text(vjust = 0.5,hjust = -0.08,face="bold",family="sans"))
  g1 <- g1 + theme(axis.text.x = element_blank(),axis.title.x=element_blank())
  
  
  fit.tl <- data.frame(index=rep(1,n),time=times,pheno=a1[(n+1):(2*n)],tpheno=phe_mean1,tpheno=phe_mean)
  fit.tl1 <- data.frame(index=rep(1,n1),time1=times1,time2=times2,time3=times3,lpheno=l2,gpheno=g2,rpheno=r2)
  
  g2 <- ggplot(fit.tl)
  g2 <- g2 + geom_line(data=fit.tl1,aes(x=time1,y=lpheno,group=index),colour="#838B83",size=0.5)
  g2 <- g2 + geom_line(data=fit.tl1,aes(x=time2,y=gpheno,group=index),colour="#6959CD",size=0.5)
  g2 <- g2 + geom_line(data=fit.tl1,aes(x=time3,y=rpheno,group=index),colour="#EEEE00",size=0.5)
  g2 <- g2 + geom_line(aes(x=time,y=pheno,group=index),colour=color,size=0.5)
  g2 <- g2 + geom_point(data=fit.tl1,aes(x=time1,y=lpheno,group=index),fill="#838B83",colour="#838B83",size=1,shape=8)
  g2 <- g2 + geom_point(data=fit.tl1,aes(x=time2,y=gpheno,group=index),fill="#6959CD",colour="#6959CD",size=1,shape=22)
  g2 <- g2 + geom_point(data=fit.tl1,aes(x=time3,y=rpheno,group=index),fill="#EEEE00",colour="#EEEE00",size=1,shape=23)
  #g2 <- g2 + geom_point(aes(x=time,y=pheno,group=index),fill=color,size=1)
  g2 <- g2 + geom_point(aes(x=time,y=tpheno,group=index),fill=color,size=1)
  g2 <- g2 + scale_x_continuous(limits=c(0,63),breaks=seq(3,63,10),labels=seq(3,63,10))
  g2 <- g2 + geom_hline(yintercept=0,linetype=1,size=0.3)
  g2 <- g2 + scale_y_continuous(limits=c(0,100),breaks=seq(0,100,25),labels=seq(0,100,25))
  g2 <- g2 + xlab("Time (day)")+ylab("Taproot length (mm)") + theme_zg()
  #g2 <- g2 + ggtitle("B")
  g2 <- g2 + theme(plot.title = element_text(vjust = 0.5,hjust = -0.08,face="bold",family="sans"))
  
  g <- g1 + g2 + plot_layout(ncol = 1)
  
  tiff("stem-taproot.tif",height=1200,width=900,res=300)
  print(g)
  dev.off()
  
  
}


pairs.fit1 <- function(dat3,dat4,color,init.par,lpar1,gpar1,rpar1,lpar2,gpar2,rpar2){
  
  X1 <- dat3$pheno
  Y1 <- dat4$pheno
  times<-dat3$sample_times
  n<-length(times)
  
  phe_mean<-colMeans(X1,na.rm=TRUE)
  phe_mean1<-colMeans(Y1,na.rm=TRUE)
  a1 <- ind.get_mu(init.par,times,x1=phe_mean[1],x2=phe_mean1[1])
  a2<- ind.get_mu1(init.par[c(1,2,3,5,6,7)],times,x1=phe_mean[1],x2=phe_mean1[1])
  
  l1<-lpar1[1]/(1+lpar1[3]*exp(-lpar1[2]*times))
  g1<-gpar1[1]*exp(-gpar1[3]*exp(-gpar1[2]*times))
  r1<-rpar1[1]*(1-rpar1[3]*exp(-rpar1[2]*times))^(1/(1-rpar1[4]))
  l2<-lpar2[1]/(1+lpar2[3]*exp(-lpar2[2]*times))
  g2<-gpar2[1]*exp(-gpar2[3]*exp(-gpar2[2]*times))
  r2<-rpar2[1]*(1-rpar2[3]*exp(-rpar2[2]*times))^(1/(1-rpar2[4]))
  
  fit.tn <- data.frame(index=rep(1,n),time=times,pheno=a1[1:n],tpheno=phe_mean,lpheno=l1,gpheno=g1,rpheno=r1)
  
  g1 <- ggplot(fit.tn)
  g1 <- g1 + geom_point(aes(x=time,y=tpheno,group=index),fill=color,size=1)
  g1 <- g1 + geom_line(aes(x=time,y=pheno,group=index),colour=color,size=0.5)
  g1 <- g1 + geom_line(aes(x=time,y=lpheno,group=index),colour="#838B83",size=0.5)
  g1 <- g1 + geom_line(aes(x=time,y=gpheno,group=index),colour="#6959CD",size=0.5)
  g1 <- g1 + geom_line(aes(x=time,y=rpheno,group=index),colour="#EEEE00",size=0.5)
  g1 <- g1 + scale_x_continuous(limits=c(13,78),breaks=seq(13,78,13),labels=seq(13,78,13))
  g1 <- g1 + geom_hline(yintercept=0,linetype=1,size=0.3)
  g1 <- g1 + scale_y_continuous(limits=c(0,60),breaks=c(0,seq(0,60,20)),labels=c(0,seq(0,60,20)))
  g1 <- g1 + xlab("Time (day)")+ylab("Root number (count)") + theme_zg()
  #g1 <- g1 + ggtitle("A")
  #g1 <- g1 + theme(plot.title = element_text(vjust = 0.5,hjust = -0.08,face="bold",family="sans"))
  g1 <- g1 + theme(axis.text.x = element_blank(),axis.title.x=element_blank())
  
  
  fit.tl <- data.frame(index=rep(1,n),time=times,pheno=a1[(n+1):(2*n)],tpheno=phe_mean1,tpheno=phe_mean,lpheno=l2,gpheno=g2,rpheno=r2)
  
  g2 <- ggplot(fit.tl)
  g2 <- g2 + geom_point(aes(x=time,y=tpheno,group=index),fill=color,size=1)
  g2 <- g2 + geom_line(aes(x=time,y=pheno,group=index),colour=color,size=0.5)
  g2 <- g2 + geom_line(aes(x=time,y=lpheno,group=index),colour="#838B83",size=0.5)
  g2 <- g2 + geom_line(aes(x=time,y=gpheno,group=index),colour="#6959CD",size=0.5)
  g2 <- g2 + geom_line(aes(x=time,y=rpheno,group=index),colour="#EEEE00",size=0.5)
  g2 <- g2 + scale_x_continuous(limits=c(13,78),breaks=seq(13,78,13),labels=seq(13,78,13))
  g2 <- g2 + geom_hline(yintercept=0,linetype=1,size=0.3)
  g2 <- g2 + scale_y_continuous(limits=c(0,75),breaks=seq(0,75,25),labels=seq(0,75,25))
  g2 <- g2 + xlab("Time (day)")+ylab("Root length (cm)") + theme_zg()
  #g2 <- g2 + ggtitle("B")
  g2 <- g2 + theme(plot.title = element_text(vjust = 0.5,hjust = -0.08,face="bold",family="sans"))
  
  g <- g1 + g2 + plot_layout(ncol = 1)
  
  tiff("rootnumber-rootlength.tif",height=1200,width=900,res=300)
  print(g)
  dev.off()
  
  
}



H2.plot <- function(dat,all_sigqtl,cutree_num,filename="heatmap_qtl.tif"){
  
  allsnp<-dat$snps[all_sigqtl,]
  y<-dat$pheno
  dd <- H2(y,allsnp)
  h2 <- c()
  for(ii in 1:dim(allsnp)[1]){
    h2 <- rbind(h2,cumsum(dd$va[,ii]+dd$vd[,ii])/cumsum(dd$vs)*100)
  }
  rownames(h2)<-paste("Q",all_sigqtl,sep="")
  colnames(h2)<-paste("day",dat$sample_times,sep="")
  #bk <- c(seq(-3,3,by=1.2))
  scale_h2<-t(apply(h2,1,scale))
  hc<-hclust(dist(scale_h2))
  
  tiff(filename=paste("hclust",filename,sep=""),width=400,height=1000,res=300)
  p1<-fviz_dend(hc, k = cutree_num, 
                k_colors = c("#E7B800", "#00AFBB", "#FFAEB9")[1:cutree_num],#c(yellow,green,pink)
                show_labels = FALSE,
                rect = TRUE,
                horiz = T,
                main="",
                xlab="")
  print(p1)
  dev.off()
  
  
  if(cutree_num==2)
    ss<-c(2,1)
  if(cutree_num==3)
    ss<-c(2,1,3)
  
  nh2<-c()
  group<-c()
  gaps<-c()
  for(i in ss){
    pos<-which(cutree(hc, k=cutree_num)==i)
    nh2<-rbind(nh2,h2[pos,])
    group<-c(group,rep(paste("type",i,sep=""),length(pos)))
    gaps<-c( gaps,length(pos))
  }
  group<-data.frame(group)
  rownames(group)<-rownames(nh2)
  if(cutree_num==3)
    colors<-list(group=c(type1="#00AFBB", type2= "#FFAEB9",type3="#E7B800"))
  if(cutree_num==2)
    colors<-list(group=c(type1="#00AFBB",type2= "#E7B800"))
  gaps_row<-cumsum(gaps[-length(gaps)])
  
  tiff(filename=filename,width=1000,height=1000,res=300)
  p<-pheatmap(nh2,scale = "row",display_numbers=F,
              cluster_cols=F,cluster_rows = F,
              show_colnames=T,show_rownames=F,
              cutree_rows=cutree_num,
              annotation_row=group,
              annotation_colors=colors,
              gaps_row = gaps_row,
              #breaks=bk,
              angle_col=45,
              border_color = "NA",
              annotation_legend=FALSE,
              color=viridis_pal(option = "D")(50),
              silent=T)
  #color=colorRampPalette(brewer.pal(11, "PuOr"))(50))
  #color=colorRampPalette(brewer.pal(11, "RdYlBu"))(50))
  #color=colorRampPalette(c("navy", "white", "firebrick3"))(50))
  #color=colorRampPalette(brewer.pal(11, "Spectral"))(50))
  #             color = redgreen(100))
  #             color = bluered(100))
  #             color=colorRampPalette(brewer.pal(9, "OrRd"))(50))
  #color = colorRampPalette(colors = c("blue","white","red"))(100))
  require(ggplotify)
  g <- as.ggplot(p)
  g<- g + theme(plot.margin = margin(t = 1,  r = 1, b = 1,  l = 7)) # c(top,right.bottom,left)
  print(g)
  dev.off()
  
  
  
}


H2.plot11 <- function(dat,all_sigqtl,cutree_num,filename="heatmap_qtl.tif"){
  
  allsnp<-dat$snps[all_sigqtl,]
  y<-dat$pheno
  dd <- H2(y,allsnp)
  h2 <- c()
  for(ii in 1:dim(allsnp)[1]){
    h2 <- rbind(h2,cumsum(dd$va[,ii]+dd$vd[,ii])/cumsum(dd$vs)*100)
  }
  rownames(h2)<-paste("Q",all_sigqtl,sep="")
  colnames(h2)<-paste("day",dat$sample_times,sep="")
  #bk <- c(seq(-3,3,by=1.2))
  scale_h2<-t(apply(h2,1,scale))
  hc<-hclust(dist(scale_h2))
  
  tiff(filename=paste("hclust",filename,sep=""),width=400,height=1000,res=300)
  p1<-fviz_dend(hc, k = cutree_num, 
                k_colors = c("#E7B800", "#00AFBB", "#FFAEB9")[1:cutree_num],#c(yellow,green,pink)
                show_labels = FALSE,
                rect = TRUE,
                horiz = T,
                main="",
                xlab="")
  print(p1)
  dev.off()
  
  
  if(cutree_num==2)
    ss<-c(2,1)
  if(cutree_num==3)
    ss<-c(2,3,1)
  
  
  nh2<-c()
  group<-c()
  gaps<-c()
  for(i in ss){
    pos<-which(cutree(hc, k=cutree_num)==i)
    nh2<-rbind(nh2,h2[pos,])
    group<-c(group,rep(paste("type",i,sep=""),length(pos)))
    gaps<-c( gaps,length(pos))
  }
  group<-data.frame(group)
  rownames(group)<-rownames(nh2)
  if(cutree_num==3)
    colors<-list(group=c(type1="#E7B800", type2="#FFAEB9",type3= "#00AFBB"))
  if(cutree_num==2)
    colors<-list(group=c(type1="#00AFBB",type2= "#E7B800"))
  gaps_row<-cumsum(gaps[-length(gaps)])
  
  tiff(filename=filename,width=1000,height=1000,res=300)
  p<-pheatmap(nh2,scale = "row",display_numbers=F,
              cluster_cols=F,cluster_rows = F,
              show_colnames=T,show_rownames=F,
              cutree_rows=cutree_num,
              annotation_row=group,
              annotation_colors=colors,
              gaps_row = gaps_row,
              #breaks=bk,
              angle_col=45,
              border_color = "NA",
              annotation_legend=FALSE,
              color=viridis_pal(option = "D")(50),
              silent=T)
  require(ggplotify)
  g <- as.ggplot(p)
  g<- g + theme(plot.margin = margin(t = 1,  r = 1, b = 1,  l = 7)) # c(top,right.bottom,left)
  print(g)
  dev.off()
  
  
  
}


H2.plot1 <- function(dat,all_sigqtl,gtitle="a             Stem",cutree_num,filename="heatmap_qtl.tif"){
  
  allsnp<-dat$snps[all_sigqtl,]
  y<-dat$pheno
  dd <- H2(y,allsnp)
  h2 <- c()
  for(ii in 1:dim(allsnp)[1]){
    h2 <- rbind(h2,cumsum(dd$va[,ii]+dd$vd[,ii])/cumsum(dd$vs)*100)
  }
  rownames(h2)<-paste("Q",all_sigqtl,sep="")
  colnames(h2)<-paste("day",dat$sample_times,sep="")
  #bk <- c(seq(-3,3,by=1.2))
  scale_h2<-t(apply(h2,1,scale))
  hc<-hclust(dist(scale_h2))
  
  tiff(filename=paste("hclust",filename,sep=""),width=400,height=1000,res=300)
  p1<-fviz_dend(hc, k = cutree_num, 
                #k_colors = c("#E7B800", "#00AFBB", "#2E9FDF")[1:cutree_num],#c(yellow,green,blue)
                k_colors = c("#E7B800", "#00AFBB", "#FFAEB9")[1:cutree_num],#c(yellow,green,pink)
                show_labels = FALSE,
                rect = TRUE,
                horiz = T,
                main="",
                xlab="")
  print(p1)
  dev.off()
  
  nh2<-c()
  group<-c()
  gaps<-c()
  
  if(cutree_num==2)
    ss<-c(2,1)
  if(cutree_num==3)
    ss<-c(1,3,2)
  
  for(i in ss){
    pos<-which(cutree(hc, k=cutree_num)==i)
    nh2<-rbind(nh2,h2[pos,])
    group<-c(group,rep(paste("type",i,sep=""),length(pos)))
    gaps<-c( gaps,length(pos))
  }
  group<-data.frame(group)
  rownames(group)<-rownames(nh2)
  if(cutree_num==3)
    colors<-list(group=c(type1="#FFAEB9", type2="#E7B800",type3= "#00AFBB"))
  if(cutree_num==2)
    colors<-list(group=c(type1="#E7B800",type2= "#00AFBB"))
  gaps_row<-cumsum(gaps[-length(gaps)])
  
  tiff(filename=filename,width=1000,height=1000,res=300)
  p<-pheatmap(nh2,scale = "row",display_numbers=F,
              cluster_cols=F,cluster_rows = F,
              show_colnames=T,show_rownames=F,
              cutree_rows=cutree_num,
              annotation_row=group,
              annotation_colors=colors,
              gaps_row = gaps_row,
              #breaks=bk,
              angle_col=45,
              border_color = "NA",
              annotation_legend=FALSE,
              color=viridis_pal(option = "D")(50),
              silent=T)
  
  require(ggplotify)
  g <- as.ggplot(p)
  g<- g + theme(plot.margin = margin(t = 1,  r = 1, b = 1,  l = 7)) # c(top,right.bottom,left)
  print(g)
  dev.off()
  
  
}




H2 <- function(y,allsnp){
  
  a <- Impute(t(allsnp), impute.method="bestguess")
  nnsnp <- dim(a)[2]
  
  asnp <- a-1
  dsnp <- 1-abs(asnp)
  colnames(asnp) <- paste("A",1:nnsnp,sep="")
  colnames(dsnp) <- paste("D",1:nnsnp,sep="")
  vs <- c()
  va <- c()
  vd <- c()
  ad.e <- c()
  for(i in 1:dim(y)[2]){
    gsnp <- as.data.frame(cbind(phen=y[,i],asnp,dsnp))
    
    fmla <- as.formula(paste("phen ~ ", paste(colnames(gsnp)[-1], collapse= "+")))
    
    glm.ad <- glm(fmla, family = gaussian(),data=gsnp)
    ad.e <- rbind(ad.e,glm.ad$coefficients[-1]) 
    
    V <- anova(glm.ad)
    snpv <- V$Deviance[-1]
    vs <- c(vs,V$`Resid. Dev`[1])
    va <- rbind(va,snpv[1:nnsnp])
    vd <- rbind(vd,snpv[(nnsnp+1):length(snpv)])
  }
  
  return(list(vs=vs,va=va,vd=vd,ad.e=ad.e))
  
}



Impute<-function(Z, impute.method){
  
  p<-dim(Z)[2]
  
  if(impute.method =="random"){
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX,i])/2
        Z[IDX,i]<-rbinom(length(IDX),2,maf1)
      }
    }
  } else if(impute.method =="fixed"){
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX,i])/2
        Z[IDX,i]<-2 * maf1
      }
    }
  } else if(impute.method =="bestguess") {
    
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX,i])/2
        Z[IDX,i]<-round(2 * maf1)
      }
    }
    
  } else {
    stop("Error: Imputation method shoud be \"fixed\", \"random\" or \"bestguess\" ")
  }
  
  return(Z)
}



pieplot<-function(sigqtl,filename){
  
  data <- data.frame(name=names(table(sigqtl[,2])),percent=as.numeric(table(sigqtl[,2])))
  #lgname<-paste("lg",sort(as.numeric(substr(data1[,1],3,4))),sep="")
  #data<-c()
  #for(i in 1:length(lgname)){
  #  data<-rbind(data,data1[which(data1[,1]==lgname[i]),])
  #}
  #data<-data.frame(data)
  
  pieval <- data$percent/sum(data$percent)
  pielabels <- paste(data$name, "(", data$percent, ")",sep="")
  
  tiff(filename  ,height=1500,width=1500,res=300)
  par(family = "serif")
  pie3D(pieval, radius=0.9,explode=0.1, labels=pielabels, labelcex=0.8)
  dev.off()
}




multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),byrow=TRUE,
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


bar.plot<-function(twotrait=c("Stem","Taproot"),net1,net2){
  
  trait <- c(rep( twotrait[1], 2) , rep(twotrait[2] , 2) )
  link <- rep(c("inhibitation"  , "activation") , 2)
  count1 <- c(as.numeric(table(net1[,3])),as.numeric(table(net2[,3])))
  node<- rep(c("outgoing"  , "incoming") , 2)
  count2 <- c(length(table(net1[,1])),length(table(net1[,2])),length(table(net2[,1])),length(table(net2[,2])))
  data <- data.frame(trait,link,count1,node,count2)
  
  bar1 <- ggplot(data, aes(fill=link, y=count1, x=trait)) +
    geom_bar(position="dodge", stat="identity",width = 0.8)
  bar1 <- bar1 + scale_fill_brewer(palette = 'Pastel1')
  bar1 <- bar1 + theme_classic(base_family="serif")
  bar1 <- bar1 + labs(y = 'Link type') 
  bar1 <- bar1 + xlab(NULL) 
  bar1 <- bar1 + geom_text(aes(label=count1),size=4,
                           position = position_dodge(width = 0.8), 
                           vjust=0.4,hjust=0.4,family="serif")
  bar1 <- bar1 + theme( plot.margin=unit(c(0.1,0.1,0.1,0.1), 'lines'),#c(top,right,bottom,left)
                        axis.title.x = element_text(vjust=0,size=12),
                        axis.title.y = element_text(vjust=-0.2,size=12),
                        axis.text.x=element_text(size=10),
                        axis.text.y=element_text(size=10),
                        legend.text=element_text(size=10),
                        legend.title = element_blank(),
                        #legend.title = element_text(size=8,hjust=0.5),
                        legend.position="bottom",
                        legend.box = "horizontal")
  bar1 <- bar1 + guides(fill  = guide_legend(nrow = 1), byrow = TRUE)
  bar1 <- bar1 + labs(fill="Link type")
  bar1 <- bar1 +  coord_flip() 
  
  
  bar2 <- ggplot(data, aes(fill=node, y=count2, x=trait)) +
    geom_bar(position="dodge", stat="identity",width = 0.8)
  bar2 <- bar2 + scale_fill_brewer(palette = 'Set2')
  bar2 <- bar2 + theme_classic(base_family="serif")
  bar2 <- bar2 + labs(y = 'Node type') 
  bar2 <- bar2 + xlab(NULL) 
  bar2 <- bar2 + geom_text(aes(label=count2),size=4,
                           position = position_dodge(width = 0.8), 
                           vjust=0.4,hjust=0.4,family="serif")
  bar2 <- bar2 + theme( plot.margin=unit(c(0.1,0.1,0.1,0.1), 'lines'),#c(top,right,bottom,left)
                        axis.title.x = element_text(vjust=0,size=10),
                        axis.title.y = element_text(vjust=-0.2,size=10),
                        axis.text.x=element_text(size=10),
                        axis.text.y=element_text(size=10),
                        legend.text=element_text(size=10),
                        #legend.title = element_blank(),
                        legend.title = element_text(size=12,hjust=0.5),
                        legend.position="bottom",
                        legend.box = "horizontal")
  bar2 <- bar2 + guides(fill  = guide_legend(nrow = 1), byrow = TRUE)
  bar2 <- bar2 + labs(fill="Node type")
  bar2 <- bar2 +  coord_flip() 
  
  tiff("barplot.tif",width=2000,height=500,res=300)
  multiplot(plotlist = list(bar1,bar2),cols =2)
  dev.off()
  
}



theme_zg <- function(..., bg='white'){
  require(grid)
  #theme_classic(...,base_family="sans") +
  theme_classic(...,base_family="serif") +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(c(0.1,0.1,0.1,0.1), 'lines'),#c(top,right,bottom,left)
          panel.background=element_rect(fill='transparent', color='black',size=0.5),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),
          plot.title=element_text(size=15,vjust=1),
          #axis.title = element_text(color='black',size=28),
          axis.title.x = element_text(vjust=2.5,size=13),
          axis.title.y = element_text(vjust=-0.2,size=13),
          #axis.ticks.length = unit(-0.4,"lines"),
          axis.text.x=element_text(size=10,vjust=2),
          axis.text.y=element_text(size=10),
          axis.ticks = element_line(color='black'),
          #axis.ticks.margin = unit(0.8,"lines"),
          legend.position='none',
          legend.title=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'),
          strip.background=element_rect(fill='transparent', color='transparent'),
          strip.text=element_text(color='black',vjust=-3,hjust=0.06,size=16))
  
}




#with lengend
theme_zg1 <- function(..., bg='white'){
  require(grid)
  theme_classic(...,base_family="serif") +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(c(0.5,0.5,0.6,1.2), 'lines'),#c(top,right,bottom,left)
          panel.background=element_rect(fill='transparent', color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),
          plot.title=element_text(size=30,vjust=1),
          axis.title = element_text(color='black',size=28),
          axis.title.x = element_text(vjust=0,size=25),
          axis.title.y = element_text(vjust=1.8,size=25),
          #axis.ticks.length = unit(-0.4,"lines"),
          axis.text.x=element_text(size=18),
          axis.text.y=element_text(size=20),
          axis.ticks = element_line(color='black'),
          #axis.ticks.margin = unit(0.8,"lines"),
          legend.key=element_rect(fill='transparent', color='transparent'),
          strip.background=element_rect(fill='transparent', color='transparent'),
          strip.text=element_text(color='black',vjust=-3,hjust=0.06,size=16))
  
}
