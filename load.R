load.data <- function(datt,pheno,nstart,nlen,times){
 
 
  mergedat <- fg_read_data(geno.dat=datt,phe=pheno,nstart, nlen, times)
  
  final <- get_geno_code(d.gen=mergedat$gen,d.phe=mergedat$phe,d.ids=mergedat$ids,nstart,nlen)
  
  info<- mergedat$snp.info
  nobs <- mergedat$n.obs
  gene_table <- matrix(final$gen0,nrow=nobs,byrow=T) 
  colnames(gene_table) <- mergedat$snp.info[,1]
  
  snp <- t(gene_table)
  new <- final$phe0
  new[which(is.na(new))]<-0
  dattt <- list(pheno=new[,1:length(times)],snps=snp,sample_times=times,snpinfo=info)
  return(dattt)
  
}



fg_read_data<-function(geno.dat,phe, nstart, nlen, times)
{
  tb.gen<-geno.dat
  tb.ids<-colnames(geno.dat);
  tb.ids<-tb.ids[-c(1:3,546:557)]
  tb.ids<-as.numeric(substring(tb.ids,2,4))
  tb2.ids<-phe[,1]
  tb.idx<-c()
  for(i in 1:length(tb2.ids)){
    tb.idx<-c(tb.idx,which(tb.ids==tb2.ids[i]))
  }
  
  if(length(tb.idx)>0){
    tb.gen<-tb.gen[,c(1:3,tb.idx+3)]
  }

  tb.phe<-phe

  nstop<-(nstart+nlen-1);
  if (nstop > dim(tb.gen)[1])
    nstop <- dim(tb.gen)[1];
  
  if (nstart==0 && nlen==0)
  {
    nstop <- dim(tb.gen)[1];
    nstart<- 1;
  }		
  
  n.obs <- dim(tb.phe)[1];
  times <- seq(1,57,4); 
  
  snp.info <- tb.gen[c(nstart:nstop),c(1:3)];
  gen   <- tb.gen[c(nstart:nstop), c(4:(n.obs+3))]
  n.snp <- dim(gen)[1];
  
  return(list(ids=tb.phe[,1], n.obs=n.obs, n.snp=n.snp, times=times, 
              snp.info=snp.info, gen=gen, phe=tb.phe[,-1] ))
}

get_geno_code<-function( d.gen, d.phe, d.ids,nstart,nlen)
{
  snpall<-matrix(NA,nrow=nlen,ncol=2)
  for(i in nstart:(nstart+nlen-1)){
    character <- paste(unlist(d.gen[i,]),collapse="")
    name <- names(table(substring(character,1:20,1:20)))
    type <- name[which(name!="-")]
    snpall[i,]<- type
    
  }
  d.gen2 <- as.character(unlist(d.gen));
  
  snpB <- as.character(snpall[,1]);
  snpA <- as.character(snpall[,2]);
  
  QQ2<- paste(snpB, snpB, sep=""); 
  qq0<- paste(snpA, snpA, sep="");
  Qq1<- list(paste(snpA, snpB, sep=""), paste( snpB, snpA, sep="") ) ;
  
  d.g <- array( NA, length(d.ids) );
  d.g[which(d.gen2==QQ2)]<-2;
  d.g[which(d.gen2==qq0)]<-0;
  d.g[which(d.gen2==c(Qq1[[1]]))]<-1;
  d.g[which(d.gen2==c(Qq1[[2]]))]<-1;
  
  return(list(gen0=d.g, phe0=d.phe, ids0=d.ids))
}
