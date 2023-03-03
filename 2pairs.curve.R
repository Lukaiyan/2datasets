pairs.H0 <- function(dat,dat1,parin2){
  
  allpheno <- cbind(dat$pheno,dat1$pheno)
  times <- dat$sample_times
  n<-length(times)
  
  mpheno <- as.numeric(colMeans(allpheno, na.rm=TRUE))
  #if(mpheno[1]==0)
  #  mpheno[1]<-0.1
  loop_k <- 1
  max_iter <- 100
  epsi <- 10^-4
  max_err <- 1
  
  while(loop_k<max_iter && max_err>epsi){
    
    
    oldpar <-c(parin2);
    mle.covar1 <- function(npar){
      
      nnpar <- c(npar,parin2[6:13])
      AA <- curve.mle(par=nnpar,y=allpheno,time.std =times,x1=mpheno[1],x2=mpheno[n+1])
      AA
    }
    r1.covar <- optim(parin2[1:5],mle.covar1,method = "BFGS",control=list(maxit=32000))
    new.covar1 <- r1.covar$par
    #cat("new.coavr1:",unlist( new.covar1), "\n");
    
    mle.1 <- function(npar){
      
      nnpar <- c(new.covar1,npar)
      AA <- curve.mle(nnpar,y=allpheno,time.std =times,x1=mpheno[1],x2=mpheno[n+1])
      AA
    }
    r1 <- optim(c(parin2[6:13]),mle.1,method = "BFGS",control=list(maxit=32000))
    #r1 <- try(optim(c(parin2[6:13]),mle.1,method = "BFGS",control=list(maxit=32000)),TRUE)
    #if (class(r1) == "try-error") 
    #  r1 <- optim(c(parin2[6:13]),mle.1,method = "Nelder-Mead",control=list(maxit=2000))
    
    new1 <- r1$par
  
    
    nparin <- c(new.covar1,new1)
    
    newpar <- c(nparin)
    #cat("newpar:", newpar, "\n");
    
    max_err <- max(abs( oldpar - newpar) );
    
    parin2 <- nparin
    #cat(loop_k, "max.err=", max_err, "allpar", newpar,"\n");
    loop_k <- loop_k+1; 
  }
  LL <- curve.mle(parin2,y=allpheno,time.std =times,x1=mpheno[1],x2=mpheno[n+1])
  return(c(LL,parin2))
}



curve.mle <-function( par,y,time.std,x1,x2)
{
  len.cov <- 5
  par.covar <- par[1:len.cov]
  #sig.inv3 <- SAD3.get_inv_mat(par.covar,time.std, 2)
  sigma <- SAD3.get_mat(par.covar,time.std, 2)
  #solve( sig.inv3 )
  
  #if(x1==0){
  #  x1<-0.1
  #}
  
  curve.par <- par[(len.cov+1):(len.cov+8)]
  mu <- ind.get_mu(curve.par,time.std,x1,x2)
  
  yy <- y
  fy <- dmvnorm(yy,mu,sigma)
  #fy[which(fy<=.Machine$double.eps)] <- .Machine$double.eps
  A <- -sum(log(fy))
  #cat("LL=",A,"\n")
  return(A)
}
