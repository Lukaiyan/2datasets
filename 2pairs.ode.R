
COMP.f1 <-function( parameters, state, times ){
  Lorenz<-function(t, state, parameters) 
  {
    with( as.list(c(state, parameters)),
          {
            dX <- r1*X*(1-(X/K1)^s1)
            dY <- r2*Y*(1-(Y/K2)^s2)
            
            list(c(dX, dY))
          }
          
    ) # end with(as.list ...
  }
  
  out <- try(rk(y = state, times = times, func = Lorenz, parms = parameters,method="rk4") );
  out;
}


COMP.f <-function( parameters, state, times ){
  Lorenz<-function(t, state, parameters) 
  {
    with( as.list(c(state, parameters)),
          {

            #2
            dX <- r1*X*(1-(X/K1)^s1)+r1*X*a/K1*Y
            dY <- r2*Y*(1-(Y/K2)^s2)+r2*Y*b/K2*X
            list(c(dX, dY))
          }
          
    ) # end with(as.list ...
  }
  
  out <- try(rk(y = state, times = times, func = Lorenz, parms = parameters,method="rk4") );
  out;
}


ind.mle <- function(s.par,s.y,s.t,x1,x2){
  
  A <- sum((s.y - ind.get_mu(s.par,s.t,x1,x2))^2 )
  
  A
}


ind.get_mu <- function(par, times,x1,x2)
{
  par0 <- par;
  if (class(par0)!="list")
  {
    par0 <- list(
      K1 = par[1],
      r1 = par[2],
      s1 = par[3],
      a = par[4],
      K2 = par[5],
      r2 = par[6],
      s2 = par[7],
      b = par[8]);
  }
  
  state0 <- c(X=x1, Y=x2);
  y <- COMP.f( par0, state0, times );
  

  index.time <- 1:length(times)
  return ( c(y[index.time, 2:3] ) );
}


ind.get_mu1 <- function(par, times,x1,x2)
{
  par0 <- par;
  if (class(par0)!="list")
  {
    par0 <- list(
      K1 = par[1],
      r1 = par[2],
      s1 = par[3],
      K2 = par[4],
      r2 = par[5],
      s2 = par[6]);
  }
  
  state0 <- c(X=x1, Y=x2);
  y <- COMP.f1( par0, state0, times );
  
  
  index.time <- 1:length(times)
  return ( c(y[index.time, 2:3] ) );
}
