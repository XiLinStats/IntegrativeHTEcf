#' Integrative analysis of the HTE and confounding functions (IntHTEcf)
#'
#' Implements integrative HTE analysis combing randomized trial and real-world data.
#' @param A is a vector of treatment  ( (n+m) x 1).
#' @param X is a matrix of covariates without intercept  ( (n+m) x p).
#' @param X.hte is a matrix of covariates for the heterogeneous treatment effect function without intercept  ( (n+m) x p1).
#' @param X.cf is a matrix of covariates for the confounding function without intercept  ( (n+m) x p2).
#' @param Y is a vector of outcome  ( (n+m) x 1).
#' @param S is a vector of the binary indicator of belonging to the randomized trial (RT) sample;
#' i.e., 1 if the unit belongs to the RT sample, and 0 otherwise  ( (n+m) x 1).
#' @param nboots is the number of bootstrap samples.
#'
#' @return
#'
#' \itemize{
#'
#' \item \code{est.meta}: the HTE estimator based on a meta analysis
#'
#' \item \code{est.rct}: the HTE estimator using the semiparametric efficient estimation based only on RCT data
#'
#' \item \code{est.int}: the HTE estimator using the semiparametric efficient estimation based on the combined RCT and RWD
#'
#' \item \code{att.meta}: the ATT estimator based on a meta analysis
#'
#' \item \code{att.rct}: the ATT estimator based only on RCT data
#'
#' \item \code{att.int}: the ATT estimator based on the combined RCT and RWD
#'
#' \item \code{ve.meta}: variance estimates of est.meta
#'
#' \item \code{ve.rct}:  variance estimates of est.rct
#'
#' \item \code{ve.int}:  variance estimates of est.int
#'
#' \item \code{ve.att.meta}:  variance estimates of att.meta
#'
#' \item \code{ve.att.rct}:  variance estimates of att.rct
#'
#' \item \code{ve.att.int}:  variance estimates of att.int
#' }
#'
#' @details
#' Details will be provided in the reference paper.
#'
#' @import MASS  rootSolve mgcv stats
#'
#' @references
#'
#'A reference paper will come up soon.
#'
#' @examples
#'
#'library(mgcv)
#'library(MASS)
#'library(rootSolve)
#'set.seed(1)
#'n=500
#'m=1000
#'## RCT: A,X,U
#'A1 <- rbinom(n,1,0.5)
#'X1 <- runif(n,-1,1)
#'U1 <- rnorm(n,1,1)
#'## RWE: A,X,U
#'A2 <- rbinom(m,1,0.5)
#'Sigma0 <- Sigma1 <- matrix(1,2,2)
#'Sigma0[1,2] <- Sigma0[2,1] <- -0.5
#'Sigma1[1,2] <- Sigma1[2,1] <-  0.5
#'XU0 <- mvrnorm(n = m, c(1,0), Sigma0)
#'XU1 <- mvrnorm(n = m, c(1,0), Sigma1)
#'X2 <- XU0[,1]*(1-A2)+XU1[,1]*A2
#'U2 <- XU0[,2]*(1-A2)+XU1[,2]*A2
#'## Combined data
#'A <- c(A1,A2)
#'X <- c(X1,X2)
#'U <- c(U1,U2)
#'S <- c(rep(1,n),rep(0,m))
#'Y <- 1+X+0.5*X^2+0.5*U+0.5*rnorm(n+m)+A*(1+0.75*X^2)
#'X.hte<-cbind(X,X^2)
#'X.cf <-cbind(X)
#'## True parameter values
#'psi<-c(1,0,0.75)
#'att.true<-2.5
#'X<-as.matrix(X)
#'X.hte<-as.matrix(X.hte)
#'X.cf<-as.matrix(X.cf)
#'IntHTEcf( A, X, X.hte, X.cf, Y, S, nboots=50)
#'
#' @export

IntHTEcf <- function( A, X, X.hte, X.cf, Y, S, nboots){
  options(warn=-1)
  X<-as.matrix(X)
  X.hte<-as.matrix(X.hte)
  X.cf <-as.matrix(X.cf)
  p1<-dim(X.hte)[2]+1
  p2<-dim(X.cf)[2]+1

  out<-allestimators(X.hte, X.cf, A, X, Y, S)

  est.meta <- out$est.meta
  att.meta <- out$att.meta
  est.rct <- out$est.rct
  att.rct <- out$att.rct
  est.int <- out$est.int
  att.int <- out$att.int

  if(nboots>0){
    boot.est.meta <- matrix(NA,nboots,length(est.meta))
    boot.att.meta <- matrix(NA,nboots,length(att.meta))
    boot.est.rct  <- matrix(NA,nboots,length(est.rct))
    boot.att.rct  <- matrix(NA,nboots,length(att.rct))
    boot.est.int  <- matrix(NA,nboots,length(est.int))
    boot.att.int  <- matrix(NA,nboots,length(att.int))

    for(kk in 1:nboots){

      sample1<-sample( which(S==1),size= sum(S==1),replace=TRUE)
      sample2<-sample( which(S==0),size= sum(S==0),replace=TRUE)
      bootsample<-c(sample1,sample2)
      A.boot<-A[bootsample]
      X.boot<-X[bootsample,]
      X.hte.boot<-X.hte[bootsample,]
      X.cf.boot<-X.cf[bootsample,]
      S.boot<-S[bootsample]
      Y.boot<-Y[bootsample]

      out<-allestimators(X.hte.boot, X.cf.boot, A.boot, X.boot, Y.boot, S.boot)
      boot.est.meta[kk,] <- out$est.meta
      boot.att.meta[kk,] <- out$att.meta
      boot.est.rct[kk,] <- out$est.rct
      boot.att.rct[kk,] <- out$att.rct
      boot.est.int[kk,] <- out$est.int
      boot.att.int[kk,] <- out$att.int

    }
    ve.meta     <- apply(boot.est.meta,2,var)
    ve.att.meta <- apply(boot.att.meta,2,var)
    ve.rct      <- apply(boot.est.rct,2,var)
    ve.att.rct  <- apply(boot.att.rct,2,var)
    ve.int      <- apply(boot.est.int,2,var)
    ve.att.int  <- apply(boot.att.int,2,var)

    names(ve.meta)<- paste("psi",(1:p1)-1,sep="")
    names(ve.rct )<- paste("psi",(1:p1)-1,sep="")
    names(ve.int )<-c( paste("psi",(1:p1)-1,sep=""), paste("phi",(1:p2)-1,sep=""))
  }

  if(nboots==0){
    ve.meta     <- NA
    ve.att.meta <- NA
    ve.rct      <- NA
    ve.att.rct  <- NA
    ve.int      <- NA
    ve.att.int  <- NA
  }

  return(list(est.meta=est.meta,
              att.meta=att.meta,
              est.rct=est.rct,
              att.rct=att.rct,
              est.int=est.int,
              att.int=att.int,
              ve.meta=ve.meta,
              ve.att.meta=ve.att.meta,
              ve.rct=ve.rct,
              ve.att.rct=ve.att.rct,
              ve.int=ve.int,
              ve.att.int=ve.att.int))
}

####################################################################################
## function list
####################################################################################

#' allestimators
#'
#' A function calculates all estimators.
#' @param X.hte cov
#' @param X.cf  cov
#' @param A treatment
#' @param X cov
#' @param Y outcome
#' @param S sample indicator
#' @return all estimators
#' @import MASS  rootSolve mgcv stats
#' @export

allestimators<-function(X.hte, X.cf, A, X, Y, S){
  loc.rct <- which(S==1)
  loc.rwe <- which(S==0)
  n<-length(loc.rct)
  m<-length(loc.rwe)
  X<-as.matrix(X)
  X.hte<-as.matrix(X.hte)
  X.cf <-as.matrix(X.cf)
  p1<-dim(X.hte)[2]+1
  p2<-dim(X.cf)[2]+1
  ## RCT
  A1<-A[loc.rct]
  X1<-X[loc.rct,]
  X1.hte<-X.hte[loc.rct,]
  X1.bs <-X.cf[loc.rct,]
  Y1<-Y[loc.rct]
  ## RWD
  A2<-A[loc.rwe]
  X2<-X[loc.rwe,]
  X2.hte<-X.hte[loc.rwe,]
  X2.bs <-X.cf[loc.rwe,]
  Y2<-Y[loc.rwe]

  ## propensity score estimation
  gam1<-mgcv::gam(A1~s(X1),family = binomial() )
  eX1<- gam1$fitted.values
  gam2<-mgcv::gam(A2~s(X2),family = binomial() )
  eX2<- gam2$fitted.values
  eX<-c(eX1,eX2)
  eX[loc.rct]<-eX1
  eX[loc.rwe]<-eX2

  ## Meta analysis
  ## ipw-adjusted outcome regression based only on RCT
  Y.adj<-( (Y)*A/eX-(Y)*(1-A)/(1-eX))
  lmx<- X.hte
  jj<-glm(Y.adj~lmx)
  est.meta<-jj$coeff
  names(est.meta)<- paste("psi",(1:p1)-1,sep="")
  est.meta
  psihat <- est.meta
  tauZ <- cbind(1,X.hte)%*%psihat
  att.meta<-sum( (1-S)*tauZ)/sum(1-S)
  att.meta

  ## a preliminary estimator
  psi<-rep(0,dim(X.hte)[2]+1)
  jj<- rootSolve::multiroot(f = pre.ee1,start =psi,X.hte=X.hte,A=A,X=X,Y=Y,S=S)
  est.pre1<-jj$root
  est.pre1

  ## nuisance function estimation
  psihat <- est.pre1
  tau <- cbind(1,X.hte)%*%psihat
  H<- Y-tau*A
  thisy <- H[loc.rct]
  thisx <- X[loc.rct,]
  jj <- mgcv::gam(thisy~s(thisx))
  EHsx1 <-jj$fitted.values
  v1 <- var(jj$residuals)
  thisy <- H[loc.rwe]
  thisx <- X[loc.rwe,]
  jj <- mgcv::gam(thisy~s(thisx))
  EHsx2 <-jj$fitted.values
  v2 <- var(jj$residuals)
  EHsx <- c(EHsx1,EHsx2)
  EHsx[loc.rct]<-EHsx1
  EHsx[loc.rwe]<-EHsx2
  vh <- c(rep(v1,n),rep(v2,m))
  vh[loc.rct]<-v1
  vh[loc.rwe]<-v2

  jj<- rootSolve::multiroot(f = opt.eerct, start = psi,X.hte=X.hte, A=A,X=X,S=S,Y=Y,EHsx=EHsx)
  est.rct<-jj$root
  names(est.rct)<-paste("psi",(1:p1)-1,sep="")
  est.rct
  psihat <- est.rct[1:p1]
  tau <- cbind(1,X.hte)%*%psihat
  att.rct<-sum(S*tau)/sum(S)

  psiphi<-rep(0,dim(X.hte)[2]+dim(X.cf)[2]+2)
  jj<- rootSolve::multiroot(f = prelim.eeint, start = psiphi,
                X.hte=X.hte,X.cf=X.cf,A=A,X=X,S=S,Y=Y,eX=eX,vh=vh)

  est.intpre<-jj$root
  est.intpre

  ## nuisance function estimation

  p1<-dim(X.hte)[2]+1
  p2<-dim(X.cf )[2]+1
  psihat <- est.intpre[1:p1]
  phihat <- est.intpre[(p1+1):(p1+p2)]
  tau <- cbind(1,X.hte)%*%psihat
  lambda <- cbind(1,X.cf)%*%phihat
  H <- Y-A*tau-(1-S)*lambda*(A-eX)

  thisy <- H[loc.rct]
  thisx <- X[loc.rct,]
  jj <- mgcv::gam(thisy~s(thisx))
  EHsx1 <-jj$fitted.values
  v1 <- var(jj$residuals)
  thisy <- H[loc.rwe]
  thisx <- X[loc.rwe,]
  jj <- mgcv::gam(thisy~s(thisx))
  EHsx2 <-jj$fitted.values
  v2 <- var(jj$residuals)
  EHsx <- c(EHsx1,EHsx2)
  EHsx[loc.rct]<-EHsx1
  EHsx[loc.rwe]<-EHsx2
  vh <- c(rep(v1,n),rep(v2,m))
  vh[loc.rct]<-v1
  vh[loc.rwe]<-v2

  jj<- rootSolve::multiroot(f = opt.eeint, start = psiphi,
                X.hte=X.hte,X.cf=X.cf,A=A,X=X,S=S,Y=Y,eX=eX,vh=vh,EHsx=EHsx)

  est.int<-jj$root
  names(est.int)<-c( paste("psi",(1:p1)-1,sep=""), paste("phi",(1:p2)-1,sep=""))
  est.int
  psihat <- est.int[1:p1]
  tau <- cbind(1,X.hte)%*%psihat
  att.int<-sum(S*tau)/sum(S)

  return(list(est.meta=est.meta,
              att.meta=att.meta,
              est.rct=est.rct,
              att.rct=att.rct,
              est.int=est.int,
              att.int=att.int))
}

#' pre.ee1
#'
#' A function calculates pre.ee1.
#' @param par parameter
#' @param X.hte cov
#' @param A treatment
#' @param X cov
#' @param Y outcome
#' @param S sample indicator
#' @return pre.ee1
#' @export

pre.ee1<-function(par,X.hte,A,X,S,Y){
  loc.rct <- which(S==1)
  X<-as.matrix(X)
  X.hte<-as.matrix(X.hte)
  A1 <- A[loc.rct]
  X1 <- X[loc.rct,]
  X1.hte <- X.hte[loc.rct,]
  Y1 <- Y[loc.rct]
  psi <- par
  tau <- cbind(1,X1.hte)%*%psi
  H1 <- Y1-A1*tau
  apply(cbind(1,X1.hte) * matrix(H1*(A1-0.5),sum(S==1),dim(X1.hte)[2]+1,byrow=FALSE ),2,sum)
}

#' opt.eerct
#'
#' A function calculates opt.eerct.
#' @param par parameter
#' @param X.hte cov
#' @param A treatment
#' @param X cov
#' @param Y outcome
#' @param S sample indicator
#' @param H mimicking potential outcome
#' @param EHsx mimicking potential outcome
#' @return opt.eerct
#' @export

opt.eerct<-function(par,X.hte,A,X,S,Y,H,EHsx){
  loc.rct <- which(S==1)
  X<-as.matrix(X)
  X.hte<-as.matrix(X.hte)
  A1 <- A[loc.rct]
  X1 <- X[loc.rct,]
  X1.hte <- X.hte[loc.rct,]
  Y1 <- Y[loc.rct]
  EHsx1 <- EHsx[loc.rct]
  psi <- par
  tau <- cbind(1,X1.hte)%*%psi
  H1 <- Y1-(A1) *tau - EHsx1
  apply( cbind(1,X1.hte) * matrix(H1*(A1-0.5),sum(S==1),dim(X1.hte)[2]+1,byrow=FALSE ),2,sum)
}

#' prelim.eeint
#'
#' A function calculates prelim.eeint.
#' @param par parameter
#' @param A treatment
#' @param X.hte cov
#' @param X.cf cov
#' @param X cov
#' @param Y outcome
#' @param S sample indicator
#' @param eX propensity score
#' @param H mimicking potential outcome
#' @param vh var estimates for two samples
#' @return prelim.eeint
#' @export

prelim.eeint<-function(par,A,X.hte,X.cf,X,S,Y,H,eX,vh){
  X.hte<-as.matrix(X.hte)
  X.cf <-as.matrix(X.cf )
  X    <-as.matrix(X    )
  p1<-dim(X.hte)[2]+1
  p2<-dim(X.cf )[2]+1
  psi <- par[1:p1]
  phi <- par[(p1+1):(p1+p2)]
  tau <- cbind(1,X.hte)%*%psi
  lambda <- cbind(1,X.cf)%*%phi
  H <- (Y-(A) *tau -(1-S)*lambda*(A-eX))/vh
  c(apply(cbind(1,X.hte)*matrix(H*(A-eX),length(S),dim(X.hte)[2]+1,byrow=FALSE ),2,sum),
    apply(cbind(1,X.cf)*matrix( (1-S)*H*(A-eX),length(S),dim(X.cf)[2]+1,byrow=FALSE ),2,sum))
}

#' opt.eeint
#'
#' A function calculates opt.eeint.
#' @param par parameter
#' @param A treatment
#' @param X.hte cov
#' @param X.cf cov
#' @param X cov
#' @param Y outcome
#' @param S sample indicator
#' @param eX propensity score
#' @param H mimicking potential outcome
#' @param EHsx mimicking potential outcome
#' @param vh var estimates for two samples
#' @return opt.eeint
#' @export

opt.eeint<-function(par,A,X.hte,X.cf,X,S,Y,H,eX,EHsx,vh){
  X.hte<-as.matrix(X.hte)
  X.cf <-as.matrix(X.cf )
  X    <-as.matrix(X    )
  p1<-dim(X.hte)[2]+1
  p2<-dim(X.cf )[2]+1
  psi <- par[1:p1]
  phi <- par[(p1+1):(p1+p2)]
  tau <- cbind(1,X.hte)%*%psi
  lambda <- cbind(1,X.cf)%*%phi
  H <- (Y-(A) *tau -(1-S)*lambda*(A-eX) - EHsx)/vh
  c(apply(cbind(1,X.hte)*matrix(H*(A-eX),length(S),dim(X.hte)[2]+1,byrow=FALSE ),2,sum),
    apply(cbind(1,X.cf)*matrix( (1-S)*H*(A-eX),length(S),dim(X.cf)[2]+1,byrow=FALSE ),2,sum))
}
