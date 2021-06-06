library(matlib)
library(glmnet)
library(rlist)
library(tidyverse)

### Default options
n = 300

#p>=s+d
p = 500
s = 5
d = 50
sigmax = 2
sigmay = 2
sigma=1

###simultaneous_dataset
###Function: creates a list containing two sets of "observables", X and Y, based off the
###          parameters. Statistical matrices are also produced for future use
###INPUTS: a list (a)n, a positive integer; the sample size
###               (b)p, a positive integer; the number of covariates
###               (c)s, a positive integer; the scarcity of Beta
###               (d)d, a positive integer; the density of delta
###               (e)sigmax, a positive double; the natural error in X
###               (f)sigmay, a positive double; the natural error in Y
###               (g)sigma, a positive double; (UNSURE), used for calculating b
###Outputs: a list (a)X, nxp matrix; the p covariates of each subject n
###                (b)Y, nx1 matrix; the response of interest in each subject
###                (c)beta, px1 matrix; the amount each covariate affects Y(scarse)
###                (d)delta, px1 matrix; the amount the response(Y) affects other covariates(dense)
###                (e)e_x, nxp matrix; the error in each measurement in X
###                (f)e_y, nx1 matrix; the error in each measurement in Y
###                (g)gamma, pxp matrix; the covariance matrix between each covariate(default identity matrix)

simultaneous_dataset=function(n,p,s,d,sigmax,sigmay){

  e_x = matrix(rnorm(n*p,mean=0,sd=sigmax),n,p,byrow=TRUE)
  e_y = matrix(rnorm(n,mean=0,sd=sigmay),n,1,byrow=TRUE)
  beta = matrix(rep(c(1,0),times = c(s,p-s)),p,1,byrow = TRUE)
  delta = matrix(rep(c(0,1,0),times=c(s,d,p-d-s)),p,1,byrow=TRUE)

  gamma=diag(p)
  X = (e_y%*%t(delta)+e_x)%*%solve((diag(p)-beta%*%t(delta)))
  Y = X %*% beta+e_y
  temp=t(delta)%*%beta
  #for this case temp will always be 0
  denom=1+temp[1]
  #print(denom)
  return_list=list("X"=X,"Y"=Y,"beta"=beta,"delta"=delta,"e_x"=e_x,"e_y"=e_y,"gamma"=gamma)
  return(return_list)
}

###lin_perturbed_fail
###FUNCTION: tests Lemma 2 and equation 7. Provides
###          further evidence that they are correct and
###          that the Perturbed Linear Model fails
###INPUTS: (1) a list: the same as the output from
###                generate_dataset
###        (2) sigma, a positive double; the error used to find b
###OUTPUTS: a double; gives the fraction of the times which the
###         fit model is greater than p^-0.5

lin_perturbed_fail=function(dset,sigma){
  #get values from dset
  X=dset$X
  Y=dset$Y
  n=dim(X)[1]
  p=dim(X)[2]
  delta=dset$delta
  beta=dset$beta
  gamma=dset$gamma
  e_y=dset$e_y

  #calculated based from paper, unsure of sigma's addition
  x_cov_inv=t((diag(p)-delta%*%t(beta)))  %*%  solve(sigma**2*(delta%*%t(delta)+gamma))  %*%  (diag(p)-delta%*%t(beta))
  b=(sigma**2)*x_cov_inv  %*% solve(diag(p)-delta%*%t(beta)) %*% delta

  #move epsilon-i over so lasso fit is direct
  Y=Y-e_y+X%*%b

  #determine coefficient(Beta+b) estimate
  fit = cv.glmnet(x=X, y=Y)
  beta_b = as.matrix((coef(fit,S =fit$lambda.min)[-1]))

  #theoretical coefficient based off known Beta and b
  beta_b_actual=dset$beta+b

  #calculate how often (Beta+b)[i] is greater than p^-0.5 (Lemma 2)
  error=0
  for(i in 1:p){
    if(beta_b[i]>=sqrt(1/p)){
      error=error+1
    }
  }
  return(error/p)
}



