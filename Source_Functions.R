library(glmnet)
library(rlist)
library(tidyverse)

### generate_dataset
### FUNCTION: Creates a list containing a single dataset
###           simulating a sample with many covariates and confounding variables
###           along with corresponding statistical values
###
### INPUT: n, a positive integer; the number of subjects in the sample
###        p, a positive integer; the number of covariates(measurements) per subject
###        s, a positive integer; the sparsity of the sample(s<<(sqrt(n)/log(p)))
###        q, a positive integer; the number of confounding variables in the sample
###        sigmaE, a positive double; the SD of the error epsilon-i
###        sigma, a positive double; the SD of Y independent of X(natural error)
###        pert, a positive double; tuning value for the perturbation(intensity of confounding)
###
### OUTPUT: a list (a) X, an nxp matrix; represents the measurements for each subject and covariate
###               (b) Y, an nx1 matrix; the measured response for each subject
###               (c) epsilon, an nx1 matrix; total error in Y
###               (d) beta, a px1 matrix; sparse coefficient vector(how X affects Y independently from H)
###               (e) b, a px1 matrix; perturbation vector induced by hidden confounding
###               (f) H, an nxq matrix; the amount of confounding for each subject
###               (g) gamma, a qxp matrix; encodes the linear effect of Hi on the covariates Xi
###               (h) E, a nxp matrix; the amount which X is independent from H(not confounded)
###               (i) nu, a nx1 matrix; the error in Y independent from X and H
###               (j) delta, a qx1 matrix; encodes the linear effect of Hi on the response Yi
###               (k) sigma, a 1x1 matrix; SD of epsilon
generate_dataset = function(n=300, p=500, s=5, q=3, sigmaE=2, sigma=2, pert=1){
  #create H and Gamma with N(0,1) values and of appropriate size. H can be tuned with pert
  H = pert*matrix(rnorm(n*q,mean=0,sd=1),n,q,byrow = TRUE)
  Gamma = matrix(rnorm(q*p,mean=0,sd=1),q,p,byrow = TRUE)
  #value of X independent from H
  E = matrix(rnorm(n*p,mean=0,sd=sigmaE),n,p,byrow = TRUE)

  #defined in eq. (2), high-dimensional measured covariates
  X = E + H %*% Gamma

  delta = matrix(rnorm(q*1,mean=0,sd=1),q,1,byrow = TRUE)

  #px1 matrix, creates beta with 1s in the first s entries and the remaining p-s as 0s
  beta = matrix(rep(c(1,0),times = c(s,p-s)),p,1,byrow = TRUE)

  #nx1 matrix with values of mean 0 and SD of sigma, error in Y independent of X
  nu = matrix(rnorm(n*1,mean=0,sd=sigma),n,1,byrow = TRUE)

  #eq. (1), the response of the Structural Equation Model
  Y = X %*% beta + H %*% delta + nu

  #eq. (3), defines perturbation vector
  b = pert^2 * solve(pert^2*t(Gamma) %*% Gamma + sigmaE^2 * diag(p)) %*% t(Gamma) %*% delta

  #eq. (3), from hidden confounding model
  epsilon = H %*% delta - X %*% b + nu

  sigma = (sigma^2 + pert^2 * t(delta)%*% delta -pert^2 * t(b) %*% t(Gamma) %*% delta)^0.5

  return_list = list("X"= X,"Y"= Y, "epsilon"= epsilon,"beta"= beta,"b"= b,"H"= H,"Gamma"= Gamma,"E"= E,"nu"= nu,"delta"= delta,"sigma"= sigma)
  return(return_list)
}

### estimate_coefficients
### FUNCTION: using a single data set, creates unbiased estimators
###           for betahat and bhat using single value decomposition
###           and a lasso fit on trimmed X and Y for betahat then
###           a lasso fit on U and the residuals for bhat
###
### INPUT:  X, an nxp matrix; represents the measurements for each subject and covariate
###         Y, an nx1 matrix; the measured response for each subject
###         rho, a double; the trim parameter for the Q transform, default is 0.5(median)
###
### OUTPUT: a list (a) betahat, an unbiased estimator for Beta_init
###                (b) bhat, an unbiased estimator for b
estimate_coefficients = function(X,Y,rho=0.5){
  #perform Single Val Decomp on X
  UDV_list = svd(X)
  U = UDV_list$u
  D = UDV_list$d
  V = t(UDV_list$v)

  #forms QX and QY(STEP 1)
  tau = quantile(D,rho)
  Dtilde = pmin(D, tau)
  Q = U %*% diag(Dtilde / D) %*%  t(U)
  Xtilde = Q %*% X
  Ytilde = Q %*% Y

  #perform a lasso fit for trimmed QX and QY(STEP 2)
  fit = cv.glmnet(x=Xtilde, y=Ytilde)

  #B-init
  betahat = as.matrix((coef(fit,S =fit$lambda.min)[-1]))

  #determine bhat
  res = Y-X %*% betahat
  fit = cv.glmnet(x=U %*% diag(D), y=res, alpha=1)
  bhat = t(V) %*% coef(fit, s=fit$lambda.min)[-1]

  return_listcoeff = list("betahat" = betahat,"bhat"=bhat)
  return(return_listcoeff)
}

### estimate_sigma
### FUNCTION: Creates an unbiased estimator of noise level based on either the default method
###           of the proposed estimator. If alt is chosen to be true, the sigmahat is found using
###           the residuals of normal lasso regression
###
### INPUT: X, an nxp matrix; represents the measurements for each subject and covariate
###        Y, an nx1 matrix; the measured response for each subject
###        rho, a double; the trim parameter for the Q transform, default is 0.5(median)
###        alt, a boolean; determines which method to use for computing sigmahat, default false
###
### OUTPUT: a double, a point approximation for the true noise error in Y
estimate_sigma = function(X,Y,rho=0.5,alt=FALSE){

  #uses both fitting estimators to create an unbiased estimator of sigma
  est_coef=estimate_coefficients(X,Y,rho)
  betahat=est_coef$betahat
  bhat=est_coef$bhat

  UDV_list = svd(X)
  U = UDV_list$u
  D = UDV_list$d
  V = t(UDV_list$v)

  tau = quantile(D,rho)
  Dtilde = pmin(D, tau)
  Q = U %*% diag(Dtilde / D) %*%  t(U)
  Xtilde = Q %*% X
  Ytilde = Q %*% Y

  #Step 7, two methods of sigmahat computation, the first has shown more robust results
  divisor=sum(diag(Q%*%Q))
  error=(norm(Ytilde-Xtilde%*%betahat,type='2'))^2
  sigmahat=(error/divisor)^0.5

  if(alt){
    residuals = Y - X %*% betahat - X %*% bhat
    sigmahat = mean(residuals^2)^0.5
  }

  return (sigmahat)
}

### find_z
### FUNCTION: computes the normalized projection direction used to construct CI;
###           P transform must be applies to X before performing this function
###
### INPUT: X, a nxp matrix with the sample observations
###        idx, an integer; the index of X(the subject) which should be used for projection
###
### OUTPUT: z, a nx1 matrix; the projection direction vector
find_z = function(X,idx){

  n = dim(X)[1]
  p = dim(X)[2]

  #Xj and X-j
  X_j = X[,idx]
  X_negj = X[,-idx]

  #regress X-j on xj, use least min lambda to estimate gamma(Step 4)
  cvfit = cv.glmnet(x=X_negj, y=X_j)
  gamma = coef(cvfit, s=cvfit$lambda.min)[-1]
  #eq. (8), residuals(Step 5)
  z = n^-0.5 *(X_j - X_negj %*% gamma)

  #variation from eq. (23) with 25% increase(read 3.6)
  V = 1.25*n^0.5*norm(z, type ="2") / (t(z) %*% X_j)

   #take first z whose variance is at most 25% larger than for the CV lambda
   for (lam in cvfit$glmnet.fit$lambda){
     gamma = coef(cvfit,s =lam)[-1]
     z = n^(-0.5)*(X_j - X_negj %*% matrix(gamma,p-1,1))
    if (n^0.5 * (norm(z,type="2")/(t(z) %*% X_j)) > V){
       break
     }
   }

  #normalize Z with 2 norm
  z = z/norm(z,type = "2")
  return(z)
}

### CI
### FUNCTION: computes the (1-alpha)%-CI for Bj, the variance of the estimator, and 
###           if the true values of beta and b are known, the bias
###
### INPUT: X, an nxp matrix; represents the measurements for each subject and covariate
###        Y, an nx1 matrix; the measured response for each subject
###        idx, a positive integer; the index, j, to be used for finding Beta-j
###        alpha, a double; the significance level of the confidence interval, default is 0.05
###        rho, a double; the trim level for the Q transform
###        rhop, a double; the trim level for the P transform
###
### OUTPUT: a list (a) point, the doubly debiased lasso estimator
###                (b) lower, lower bound of the interval
###                (b) upper, upper bound of the interval
###                (c) Variance, a double; the variance of the dblasso estimator
CI = function(X,Y,idx,alpha=0.05,rho=0.5,rhop=0.5){
  #determines parameters
  n = dim(X)[1]
  p = dim(X)[2]

  X_negj=X[,-idx]

  #single value decomposition of X (Trim transform)
  UDV_list = svd(X_negj)
  U = UDV_list$u
  D = UDV_list$d
  V = t(UDV_list$v)

  #Reduce D, then make P trim
  Dtilde = pmin(D, quantile(D,rhop))
  P = U %*% diag(Dtilde / D) %*% t(U)
  P_X = P %*% X

  #determine projection direction then estimate betahat and bhat
  z = find_z(P_X, idx)
  est = estimate_coefficients(X,Y,rho)
  betahat = est$betahat
  bhat = est$bhat

  # eq. (12) for a point estimation of Bj(Step 6)
  dblasso=t(z) %*% P%*%(Y - X_negj %*% betahat[-idx])/(t(z)%*%P%*%X[,idx])

  #eq. (23) for
  Variance=(t(z)%*%(P^2)%*%z)/(t(z)%*%P%*%X[,idx])^2

  #(Step 7)
  sigmahat = estimate_sigma(X,Y,rho)
  stddev = sigmahat*Variance^.5

  #2-sided CI critical value for alpha
  quantile = qnorm(1-alpha/2,mean=0,sd = 1,lower.tail=TRUE)
  #determine bounds of interval, from eq. (13) (Step 8)
  lower = dblasso - quantile * stddev
  upper = dblasso + quantile * stddev

  
  return_list = list("point"= dblasso, "lower" =lower, "upper" = upper, "Variance"=stddev^2)
  return(return_list)
}
