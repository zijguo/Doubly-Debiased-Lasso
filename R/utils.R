

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
  tau = quantile(D, rho)
  Dtilde = pmin(D, tau)
  Q = diag(nrow(X)) - U %*% diag(1 - Dtilde / D) %*%  t(U)
  Xtilde = Q %*% X
  Ytilde = Q %*% Y

  #perform a lasso fit for trimmed QX and QY(STEP 2)
  fit = glmnet::cv.glmnet(x=Xtilde, y=Ytilde)

  #B-init
  betahat = as.matrix((coef(fit,S =fit$lambda.min)[-1]))

  #determine bhat
  res = Y-X %*% betahat
  fit = glmnet::cv.glmnet(x=U %*% diag(D), y=res, alpha=1)
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
###        active_set_scaling: the correction for the size of the active set, scale estimate by n/(n-k)
###
### OUTPUT: a double, a point approximation for the true noise error in Y
estimate_sigma = function(X,Y,rho=0.5,alt=FALSE,active_set_scaling=FALSE){

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
  # The length of D has length min(n,p). In order to avoid the over-shrinkage to the last n-p eigenvector
  # when n>p, the following is necessary.
  Q =diag(nrow(X))-U %*% diag(1-Dtilde / D) %*%  t(U)
  Xtilde = Q %*% X
  Ytilde = Q %*% Y

  #Step 7, two methods of sigmahat computation, the first has shown more robust results
  divisor=sum(diag(Q%*%Q))
  error=(norm(Ytilde-Xtilde%*%betahat,type='2'))^2
  sigmahat=(error/divisor)^0.5
  if(active_set_scaling){
    sigmahat = (nrow(X)/(nrow(X)-min(nrow(X)/2, Matrix::nnzero(betahat))))^0.5*sigmahat
  }
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
###        index, an integer; the index of X(the subject) which should be used for projection
###
### OUTPUT: z, a nx1 matrix; the projection direction vector
find_z = function(X,index){

  n = dim(X)[1]
  p = dim(X)[2]

  #Xj and X-j
  X_j = X[,index]
  X_negj = X[,-index]

  #regress X-j on xj, use least min lambda to estimate gamma(Step 4)
  cvfit = glmnet::cv.glmnet(x=X_negj, y=X_j)
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

