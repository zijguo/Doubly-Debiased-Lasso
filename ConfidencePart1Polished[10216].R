library(matlib)
library(glmnet)
library(rlist)
library(tidyverse)

### Default options
n = 300
p = 500
q = 3
s = 5
sigma = 2
sigmaE = 2
pert = 1

### generate_dataset
### FUNCTION: Creates a list containing a single dataset
###           representing a sample with many covariates and confounding variables
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
generate_dataset = function(n, p, s, q, sigmaE, sigma, pert){
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
### INPUT: a list (a) X, an nxp matrix; represents the measurements for each subject and covariate
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
###
### OUTPUT: a list (a) betahat, an unbiased estimator for Beta-j
###                (b) bhat, an unbiased estimator for b
estimate_coefficients = function(dset){
  X = dset$X
  Y = dset$Y

  #perform Single Val Decomp on X
  UDV_list = svd(X)
  U = UDV_list$u
  D = UDV_list$d
  V = t(UDV_list$v)

  #forms QX and QY
  tau = median(D)
  Dtilde = pmin(D, tau)
  F = U %*% diag(Dtilde / D) %*%  t(U)
  Xtilde = F %*% X
  Ytilde = F %*% Y

  #perform a lasso fit for trimmed QX and QY
  fit = cv.glmnet(x=Xtilde, y=Ytilde)
  betahat = as.matrix((coef(fit,S =fit$lambda.min)[-1]))

  res = Y-X %*% betahat
  #determine B-init hat
  fit = cv.glmnet(x=U %*% diag(D), y=res, alpha=1)

  bhat = t(V) %*% coef(fit, s=fit$lambda.min)[-1]

  return_listcoeff = list("betahat" = betahat, "bhat"=bhat)
  return(return_listcoeff)
}

### estimate_sigma
### FUNCTION: Creates a list containing a single dataset
###           representing a sample with many covariates and confounding variables
###           along with
###
### INPUT: a list (a) X, an nxp matrix; represents the measurements for each subject and covariate
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
###
### OUTPUT: a double, a point approximation for the true standard error in Y
estimate_sigma = function(dset, betahat, bhat){
  #uses both fitting estimators to create an unbiased estimator of sigma
  residuals = dset$Y - dset$X %*% betahat - dset$X %*% bhat
  sigmahat = mean(residuals^2)^0.5

  return (sigmahat)
}

### find_z
### FUNCTION: computes the projection direction used to construct CI
###
###
###
### INPUT: X, a nxp matrix with the sample observations
###        idx, the index of X(the subject) which should be used for projection
###
### OUTPUT: z, a nx1 matrix; the projection direction vector
find_z = function(X, idx){
  #determines parameters
  n = dim(X)[1]
  p = dim(X)[2]

  #Xj and X-j
  Y = X[,idx]
  X = X[,-idx]

  #regress Xj on x-j, use least min lambda to estimate gamma
  cvfit = cv.glmnet(x=X, y=Y )
  gamma = coef(cvfit, s=cvfit$lambda.min)[-1]

  #eq. (8)
  z = n^(-0.5)*(Y - X %*% gamma)

  #variation from eq. (23) with 25% increase(read 3.6)
  V = 1.25 * n^0.5 * norm(z, type ="F") / (t(z) %*% Y)

  #The Z with the highest magnitude will be used
  for (lam in cvfit$glmnet.fit$lambda){
    gamma = coef(cvfit,s =lam)[-1]
    z = n^(-0.5)*(Y - X %*% matrix(gamma,p-1,1))
    if (n^0.5 * (norm(z)/(t(z) %*% Y)) > V){
      break
    }
  }
  #normalize Z with Frobenius norm
  z = z/norm(z,type = "F")
  return(z)
}
### CI
### FUNCTION: computes the (1-alpha)%-CI for Bj
###
###
###
### INPUT: a list (a) X, an nxp matrix; represents the measurements for each subject and covariate
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
###        idx, a positive integer, the index, j, to be used for finding Beta-j
###        alpha, a double, the significance level of the confidence interval, default is 0.05
###
### OUTPUT: a list (a) a vector of length one with the lower bound of the interval as a double
###                (b) a vector of length one with the upper bound of the interval as a double
CI = function(dset, idx, alpha=0.05){
  #determines parameters
  X = dset$X
  Y = dset$Y
  n = dim(X)[1]
  p = dim(X)[2]

  #single value decomposition of X (Trim transform)
  UDV_list = svd(X)
  U = UDV_list$u
  D = UDV_list$d
  V = t(UDV_list$v)

  #Reduce D, then trim X and Y(QX and QY)
  Dtilde = pmin(D, median(D))
  F = U %*% diag(Dtilde / D) %*% t(U)
  X = F %*% X
  Y = F %*% Y

  #determine projection direction then estimate betahat and bhat
  z = find_z(X, idx)
  est = estimate_coefficients(dset)
  betahat = est$betahat
  bhat = est$bhat

  # eq. (12) for a point estimation of Bj
  dblasso = betahat[idx] + t(z) %*% (Y - X %*% betahat) / (t(z) %*% X[, idx])

  #eq. (23) for
  V = n^0.5 * norm(F %*% z, type ="F") / (t(z) %*% X[, idx])
  sigmahat = estimate_sigma(dset, betahat, bhat)

  #2-sided CI critical value for alpha
  quantile = qnorm(1-alpha/2,mean=0,sd = 1,lower.tail=TRUE)

  #determine bounds of interval, from eq. (13)
  lower = dblasso - quantile * sigmahat * V * n^(-0.5)
  upper = dblasso + quantile * sigmahat * V * n^(-0.5)

  return_list = list("lower" =as.vector(lower), "upper" = as.vector(upper))
  return(return_list)
}

### FUNCTION: for a list of datasets, it computes the width,
###           coverage and bias of the corresponding CIs
###
### INPUT: A list of N lists each with the objects:
###                (a) X, an nxp matrix; represents the measurements for each subject and covariate
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
###
### OUTPUT: a list, (a) an atomic vector; the lower bound of the confidence interval of Bj
###                 (b) an atomic vector; the upper bound of the confidence interval of Bj
test_CI = function(dsets){
  #define variables
  width = 0
  coverage = 0
  bias = 0
  idx = 1
  N = length(dsets)
  for (dset in dsets){
    #determine average width and bias among all confidence intervals
    ci = CI(dset, idx)
    width = width + (ci$upper - ci$lower)/N
    bias = bias + ((ci$upper + ci$lower)/2 - dset$'beta'[idx])/N

    #determine if the the true value of beta lies within each interval(coverage)
    if(ci$lower < dset$'beta'[idx] & ci$upper> dset$'beta'[idx]){
      coverage = coverage + 1/N
    }
  }
  return_list = list("width" =width, "coverage" = coverage, "bias" = bias)
  return(return_list)
}
### compute_bias
### FUNCTION: Creates a list containing a single dataset
###           representing a sample with many covariates and confounding variables
###           along with
###
### INPUT: A list of N lists each with the objects:
###                (a) X, an nxp matrix; represents the measurements for each subject and covariate
###               (b) Y, an nx1 matrix; the measured response for each subject
###               (c) epsilon, an nx1 matrix; total error in Y
###               (d) beta, a px1 matrix; sparse coefficient vector(how X affects Y independently from H)
###               (e) b, a px1 matrix; perturbation vector induced by hidden confounding
###               (f) H, an nxq matrix; the amount of confounding for each subject
###               (g) gamma, a qxp matrix; encodes the linear effect of Hi on the covariates Xi
###               (h) E, a nxp matrix; the amount which X is independent from H
###               (i) nu, a nx1 matrix; the error in Y independent from X and H
###               (j) delta, a qx1 matrix; encodes the linear effect of Hi on the response Yi
###               (k) sigma, a 1x1 matrix; SD of epsilon
###        idx, the index, j, of interest
###
### OUTPUT: a list, (a) a double, average bias term for Beta
###                 (b) a double, average bias term for b
###                 (c) a double, average bias term for b if oracle value is known
###                 (d) a double, average bias term for b if bound
###                 (e) a double, average bias term for b if corrected

#computes the projection direction used to construct CI
#dsets should be in the form of list
compute_bias = function(dsets){
  #define variables
  ret_B_beta = 0
  ret_B_oracle = 0
  ret_B_bound = 0
  ret_B_b = 0
  ret_B_corrected = 0
  N = length(dsets)

  for (dset in dsets){
    #determine parameters from dset
    X = dset$X
    Y = dset$Y
    n = dim(X)[1]
    p = dim(X)[2]

    #performs SVD on X
    UDV_list = svd(dset$X)
    U = UDV_list$u
    D = UDV_list$d
    V = t(UDV_list$v)

    #trim X and Y (QX and QY)
    Dtilde = pmin(D, median(D))
    F = U %*% diag(Dtilde / D) %*%  t(U)
    X = F %*% X
    Y = F %*% Y

    idx = 1
    z = find_z(X, idx)
    est = estimate_coefficients(dset)
    betahat = est$betahat
    bhat = est$bhat

    #determine the error between approximation and true value
    row_del = dset$beta - betahat

    #defined bottom of page 19. bias terms for beta and b
    B_beta =c(t(z) %*% X[,-idx] %*% row_del[-idx,]/ norm(F %*% z, type = "F"))
    B_b = c((t(z) %*% X %*% dset$b) /norm(F %*% z, type = "F"))

    #alternate measurements of bias for b term dependent on situation
    B_b_oraclebound = norm(t(z) %*% X, type = "F") / norm(F %*% z, type = "F") * norm(dset$b, type = "F")
    B_b_bound = norm(t(z) %*% X, type = "F") / norm(F %*% z, type = "F") * log(p) * p^-0.5
    B_b_corrected = c((t(z) %*% X %*% (dset$b - bhat)) /norm(F %*% z, type = "F"))

    #averages each of the statistics over each data set
    ret_B_beta = ret_B_beta + abs(B_beta) / N
    ret_B_b = ret_B_b+abs(B_b) / N
    ret_B_oracle =ret_B_oracle+ abs(B_b_oraclebound) / N
    ret_B_bound =ret_B_bound+ abs(B_b_bound) / N
    ret_B_corrected = ret_B_corrected+abs(B_b_corrected) / N

  }
  return_list = list("B_beta"=ret_B_beta, "B_b" = ret_B_b, "B_oracle" =ret_B_oracle,"B_bound" =ret_B_bound,
                     "B_corrected"=ret_B_corrected)
  return(return_list)
}

