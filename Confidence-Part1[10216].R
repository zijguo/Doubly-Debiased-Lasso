library(matlib)
library(glmnet)
library(rlist)
library(tidyverse)


generate_dataset = function(n, p, s, q, sigmaE, sigma, pert){ #sample size n, measurements per sample p, confounding variables q, mean error sigmaE, error sigma,
  H = pert*matrix(rnorm(n*q,mean=0,sd=1),n,q,byrow = TRUE) #n*q data values are made with mean 0, sd 1, and placed in n x q matrix. These are then multiplied by pert, some amplitude for each confounder
  Gamma = matrix(rnorm(q*p,mean=0,sd=1),q,p,byrow = TRUE)#qxp matrix same vals as above
  E = matrix(rnorm(n*p,mean=0,sd=sigmaE),n,p,byrow = TRUE)#nxp matrix with values of sd sigmaE
  X = E + H %*% Gamma #defined in equation (2)

  delta = matrix(rnorm(q*1,mean=0,sd=1),q,1,byrow = TRUE) #qx1 matrix
  beta = matrix(rep(c(1,0),times = c(s,p-s)),p,1,byrow = TRUE)#s 1s and p-s 0s in nx1 matrix
  nu = matrix(rnorm(n*1,mean=0,sd=sigma),n,1,byrow = TRUE)#nx1 matrix

  Y = X %*% beta + H %*% delta + nu #defined in equation (1)
  b = pert^2 * solve(pert^2*t(Gamma) %*% Gamma + sigmaE^2 * diag(p)) %*% t(Gamma) %*% delta #from equation (3), does solve need second argument?

  epsilon = H %*% delta - X %*% b + nu #equation (3)
  sigma = (sigma^2 + pert^2 * t(delta)%*% delta -pert^2 * t(b) %*% t(Gamma) %*% delta)^0.5

  return_list = list("X"= X,"Y"= Y, "epsilon"= epsilon,"beta"= beta,"b"= b,"H"= H,"Gamma"= Gamma,"E"= E,
                     "nu"= nu,"delta"= delta,"sigma"= sigma)
 return(return_list)
}

## Default options
n = 300
p = 500
q = 3
s = 5
sigma = 2
sigmaE = 2
pert = 1
default_options = list("beta:oracle"=TRUE,
                       "beta:transformation"=TRUE,
                       "b:transformation"=FALSE,
                       "sigma:oracle"=FALSE,
                       "z:peter"=TRUE,
                       "CI:transformation"=TRUE,
                       "bias:add_b"=FALSE,
                       "bias:widen"=FALSE,
                       "bias:oracle"=FALSE)


options_description = "
beta:oracle = choose penalty lambda for which the estimate on the path is the closest
beta:transformation = do we apply trim transformation for estimating beta
b:transformation = if true, we estimate b from the residuals by using L1 penalty in
                   the basis of singular vectors of X
                   if false, we estimate b from the residuals by ridge
sigma:oracle = do not estimate sigma from the data but use the true value
z:peter = in order to get the projection direction, use peter's hdi method with
          largest lambda with at most 25% higher variance than cv, cv used otherwise
CI:transformation = trim variables in order to construct CIs
bias:add_b = estimate b in addition to beta and use it to produce residuals from
             which we get CI
bias:widen = widen the CI by an upper bound of the bias term induced by b
bias:oracle = widen the CI by the oracle value of the bias term induced by b
"


#estimates beta and b in the given dataset

estimate_coefficients = function(dset, options){
  X = dset$X
  Y = dset$Y

  UDV_list = svd(X) #spectral decomposition(forming Q)
  U = UDV_list$u
  D = UDV_list$d
  V = t(UDV_list$v)
  tau = median(D)
  Dtilde = pmin(D, tau) #Dtilde=D except all elements greater than median(D) become the median.
  F = U %*% diag(Dtilde / D) %*%  t(U)
  Xtilde = F %*% X #both X and Y are reduced
  Ytilde = F%*%Y

  if(options$'beta:oracle'){
    if(options$'beta:transformation'){
      fit = glmnet(x=Xtilde, y=Ytilde, intercept = FALSE)} #lasso regression with tilde
    else{fit = glmnet(x=X, y=Y, intercept = FALSE)}#lasso regression w/o tilde
    value = Inf
    for(lam in fit$lambda){
      bh = coef(fit, s= lam) #for each possible lambda model, determine beta hat, the coef of fit
      bh = matrix(bh[-1],dim(dset$beta)[1],dim(dset$beta)[2], byrow = TRUE) #removed intercept from the coeffs #check
      if (value > norm(dset$beta - bh, type = "1")) {
        value = norm(dset$beta - bh, type = "1")#determine which beta hat has the least error(lowest norm with beta from data set)
        betahat = bh #beta hat is estimated to be the best fit among lambda from lasso
    }
    }
  }else{
    if(options$'beta:transformation'){
      fit = cv.glmnet(x=Xtilde, y=Ytilde)}
    else{fit = cv.glmnet(x=X, y=Y)}
    betahat = as.matrix((coef(fit,S =fit$lambda.min)[-1]))}
  if( options$'b:transformation'){
    res = Y-X%*%betahat    #residual from betahat
    fit = cv.glmnet(x=U%*%diag(D), y=res, alpha=1)
    bhat = t(V)%*%coef(fit, s=fit$lambda.min)[-1] #calculating bhat from the residuals, what equation is this from?
  }else{res = Y - X%*%betahat
  fit = cv.glmnet(x=X, y=res, alpha=0)
  bhat = coef(fit,s=fit$lambda.min)[-1]}
  return_listcoeff = list("betahat" = betahat, "bhat"=bhat)
  return(return_listcoeff)
}

#estimates sigma in the given dataset


estimate_sigma = function(dset, betahat, bhat, options){
  residuals = dset$Y - dset$X %*% betahat - dset$X %*% bhat
  sigmahat = mean(residuals^2)^0.5 #definition of sigma hat, should it be divided by n-1 instead of n though?
  if (options$'sigma:oracle'){
    sigmahat = dset$sigma
  }
  return (sigmahat)
}

#computes the projection direction used to construct CI

find_z = function(X, idx, options){
  n = dim(X)[1] #X is nxp matrix, so this finds n from X
  p = dim(X)[2] #same as above^, but second finds p from X
  Y = X[,idx]
  X = X[,-idx]

  cvfit = cv.glmnet(x=X, y=Y )
  gamma = coef(cvfit, s=cvfit$lambda.min)[-1]
  z = n^(-0.5)*(Y - X%*%gamma) #equation (8)

  if (options$'z:peter'){
    #take first z whose variance is at most 25% larger than for the CV lambda
    V = 1.25 * n^0.5 * norm(z, type ="F") / (t(z) %*% Y)
    for (lam in cvfit$glmnet.fit$lambda){
      gamma = coef(cvfit,s =lam)[-1] #equation (9)
      z = n^(-0.5)*(Y - X%*%matrix(gamma,p-1,1))
      if (n^0.5 * (norm(z)/(t(z) %*% Y)) > V){break}
    }
  }
  z = z/norm(z,type = "F") #normalize Z with frobenius norm
  return(z)
}

#computes the (1-alpha)%-CI for the variable indexed by idx

CI = function(dset, idx, options, alpha=0.05){
  X = dset$X
  Y = dset$Y
  n = dim(X)[1]
  p = dim(X)[2]

  UDV_list = svd(X)
  U = UDV_list$u
  D = UDV_list$d
  V = t(UDV_list$v)
  if (options$'CI:transformation'){
    Dtilde = pmin(D, median(D))
  }else{Dtilde = D}

  F = U %*% diag(Dtilde / D) %*%  t(U)
  X = F %*% X
  Y = F%*%Y

  z = find_z(X, idx, options)
  est = estimate_coefficients(dset, options)  #we use original variables
  betahat = est$betahat
  bhat = est$bhat

  if(options$'bias:add_b'){
    dblasso = betahat[idx] + t(z) %*% (Y - X%*%betahat - X%*%bhat)/ (t(z) %*% X[, idx]) #using equation (11)
  }else{
    dblasso = betahat[idx] + t(z) %*% (Y - X%*%betahat) / (t(z) %*% X[, idx])#using equation (12)
  }

  V = n^0.5 * norm(F %*% z, type ="F") / (t(z) %*% X[, idx]) #based from equation (23)

  sigmahat = estimate_sigma(dset, betahat, bhat, options)

  quantile = qnorm(1-alpha/2,mean=0,sd = 1,lower.tail=TRUE) #critical value
  lower = dblasso - quantile * sigmahat * V * n^(-0.5)
  upper = dblasso + quantile * sigmahat * V * n^(-0.5) #determines bounds

  if (options$'bias:widen'){
    bias_bound = norm(t(z) %*% X, type = "F")/(t(z) %*% X[, idx]) * (log(p) / p^0.5) #where is this from?
    if (options$'bias:oracle'){
      bias_bound = (t(z) %*% X %*% dset$b) / (t(z) %*% X[, idx])
    }
    lower = lower + bias_bound
    upper = upper - bias_bound #width is decreased due to removal of bias
  }

  if (runif(1) < 0.05){
    print(paste("betahat:", betahat[idx], "dblasso:", dblasso, "V:", V,
                "sigmahat:", sigmahat, "sigma:", dset['sigma'], "lo:", lower,
                "hi:", upper))
  }

  return_list = list("lower" =as.vector(lower), "upper" = as.vector(upper))
  return(return_list)
}

#for a list of datasets, it computes the width, coverage and bias of the corresponding CIs

test_CI = function(dsets, options){
  width = 0
  accuracy = 0
  bias = 0
  idx = 1 #arbitrary choice
  N = length(dsets)

  for (dset in dsets){
    ci = CI(dset, idx, options)
    width = width + (ci$upper - ci$lower)/N
    bias = bias + ((ci$upper + ci$lower)/2 - dset$'beta'[idx])/N #both the width and bias are averaged over all data sets

    if(ci$lower < dset$'beta'[idx] & ci$upper> dset$'beta'[idx]){ #for each data set determine if beta is within the calculated interval
      accuracy = accuracy + 1/N #will give a val between 0 and 1 for the proportion of intervals which are accurate, we expect this to be 1-alpha in the limit
    }
  }
  return_list = list("width" =width, "accuracy" = accuracy, "bias" = bias)
  return(return_list)
}

#Estimating b
#Estimating b well helps decreasing the corresponding bias term B_b.

#We can do it from the residuals by using ridge or Lasso in the basis of singular vectors of X

N = 30
ratio = 0

# for (i in seq(30)){
#   dset = generate_dataset(n, p, s, q, sigmaE, sigma, pert)
#   options = options
#   options$'b:transformation' = TRUE
#   est = estimate_coefficients(dset, options)
#   betahat = est$betahat
#   bhat = est$bhat
#
#   UDV_list = svd(dset$X)
#   U = UDV_list$u
#   D = UDV_list$d
#   V = t(UDV_list$v) #check
#   Dtilde = pmin(D, median(D))
#   F = U %*% diag(Dtilde / D) %*%  t(U)
#
#   if(runif(1) < 0.1){
#     print(paste('betahat',betahat[1], betahat[2], betahat[3]))
#     print(paste(bhat[1]/dset$'b'[1], bhat[2]/dset$'b'[2], bhat[3]/dset$'b'[3],bhat[10]/dset$'b'[10], bhat[11]/dset$'b'[11], bhat[12]/dset$'b'[12]))
#     print(norm(F%*%dset$'X'%*%(bhat - dset$'b'), type = "F") / norm(F%*%dset$'X'%*%dset$'b', type = "F"))
#     print(norm(dset$'X'%*%(bhat - dset$'b'),type ="F") /norm(dset$'X'%*%dset$'b', type = "F"))
#   }
#   ratio = ratio + norm(F%*%dset$'X'%*%(bhat - dset$'b'), type = "F") / norm(F%*%dset$'X'%*%dset$'b')/N
# }


#dsets should be in the form of list
compute_bias = function(dsets, options){
  ret_B_beta = 0
  ret_B_oracle = 0
  ret_B_bound = 0
  ret_B_b = 0
  ret_B_corrected = 0
  N = length(dsets)
  for (dset in dsets){
    X = dset$X
    Y = dset$Y
    n = dim(X)[1]
    p = dim(X)[2]

    UDV_list = svd(dset$X) #single value decomposition, non negative diagonal matrix S and two unitary matrices s.t. U*S*V=X
    U = UDV_list$u
    D = UDV_list$d
    V = t(UDV_list$v)

    if (options$'CI:transformation'){
      Dtilde = pmin(D, median(D))
    }else{Dtilde = D}

    F = U %*% diag(Dtilde / D) %*%  t(U)
    X = F %*% X
    Y = F%*%Y

    idx = 1 #1? all the time
    z = find_z(X, idx, options)

    est = estimate_coefficients(dset, options)  #we use original variables
    betahat = est$betahat
    bhat = est$bhat

    row_del = dset$beta - betahat
    B_beta =c(t(z) %*% X[,-idx] %*% row_del[-idx,]/ norm(F%*%z, type = "F"))

    B_b = c((t(z) %*% X %*% dset$b) /norm(F %*% z, type = "F"))

    B_b_oraclebound = norm(t(z) %*% X, type = "F") / norm(F %*% z, type = "F") * norm(dset$b, type = "F")

    B_b_bound = norm(t(z) %*% X, type = "F") / norm(F %*% z, type = "F") * log(p) * p^-0.5

    B_b_corrected = c((t(z) %*% X %*% (dset$b - bhat)) /norm(F %*% z, type = "F"))

    if (runif(1) < 0.1){
      print(paste('B_beta:', B_beta, 'B_b', B_b, 'Fz',norm(F %*% z, type = "F") / norm(z, type = "F")))
      print(paste('bound', B_b_bound, 'oracle_bound', B_b_oraclebound, 'B_b',B_b, 'corrected', B_b_corrected))
    }
    ret_B_beta = ret_B_beta + abs(B_beta) / N
    ret_B_b = ret_B_b+abs(B_b) / N
    ret_B_oracle =ret_B_oracle+ abs(B_b_oraclebound) / N
    ret_B_bound =ret_B_bound+ abs(B_b_bound) / N
    ret_B_corrected = ret_B_corrected+abs(B_b_corrected) / N    #all of these biases are averaged over the data sets

  }
  return_list = list("B_beta"=ret_B_beta, "B_b" = ret_B_b, "B_oracle" =ret_B_oracle,"B_bound" =ret_B_bound,
                     "B_corrected"=ret_B_corrected)
  return(return_list)
}

