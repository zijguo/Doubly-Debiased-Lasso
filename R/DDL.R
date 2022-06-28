#' Point estimation and inference for a single regression coefficient in the high-dimensional linear model with hidden confounders.
#'
#' Computes the Doubly Debiased Lasso estimator of a single regression coefficient in the high-dimensional linear model with hidden confounders. It also constructs the confidence interval for the target regression coefficient.
#'
#' @param X the covariates matrix, of dimension \eqn{n\times p}
#' @param Y the outcome vector, of length \eqn{n}
#' @param idx the vector of indexes for the regression coefficient of interest
#' @param alpha the significance level of the confidence interval, default is 0.05
#' @param rho the trim level for \eqn{X}, default is 0.5
#' @param rhop the trim level for \eqn{X_{-j}}, default is 0.5
#' @return
#' \item{idx}{the vector of indexes for the regression coefficient of interest}
#' \item{point}{The vector of the Doubly Debiased Lasso estimator of the target regression coefficient}
#' \item{se}{The vector of the standard error of the Doubly Debiased Lasso estimator}
#' \item{CI}{The matrix of the confidence interval for the target regression coefficient}
#' \item{beta0}{The vector of the spectral deconfounding estimator of the whole regression vector}
#' @examples
#' idx = c(1,2,10)
#' n=100
#' p=200
#' s=5
#' q=3
#' sigmaE=2
#' sigma=2
#' pert=1
#'
#' H = pert*matrix(rnorm(n*q,mean=0,sd=1),n,q,byrow = TRUE)
#' Gamma = matrix(rnorm(q*p,mean=0,sd=1),q,p,byrow = TRUE)
#' #value of X independent from H
#' E = matrix(rnorm(n*p,mean=0,sd=sigmaE),n,p,byrow = TRUE)
#'
#' #defined in eq. (2), high-dimensional measured covariates
#' X = E + H %*% Gamma
#'
#' delta = matrix(rnorm(q*1,mean=0,sd=1),q,1,byrow = TRUE)
#'
#' #px1 matrix, creates beta with 1s in the first s entries and the remaining p-s as 0s
#' beta = matrix(rep(c(1,0),times = c(s,p-s)),p,1,byrow = TRUE)
#'
#' #nx1 matrix with values of mean 0 and SD of sigma, error in Y independent of X
#' nu = matrix(rnorm(n*1,mean=0,sd=sigma),n,1,byrow = TRUE)
#'
#' #eq. (1), the response of the Structural Equation Model
#' Y = X %*% beta + H %*% delta + nu
#'
#' result = DDL(X, Y, idx)
#' summary(result)
#' @export
#' @import stats
DDL = function(X,Y,idx,alpha=0.05,rho=0.5,rhop=0.5){
  #determines parameters
  n = dim(X)[1]
  p = dim(X)[2]
  dblasso = rep(NA,length(idx))
  stddev = rep(NA,length(idx))
  lower = rep(NA,length(idx))
  upper = rep(NA,length(idx))
  est = estimate_coefficients(X,Y,rho)
  betahat = est$betahat
  for (i in seq(length(idx))){
    X_negj=X[,-idx[i]]

    #single value decomposition of X (Trim transform)
    UDV_list = svd(X_negj)
    U = UDV_list$u
    D = UDV_list$d
    V = t(UDV_list$v)

    #Reduce D, then make P trim
    Dtilde = pmin(D, quantile(D,rhop))
    P = diag(nrow(X_negj)) - U %*% diag(1 - Dtilde / D) %*% t(U)
    P_X = P %*% X

    #determine projection direction then estimate betahat and bhat(the z here has been multiplied with P)
    z = find_z(P_X, idx[i])
    # bhat = est$bhat

    # eq. (12) for a point estimation of Bj(Step 6)
    dblasso[i] = t(z) %*% P %*% (Y - X_negj %*% betahat[-idx[i]]) / (t(z) %*% P %*% X[,idx[i]])

    #eq. (23) for
    Variance = (t(z) %*% (P^2) %*% z)/(t(z) %*% P %*% X[,idx[i]]) ^ 2

    #(Step 7)
    sigmahat = estimate_sigma(X,Y,rho)
    stddev[i] = sigmahat*Variance^.5

    #2-sided CI critical value for alpha
    quantile = qnorm(1-alpha/2,mean=0,sd = 1,lower.tail=TRUE)
    #determine bounds of interval, from eq. (13) (Step 8)
    lower[i] = dblasso[i] - quantile * stddev[i]
    upper[i] = dblasso[i] + quantile * stddev[i]
    # B_b=t(z) %*% P_X %*% b/(t(z) %*% P_X[,idx]*stddev)
    # B_beta= t(z) %*% P_X[,-idx] %*% (beta[-idx]-betahat[-idx])/(t(z) %*% P_X[,idx]*stddev)
  }
  obj = list(idx = idx,
             beta0 = betahat[idx],
             point= dblasso,
             se = stddev,
             CI = cbind(lower,upper))
  class(obj) = "DDL"
  obj
}


