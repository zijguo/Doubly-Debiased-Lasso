#' Summarizing DDL
#' @description 'summary' method for class 'DDL'
#' @param object An object of class 'DDL'
#' @param ... Ignored
#' @return The function 'summary.DDL' returns a list of summary statistics of DDL given 'DDL'
#' @export
#' @example
#' index = 1
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
#' result = DDL(X, Y, index)
#'
#' summary(result)
summary.DDL = function(object,...){
  est_ddl = object$est_ddl
  se = object$se
  index = object$index
  
  n.loading = length(se)
  output.est = data.frame(matrix(NA,nrow = n.loading, ncol = 5))
  output.est[,1] = index
  output.est[,2] = est_ddl
  output.est[,3] = se
  output.est[,4] = est_ddl / se
  output.est[,5] = 2 * pnorm(-abs(output.est[,4]))
  colnames(output.est) = c("Index", "est_ddl","Std. Error","z value", "Pr(>|z|)")
  obj = list(output.est = output.est)
  class(obj) = "summary.DDL"
  obj
}
#' Summarizing DDL
#' @description 'summary' method for class 'DDL'
#' @param x An object of class 'summary.DDL'
#' @param ... Ignored
#' @export
print.summary.DDL = function(x,...){
  cat("Call: \n Estimation and Inference for the index Coefficient \n\n")
  print(x$output.est, quote = F, row.names = F)
}


ci = function(x, alpha = 0.05, alternative = c("two.sided","less","greater")){
  UseMethod("ci")
}
ci.default = function(x, alpha = 0.05, alternative = c("two.sided","less","greater")){
  cat("This is a method of class DDL")
}
#' Computing confidence intervals
#' @description 'ci' method for class 'DDL'
#' @param x An object of class 'DDL'
#' @param alpha alpha Level of significance to construct confidence interval
#' @param alternative indicates the alternative hypothesis to construct confidence interval and must be one of "two.sided" (default), "less", or "greater".
#' @export
#' @example
#' index = 1
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
#' result = DDL(X, Y, index)
#' # default alpha is 0.05
#' ci(result, alpha = 0.05)
#' ci(result, alpha = 0.05, alternative = "less")
#' ci(result, alpha = 0.05, alternative = "greater")
ci.DDL = function(x, alpha = 0.05, alternative = c("two.sided","less","greater")){
  alternative = match.arg(alternative)
  se = x$se
  est_ddl = x$est_ddl
  index = x$index
  if(alternative=="two.sided"){
    output.ci = cbind(est_ddl - qnorm(1 - alpha / 2)*se, est_ddl + qnorm(1 - alpha / 2)*se)
    output.ci = data.frame(cbind(index, output.ci))
  }else if(alternative=="less"){
    output.ci = cbind(-Inf, est_ddl + qnorm(1 - alpha)*se)
    output.ci = data.frame(cbind(index, output.ci))
  }else if(alternative=="greater"){
    output.ci = cbind(est_ddl - qnorm(1 - alpha)*se, Inf)
    output.ci = data.frame(cbind(index, output.ci))
  }
  colnames(output.ci) = c("index","lower","upper")
  rownames(output.ci) = paste(index, sep = " ")
  cat("Call: \n Confidence Intervals Construction for the index coefficient \n\n")
  cat("Confidence Intervals: \n")
  print(output.ci, quote = F, row.names =F)
  invisible(x)
}
