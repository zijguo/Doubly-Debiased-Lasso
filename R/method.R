#' Summarizing DDL
#' @description 'summary' method for class 'DDL'
#' @param object An object of class 'DDL'
#' @param ... Ignored
#' @return The function 'summary.DDL' returns a list of summary statistics of DDL given 'DDL'
#' @export
#' @example
#' idx = 1
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
#'
#' summary(result)
summary.DDL = function(object,...){
  point = object$point
  se = object$se
  CI = object$CI
  beta0 = object$beta0
  idx = object$idx
  lower = (object$CI)[,1]
  upper = (object$CI)[,2]


  n.loading = length(se)
  output.est = data.frame(matrix(NA,nrow = n.loading, ncol = 5))
  output.est[,1] = idx
  output.est[,2] = beta0
  output.est[,3] = point
  output.est[,4] = se
  output.est[,5] = point / se
  output.est[,6] = 2 * pnorm(-abs(output.est[,5]))
  colnames(output.est) = c("id","est_spectral_deconfounding","est_ddl","Std. Error","z value", "Pr(>|z|)")
  output.ci = cbind(lower, upper)
  output.ci = data.frame(cbind(idx, output.ci))
  colnames(output.ci) = c("id","lower","upper")
  rownames(output.ci) = paste(idx, sep = " ")
  obj = list(output.est = output.est, output.ci = output.ci)
  class(obj) = "summary.DDL"
  obj
}
#' Summarizing DDL
#' @description 'summary' method for class 'DDL'
#' @param x An object of class 'summary.DDL'
#' @param ... Ignored
#' @export
print.summary.DDL = function(x,...){
  cat("Call: \n Estimation and Inference for the idx Coefficient \n\n")
  cat("Estimators: \n")
  print(x$output.est, quote = F, row.names = F)
  cat("\n Confidence Intervals: \n")
  print(x$output.ci, quote = F, row.names =F)
  invisible(x)
}
