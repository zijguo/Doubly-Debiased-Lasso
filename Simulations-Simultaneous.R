library(matlib)
library(glmnet)
library(rlist)
library(tidyverse)

#####test expectation values
idx = 1
n = 3000
p = 500
s = 5
d = 10
sigmax = 1
sigmay = 1
N=100
e_y = matrix(rnorm(n,mean=0,sd=sigmay),n,1,byrow=TRUE)
beta = matrix(rep(c(1,0),times = c(s,p-s)),p,1,byrow = TRUE)
delta = matrix(rep(c(0,1,0),times=c(s,d,p-d-s)),p,1,byrow=TRUE)
gamma_test1=0
gamma_test2=0
gamma_test_total=0
variance1=0
variance2=0
variance3=0
for(i in 1:N){
  e_x = matrix(rnorm(n*p,mean=0,sd=sigmax),n,p,byrow=TRUE)
  X=e_x
  Y=X%*%beta+e_y
  X=Y%*%t(delta)+e_x
  gamma_test1=gamma_test1+X[,idx][-1]%*%t(X[,idx][-1])/N #(p-1)x(p-1)
  gamma_test2=gamma_test2+t(X[idx][1]%*%X[,idx][-1])/N #(p-1)x1
  variance1=variance1+t(X[idx][1]%*%X[,idx][-1])/N     #(p-1)x1
  variance2=variance2+X[idx][1]%*%X[,idx][-1]/N        #1x(p-1)
  variance3=variance3+(X[,idx][-1])%*%t(X[,idx][-1])/N #(p-1)x(p-1)
}
#X[,idx][-1]   (p-1)x1
#print(variance2%*%solve(variance3)%*%variance1)
print(gamma_test1%*%gamma_test2)
#####

#####test perturbed linear model approach
n=30000
p=1000
s=5
d=995
e_x = matrix(rnorm(n*p,mean=0,sd=sigmax),n,p,byrow=TRUE)
e_y = matrix(rnorm(n,mean=0,sd=sigmay),n,1,byrow=TRUE)
beta = matrix(rep(c(1,0),times = c(s,p-s)),p,1,byrow = TRUE)
delta = matrix(rep(c(0,1,0),times=c(s,d,p-d-s)),p,1,byrow=TRUE)
gamma = diag(p)
X=e_x
Y=X%*%beta+e_y
X=Y%*%t(delta)+e_x
Yalt=Y-e_y
fitalt=cv.glmnet(x=X,y=Yalt)
beta_alt=as.matrix((coef(fit,s =fit$lambda.min)[-1]))

x_cov_inv=t((diag(p)-delta%*%t(beta)))  %*%  solve(sigma**2*(delta%*%t(delta)+gamma))  %*%  (diag(p)-delta%*%t(beta))
b=(sigma**2)*x_cov_inv  %*% solve(diag(p)-delta%*%t(beta)) %*% delta
Yfit=Y-e_y+X%*%b
fit = cv.glmnet(x=X, y=Yfit)
beta_b = as.matrix((coef(fit,s =fit$lambda.min)[-1]))

beta_b_actual=beta+b
print(norm(beta_b_actual,type='2'))
print(norm(beta_b,type='2'))
print(p^-0.5)
#####



