# DDL
The goal of DDL is to implement the Doubly Debiased Lasso estimator proposed in  <https://arxiv.org/abs/2004.03758>. It Computes the Doubly Debiased Lasso estimator of a single regression coefficient in the high-dimensional linear model with hidden confounders and also constructs the confidence interval forthe target regression coefficient.

## Installation 
You can install the released version of DDL from CRAN with:
```R
install.packages("DDL")
```
## Example
This is a basic example which shows you how to solve a common problem

Generate the data:
```R
idx = c(1,2,10)
n = 2000
p = 200
s = 5
q = 3
sigmaE = 2
sigma = 2
pert = 1

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
```
Call 'DDL'
```R
result = DDL(X,Y,idx)
result
#> $idx
#> [1]  1  2 10
#> 
#> $beta0
#> [1] 0.8044596 0.8948753 0.0000000
#> 
#> $point
#> [1]  0.94125533  1.14088157 -0.08988171
#> 
#> $se
#> [1] 0.09172868 0.08951032 0.08310032
#> 
#> $CI
#>           lower      upper
#> [1,]  0.7614704 1.12104024
#> [2,]  0.9654446 1.31631858
#> [3,] -0.2527553 0.07299194
#> 
#> attr(,"class")
#> [1] "DDL"
```
'summary' method for 'DDL'
```R
summary(result)
#> Estimators: 
#>  id est_spectral_deconfounding     est_ddl Std. Error   z value     Pr(>|z|)
#>   1                  0.8044596  0.94125533 0.09172868 10.261298 1.052819e-24
#>   2                  0.8948753  1.14088157 0.08951032 12.745810 3.289565e-37
#>  10                  0.0000000 -0.08988171 0.08310032 -1.081605 2.794282e-01
#> 
#>  Confidence Intervals: 
#>  id      lower      upper
#>   1  0.7614704 1.12104024
#>   2  0.9654446 1.31631858
#>  10 -0.2527553 0.07299194
'''
