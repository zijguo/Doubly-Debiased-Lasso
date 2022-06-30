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
index = c(1,2,10)
n = 500
p = 300
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
result = DDL(X,Y,index)
result
#> $index
#> [1]  1  2 10
#> 
#> $est_init
#> [1] 0.7801164 0.8235082 0.0000000
#> 
#> $est_ddl
#> [1]  0.92549403  1.00838136 -0.03282794
#> 
#> $se
#> [1] 0.05655909 0.05814778 0.05246100
#> 
#> attr(,"class")
#> [1] "DDL"
```
'summary' method for 'DDL'
```R
summary(result)
#> Call: 
#>  Estimation and Inference for the index Coefficient 
#> 
#>  Index  est_ddl    Std. Error    z value     Pr(>|z|)
#>      1  0.92549403 0.05655909 16.3633131 3.495805e-60
#>      2  1.00838136 0.05814778 17.3417002 2.278949e-67
#>     10 -0.03282794 0.05246100 -0.6257589 5.314731e-01
```
'ci' method for 'DDL'
```R
#default alpha is 0.05
ci(result, alpha = 0.05)
#> Call: 
#>  Confidence Intervals Construction for the index coefficient 
#> 
#> Confidence Intervals: 
#>  index      lower      upper
#>      1  0.8146403 1.03634780
#>      2  0.8944138 1.12234892
#>     10 -0.1356496 0.06999374

ci(result, alpha = 0.05, alternative = "less")
#> Call: 
#>  Confidence Intervals Construction for the index coefficient 
#> 
#> Confidence Intervals: 
#>  index lower      upper
#>      1  -Inf 1.01852545
#>      2  -Inf 1.10402595
#>     10  -Inf 0.05346273

ci(result, alpha = 0.05, alternative = "greater")
#> Call: 
#>  Confidence Intervals Construction for the index coefficient 
#> 
#> Confidence Intervals: 
#>  index      lower upper
#>      1  0.8324626   Inf
#>      2  0.9127368   Inf
#>     10 -0.1191186   Inf
'''
