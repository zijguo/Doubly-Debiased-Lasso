### Load the source functions
```R
source('Source_Functions.R', encoding = 'UTF-8')
```

### generate data
We offer two examples to illustrate how to generate the dataset with population n, dimension p, scarcity s, confounding variable count q
#### Low Dimension Setting (n=2000, p=200, s=5, q=3, sigmaE=2, sigma=2, pert=1)
```R
n=2000
p=200
s=5
q=3
sigmaE=2
sigma=2
pert=1

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
#### High Dimension Setting (n=500, p=2000, s=5, q=3, sigmaE=2, sigma=2, pert=1)
```R
n=500
p=2000
s=5
q=3
sigmaE=2
sigma=2
pert=1

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
Or you can use the function (uses default: n=300, p=500, s=5, q=3, sigmaE=2, sigma=2, pert=1)
```R
dset_curr=generate_dataset()
```


### construct confidence intervals for beta1
```R
idx=1
# rho=0.5,rhop=0.5(default)
result = CI(dset_curr$X,dset_curr$Y,idx)

```
