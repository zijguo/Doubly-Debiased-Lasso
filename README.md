### Load the source functions
```R
source('Source_Functions.R', encoding = 'UTF-8')
```
### Example
##generate a random dataset with population n, dimension p, scarcity s, confounding variable count d

#uses default: n=300, p=500, s=5, q=3, sigmaE=2, sigma=2, pert=1


```R
dset_curr=generate_dataset()
```

## construct confidence intervals for beta1
```R
idx=1
result = CI(dset_curr$X,dset_curr$Y,idx)

```
