# Spectral-Deconfounding
## This R code is designed based on the paper titled: "Doubly Debiased Lasso: High-Dimensional Inference under Hidden Confounding" 
## This work was written and developed by Zijian Guo, Domagoj Cevid, and Peter Buhlmann
DDLConfidence Interval has several functions designed to create a simulation dataset with 
an independent and dependent variable, each with their own errors, along with hidden confounding
Using a newly developed spectral decomposition, the previously created LASSO regression is 
implemented and modified to remove multiple forms of bias from the estimator. 
This package allows the Doubly Debiased LASSO algorithm to be applied to a generated dataset 
to estimate Beta-j and provide a confidence interval for the estimator. 
Additionally, the bias of the estimator can be calculated if the true values of beta and b are known. 
