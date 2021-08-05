library(ggplot2)
#define simulation parameters

#Number of simulations to run
N=300
#Index j to estimate Beta for
idx=1

coverage=0
width=0
est_dist=data.frame(Estimate="Estimate",Variance="Variance",B_beta="B_beta",B_b="B_b")

#This loop takes a significant amount of processing power/time to complete dependent on the size of n,p,and N
for (i in 1:N){
  #produces a random dataset with population n, dimension p, scarcity s, confounding variable count d,
  #independent variability in in Y sigmaE, independent variability in X sigma, and perturbation of H pert.
  #uses default: n=300, p=500, s=5, q=3, sigmaE=2, sigma=2, pert=1
  dset_curr=generate_dataset()

  result = CI(dset_curr$X,dset_curr$Y,idx,dset_curr$beta,dset_curr$b)

  if(idx<=3){
    if(result$lower<1 && result$upper>1)
    {
      coverage=coverage+1
    }
  }
  else{
    if(result$lower<0 && result$upper>0)
    {
      coverage=coverage+1
    }
  }
  dblasso=(result$upper+result$lower)/2
  width=width+(result$upper-result$lower)
  est_dist[i,]=c(dblasso,result$Variance)

  est_dist$B_beta[i]=result$B_beta
  est_dist$B_b[i]=result$B_b
  print(paste("Estimate for Simulation",i,": ",dblasso,"Current Coverage: ",coverage/i))
}

est_dist$Estimate=as.numeric(est_dist$Estimate)
est_dist$Variance=as.numeric(est_dist$Variance)
est_dist$B_b=as.numeric(est_dist$B_b)
est_dist$B_beta=as.numeric(est_dist$B_beta)

#plots N estimates for separate simulations and plots. The estimator sample standard deviation times z-alpha is
#approximately equal to the half width of the constructed confidence interval.
#The estimator also closely follows ~N(beta_true,Estimator_Variance) which indicates that double debiasing has been successful.
ggplot(est_dist, aes(x=Estimate)) + geom_histogram(binwidth=1/(8*sqrt(n)))
print(paste("Average Estimate: ",mean(est_dist$Estimate)))
print(paste("Coverage",coverage/N))
print(paste("Average Width: ",width/N))

#external finds the variance of N estimates, internal finds the variance for each simulation and averages
print(paste("Externally Calculated Estimator Variance: ",var(est_dist$Estimate)))
print(paste("Average Internally Calculated Estimator Variance: ",mean(est_dist$Variance)))

print(paste("Average B_beta: ",mean(est_dist$B_beta),"Average Absolute B_beta: ",mean(abs(est_dist$B_beta))))
print(paste("Average B_b: ",mean(est_dist$B_b),"Average Absolute B_b: ",mean(abs(est_dist$B_b))))

