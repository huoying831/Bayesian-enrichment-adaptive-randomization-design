###################################################################################################
### Calculate (1-alpha)*100% simultaneous band based on the matrix of posterior samples MCMC_P  ###
### This function is called by bin_Bayes_penalized_splines, not directly called by the user     ###
### Parameters:                                                                                 ###
### MCMC_P: calculated posterior samples (matrix)                                               ###
### alpha: significance level                                                                   ###
###################################################################################################


jointband <- function(MCMC_P,alpha)
{
  B <- dim(MCMC_P)[1]
  
  mean_P <- apply(MCMC_P,2,mean)
  sd_P <- apply(MCMC_P,2,sd)
  
  z_P <- rep(NA,B)
  
  keep <- (!is.na(sd_P)) & (sd_P!=0) 
  
  for(j in 1:B)
  {
    z_P[j] <- max(abs((MCMC_P[j,keep]-mean_P[keep])/sd_P[keep]))
  }
  
  upr_CI <- mean_P + quantile(z_P,1-alpha)*sd_P
  lwr_CI <- mean_P - quantile(z_P,1-alpha)*sd_P
  
  return(list(upr_CI=upr_CI,lwr_CI=lwr_CI))
}


