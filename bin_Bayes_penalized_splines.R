############################################################################################################
### Implement Bayesian penalized splines to estimate the response rate for binary endpoint Y,            ###
### as a function of a continuous marker x.	                                                             ###
### Two arms (one control group and one experimental group) are assumed. 	                               ###
### Parameters:                                                                                          ###
### data = a data frame with the observations arranged by row,                                           ###
### and including the columns below (in the following order):                                            ###
###        1) y: the response outcome (=1 has response; = 0 no response) (int)                           ###
###        2) x: the value of the continuous marker (double)                                             ###
###        3) trt: the group indicator (=1 for experimental group; =0 for control group) (int)           ###
### next_x =  the marker value for the about to be assigned patients (double)                            ###
### xrange = a vector of length 2 which gives the range of marker values (vector)                        ###
### numIntKnots = the number of interior knots (should scale with number of response)                    ###
### MCMC.specs = a list of MCMC specifications including                                                 ###
###        	1) burnin: the number of burnin samples (int)                                                ###
###        	2) B: the number of posterior samples (int)                                                  ###
###        	3) thin: the thinning parameter (double)                                                     ###
###        	4) multiplier: a scaling factor for the size of random walk proposal (double)                ###
### Users might need to tune MCMC.specs$thin and MCMC.specs$multiplier to                                ###
### make sure that the acceptance probabilities and trace plots are reasonable                           ###
### lambda = shrinkage on paratmer estimation at initialization                                          ###
### ref = reference marker value in the control group (double)                                           ###
### alpha_e = the lower significance level for the credible interval                                     ###
### alpha_h = the upper significance level (1-alpha_h) for the credible interval                         ###
############################################################################################################



Bayes_penalized_splines <- function(data, 				
  next_x = NA,              
  xrange, 							   
  numIntKnots=9,						
  MCMC.specs=list(burnin=10000,B=2000,thin=10,multiplier=0.5), 
  lambda = 0.001,
  ref, 								    
  alpha_e=0.05,             
  alpha_h = 0.05,
  ){						 

  ### Code here to load in the functions that will be called later
  source("bin_Get_splines.R")
  source("bin_MCMC.R")
  source("bin_jointband.R")
  
  
  # Compute the O'Sullivan splines design matrix for the continuous marker x, and the interaction term for x and trt
  res_splines <- Get_splines(data,xrange,numIntKnots,next_x = next_x)
  
  
  # A numerical stability check for the spline design matrix
  # If not passing the stability check, stop the function call now
  if(res_splines$stabCheck)
  {
    return(list(x=sort(data$x), success=FALSE))
  }
  
  
  # If passing the stability check, perform MCMC
  else
  {
    data <- res_splines$data
    formula <- res_splines$formula
    ncolz <- res_splines$ncolz
    
    for(count in 1:25)
    {			
      res_MCMC <- MCMC(data, formula, ncolz, prior=list(alpha=0.01, beta=0.01, f=10000),
        multiplier1=MCMC.specs$multiplier, multiplier2=0.9, multiplier3=0.6,	multiplier4=1.1,	# scaling factors to adaptively adjust the size of random walk proposal
        burnin=MCMC.specs$burnin, B=MCMC.specs$B, thin=MCMC.specs$thin, lambda = lambda, numIntKnots)
      
      print(paste0("count=",count))
      
      if(res_MCMC$success) { 
        break
      } else if(count==25)
      {
        return(list(x=sort(data$x), success=FALSE))
      }
    }
    
    
    
    
    # Using the posterior samples of the spline coefficients,
    # get the posterior samples of the logit probability of response in the control group at the reference value,
    # as a function of marker value for the control group.  [g(x)- g(x0)]
    # This approach results in pinching in the credible band at the reference value. 
    data_g <- data[order(data$x),2:(ncolz+2)]                     #intercept term is just 1 and can be canceled out
    idx <- which.min(abs(data_g$x-ref))
    data_g <- sweep(data_g,2,as.numeric(data_g[idx,]))            #each row subtracts out the row with marker value closest to the reference level
    
    posterior_g <- res_MCMC$theta[,1:(ncolz+1)] %*% t(data_g)     
    simCR <- jointband(posterior_g,alpha_e)
    upr_g <- simCR$upr_CI
    lwr_g <- simCR$lwr_CI
    upr_g_pt <- apply(posterior_g, 2, function(x){quantile(x,1-alpha_h/2)})
    lwr_g_pt <- apply(posterior_g, 2, function(x){quantile(x,alpha_e/2)})
    
    
    # Using the posterior samples of the spline coefficients,
    # get the posterior samples of the difference in the logit probability of response between the two arms, as a function of marker value.
    # all patients on experiment arm vs all patients on control arm  [g(x) + h(x)*1 - g(x)]
    # Note that the group difference does not depend on the choice of the reference value
    data_h <- cbind(rep(1,nrow(data)),data[order(data$x),2:(ncolz+2)])

    
    posterior_h <- res_MCMC$theta[,-c(1:(ncolz+1))] %*% t(data_h)
    simCR <- jointband(posterior_h,alpha_e)
    upr_h <- simCR$upr_CI
    lwr_h <- simCR$lwr_CI
    upr_h_pt <- apply(posterior_h, 2, function(x){quantile(x,1-alpha_h/2)})
    lwr_h_pt <- apply(posterior_h, 2, function(x){quantile(x,alpha_e/2)})
    
    
    # Using the posterior samples of the spline coefficients,
    # get the posterior samples of the logit probability of response to the control group at the reference value,
    # as a function of marker value for the experimental group.
    posterior_exp <- posterior_g + posterior_h
    simCR <- jointband(posterior_exp,alpha_e)
    upr_exp <- simCR$upr_CI
    lwr_exp <- simCR$lwr_CI    
    upr_exp_pt <- apply(posterior_exp, 2, function(x){quantile(x,1-alpha_h/2)})
    lwr_exp_pt <- apply(posterior_exp, 2, function(x){quantile(x,alpha_e/2)})
    
    
    if(!is.na(next_x)){
      
      next_pt <- res_splines$next_pt
      next_pt_g <- next_pt[,3:(ncolz+3)]                     
  
      next_pt_g_byref <- sweep(next_pt_g,2,as.numeric(data_g[idx,]))
      
      posterior_g_next <- res_MCMC$theta[,2:(ncolz+2)] %*% t(next_pt_g_byref)     
      
      next_pt_h <- next_pt[,2:(ncolz+3)]
      
      posterior_h_next <- res_MCMC$theta[,-c(1:(ncolz+2))] %*% t(next_pt_h)
      
      posterior_exp_next <- posterior_g_next + posterior_h_next
    
      #Exponentiate all estimates
      expo_posterior_g_next <- exp(posterior_g_next)
   

      expo_posterior_exp_next <- exp(posterior_exp_next)

      
      expo_post_next_pt <- list(expo_post_g_next = as.vector(t(expo_posterior_g_next)),
        expo_post_exp_next = as.vector(t(expo_posterior_exp_next))
      )
      
    }else{
      
      expo_post_next_pt <- NA
      
    }
    # return a list with the MCMC samples, posterior mean estimate and credible band for the log hazard ratio as a function of marker value
    return(list(x=sort(data$x), success=res_MCMC$success,										# sorted marker values x and indicator of success for the function call
      theta=res_MCMC$theta, sigma2=res_MCMC$sigma2,								# MCMC samples of the model parameters in Bayesian penalized splines
      est_ctl=apply(posterior_g,2,mean), upr_ctl=upr_g, lwr_ctl=lwr_g,					# estimate and joint band for the LHR of the control group
      upr_ctl_pt=upr_g_pt, lwr_ctl_pt=lwr_g_pt,								# pointwise band for the LHR of the control group
      est_diff=apply(posterior_h,2,mean), upr_diff=upr_h, lwr_diff=lwr_h,				# estimate and joint band for the difference in LHR between two arms
      upr_diff_pt=upr_h_pt, lwr_diff_pt=lwr_h_pt,								# pointwise band for the LHR difference between two arms
      est_exp=apply(posterior_exp,2,mean), upr_exp=upr_exp, lwr_exp=lwr_exp,				# estimate and joint band for the LHR of the experimental group
      upr_exp_pt=upr_exp_pt, lwr_exp_pt=lwr_exp_pt,								# pointwise band for the LHR of the experimental group
      MCMC_ctl=posterior_g, MCMC_diff=posterior_h, MCMC_exp=posterior_exp,				# MCMC samples of the LHR as a function of marker
      acpt_ctl=res_MCMC$acpt_g, acpt_diff=res_MCMC$acpt_h, 						# acceptance rates for two MCMC chains
      heidel_theta=res_MCMC$heidel_theta, heidel_sigma2=res_MCMC$heidel_sigma2,
      expo_post_next_pt = expo_post_next_pt,
      previous_ref = as.numeric(data_g[idx,]),
      B_matrix = res_splines$B_matrix, LZ_matrix = res_splines$LZ_matrix,
      g_est = res_MCMC$theta[,1:(ncolz+1)],
      h_est = res_MCMC$theta[,-c(1:(ncolz+1))]

      ))			
  }
}



