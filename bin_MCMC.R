#####################################################################################################################
### Perform MCMC to generate posterior samples of penalized spline coefficients and regularization parameters     ###
### This function is called by bin_Bayes_penalized_splines, not directly called by the user                       ###
### 
#####################################################################################################################



MCMC <- function(data,formula,ncolz,prior,multiplier1,multiplier2,multiplier3,multiplier4,burnin,B,thin, lambda, numIntKnots)
{
  
  library(MASS)
  library(Matrix)
  library(coda)
  library(glmnet)
  
  
  ####################################### MCMC initialization ##############################
  # determine the initial values for the glm regression coefficients (including intercept)
  
  y <- data$y
  x_matrix <- model.matrix(formula, data)[,-1]

  #wihtout regularization
  # initfit <- glm(formula=formula,data=data,family = "binomial")
  # theta <- initfit$coef[-1]
  # coef <- as.numeric(as.matrix(data[,-1]) %*% theta)
  # 
  # loglik_noreg <- sum(data$y*log(exp(coef)/(1+exp(coef))) +
  #     (1-data$y)*log(1/(1+exp(coef))))
  
  #with regularization 
  initfit <- glmnet(x_matrix, y, alpha = 0, family = "binomial", lambda = lambda, intercept = FALSE)
  theta <- coef(initfit)@x
  # coef <- as.numeric(as.matrix(data[,-1]) %*% theta)
  # loglik_reg <- sum(data$y*log(exp(coef)/(1+exp(coef))) +
  #     (1-data$y)*log(1/(1+exp(coef))))
  
  # if(!is.na(loglik_reg) & loglik_reg < loglik_noreg){
  #   
  #   #wihtout regularization
  #   initfit <- glm(formula=formula,data=data,family = "binomial")
  #   theta <- initfit$coef
  #   
  #   theta_g <- theta[1:(ncolz+2)]
  #   u_g <- theta_g[-c(1:2)]
  #   
  #   theta_h <- theta[(ncolz+3):length(theta)]
  #   u_h <- theta_h[-c(1:2)]
  #   
  #   sigma2_g <- 1/rgamma(1,shape=ncolz/2+prior$alpha,rate=sum(u_g^2)/2+prior$beta)
  #   sigma2_h <- 1/rgamma(1,shape=ncolz/2+prior$alpha,rate=sum(u_h^2)/2+prior$beta)
  #   
  #   data_g <- data[,2:(ncolz+3)]
  #   data_h <- data[,(ncolz+4):ncol(data)]
  #   
  #   temp_h <- as.numeric(as.matrix(data_h) %*% theta_h) #this will treat theta_h as t(theta_h)
  #   temp_g <- as.numeric(as.matrix(data_g) %*% theta_g)
  #   
  #   data_g$h <- temp_h
  #   data_h$g <- temp_g
  #   
  #   formula1 <- as.formula(paste("y ~ ", paste(colnames(data_g), collapse="+") ))
  #   inifit_g <- glm(formula=formula1, data=cbind(y=data[,1],data_g), family = "binomial")
  #   var_g <- (multiplier1^2)*vcov(inifit_g)[-c(ncol(data_g)),-c(ncol(data_g))]
  # 
  #   formula2 <- as.formula(paste("y ~ ", paste(c(colnames(data_h),"0"), collapse="+") ))
  #   inifit_h <- glm(formula=formula2, data=cbind(y=data[,1],data_h), family = "binomial")
  #   var_h <- (multiplier1^2)*vcov(inifit_h)[-c(ncol(data_h)),-c(ncol(data_h))]
  # 
  #   theta_initial <- theta
  #   
  #   
  # }else{
    
    # initfit <- glmnet(x_matrix, y, alpha = 0, family = "binomial", lambda = lambda)
    # theta <- coef(initfit)@x
    
    theta_g <- theta[1:(ncolz+1)]
    u_g <- theta_g[-c(1)]
    
    theta_h <- theta[(ncolz+2):length(theta)]
    u_h <- theta_h[-c(1:2)]
    
    sigma2_g <- 1/rgamma(1,shape=ncolz/2+prior$alpha,rate=sum(u_g^2)/2+prior$beta)
    sigma2_h <- 1/rgamma(1,shape=ncolz/2+prior$alpha,rate=sum(u_h^2)/2+prior$beta)
    
    data_g <- data[,2:(ncolz+2)]
    data_h <- data[,(ncolz+3):ncol(data)]
    
    temp_h <- as.numeric(as.matrix(data_h) %*% theta_h) #this will treat theta_h as t(theta_h)
    temp_g <- as.numeric(as.matrix(data_g) %*% theta_g)
    
    data_g$h <- temp_h
    data_h$g <- temp_g
    
    
    skip_to_next <- FALSE
    
    
    tryCatch({
      data_g_matrix <- as.matrix(data_g)
      inifit_g <- glmnet(data_g_matrix, y, alpha = 0, family = "binomial",
              lambda = lambda, intercept = FALSE) 
      linPred_g  <- as.numeric( data_g_matrix %*% inifit_g$beta)
      predProb_g <- exp(linPred_g) / (1+exp(linPred_g))
      predProb_g[which(is.na(predProb_g))] <- 1
      w_g <- Diagonal(nrow(data), x = predProb_g)
      var_z_g <- solve(w_g)
      p1_g <- solve(t(data_g_matrix)%*%w_g%*%data_g_matrix+lambda)
      p2_g <- t(data_g_matrix)%*%w_g%*%var_z_g%*%w_g%*%data_g_matrix
      var_g_0 <- p1_g%*%p2_g%*%p1_g
      var_g <- (multiplier1^2)*(var_g_0[-c(ncol(data_g)),-c(ncol(data_g))]) ##remove variance for h term
      
      
      data_h_matrix <- as.matrix(data_h)
      inifit_h <- glmnet(data_h_matrix, y, alpha = 0, family = "binomial", lambda = lambda, intercept = FALSE)
      linPred_h  <- as.numeric(data_h_matrix %*% inifit_h$beta)
      predProb_h <- exp(linPred_h) / (1+exp(linPred_h))
      predProb_h[which(is.na(predProb_h))] <- 1
      w_h <- Diagonal(nrow(data), x = predProb_h)
      var_z_h <- solve(w_h)
      p1_h <- solve(t(data_h_matrix)%*%w_h%*%data_h_matrix+lambda)
      p2_h <- t(data_h_matrix)%*%w_h%*%var_z_h%*%w_h%*%data_h_matrix
      var_h_0 <- p1_h%*%p2_h%*%p1_h
      var_h <- (multiplier1^2)*(var_h_0[-c(ncol(data_h)),-c(ncol(data_h))])
      
    
    }, error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { 
      
      #without regularization
      
      # formula1 <- as.formula(paste("y ~ ", paste(colnames(data_g), collapse="+") ))
      # inifit_g <- glm(formula=formula1, data=cbind(y=data[,1],data_g), family = "binomial")
      # var_g <- (multiplier1^2)*vcov(inifit_g)[-c(ncol(data_g)),-c(ncol(data_g))]
      # 
      # formula2 <- as.formula(paste("y ~ ", paste(c(colnames(data_h),"0"), collapse="+") ))
      # inifit_h <- glm(formula=formula2, data=cbind(y=data[,1],data_h), family = "binomial")
      # var_h <- (multiplier1^2)*vcov(inifit_h)[-c(ncol(data_h)),-c(ncol(data_h))]
      
      
      formula1 <- as.formula(paste("y ~ 0 +", paste(colnames(data_g), collapse="+") ))
      inifit_g_noreg <- glm(formula=formula1, data=cbind(y=data[,1],data_g), family = "binomial")
      var_g_0 <- vcov(inifit_g_noreg) - 2*lambda
      var_g <- (multiplier1^2)*var_g_0[-c(ncol(data_g)),-c(ncol(data_g))]
      
      
      formula2 <- as.formula(paste("y ~ ", paste(c(colnames(data_h),"0"), collapse="+") ))
      inifit_h_noreg <- glm(formula=formula2, data=cbind(y=data[,1],data_h), family = "binomial")
      var_h_0 <- vcov(inifit_h_noreg) - 2*lambda
      var_h <- (multiplier1^2)*var_h_0[-c(ncol(data_h)),-c(ncol(data_h))]
      
      
     } 
    
    # put together the initial values and proposal variance for all coefficients,
    theta_g <- c(as.numeric(inifit_g$beta)[-length(inifit_g$beta)])
    names(theta_g) <- NULL
    u_g <- theta_g[-c(1)]
    
    theta_h <- as.numeric(inifit_h$beta)[-length(inifit_h$beta)]
    u_h <- theta_h[-c(1:2)]
    
    theta <- c(theta_g,theta_h)
    
    theta_initial <- theta
    
    
  # }
  
  
  
  ########################################### MCMC starts #####################################

  if(is.na(thin) == TRUE) {
    L <- burnin+B
  }else{
    L <- burnin+B*thin
  }
  
  theta.all <- matrix(NA,nrow=L,ncol=length(theta))
  theta.res <- matrix(NA,nrow=B,ncol=length(theta))
  
  sigma2.all <- matrix(NA,nrow=L,ncol=2)
  sigma2.res <- matrix(NA,nrow=B,ncol=2)
  
  ll <- 0
  acpt_g <- 0
  acpt_h <- 0
  
  
  # indices corresponding to the control group (i.e., without interaction terms)
  idx_g <- 1:(ncolz+1)
  
  
  # calculate initial log likelihood (only the part relevant to the acpt prob calculation)
  data_g <- as.matrix(cbind(data[,2:(ncolz+2)]))
  data_h <- as.matrix(data[,(ncolz+3):ncol(data)])
  coef_g <- as.numeric(data_g %*% theta_g)
  coef_h <- as.numeric(data_h %*% theta_h)
  
  f = prior$f
  
  for(l in 1:L)
  {
    
    ##########  update theta_g  ##########
    newtheta_g <- mvrnorm(1,theta_g,var_g)
    newu_g <- newtheta_g[-c(1)]
    
    newtheta <- theta
    newtheta[idx_g] <- newtheta_g
    
    newcoef_g <- as.numeric(data_g %*% newtheta_g)
    
    newloglik <- sum(data$y*log(exp(newcoef_g+coef_h)/(1+exp(newcoef_g+coef_h))) +
        (1-data$y)*log(1/(1+exp(newcoef_g+coef_h))))
    # newloglik2 <- sum(log(1/(1+exp(newcoef_g+coef_h)))+
    #     data$y*(newcoef_g+coef_h-log(1)))      #another way to write out the log likelihood

    coef_g <- as.numeric(data_g %*% theta_g)
    loglik <- sum(data$y*log(exp(coef_g+coef_h)/(1+exp(coef_g+coef_h))) +
        (1-data$y)*log(1/(1+exp(coef_g+coef_h))))
    
    # loglik <- sum(log(1/(1+exp(coef_g+coef_h)))+
    #       data$y*(coef_g+coef_h-log(1))) 
    

    log_acpt <- newloglik-0.5*sum(newu_g^2)/sigma2_g - 0.5*sum(newtheta_g[1:2]^2)/f^2- 
      loglik+ 0.5*sum(u_g^2)/sigma2_g + 0.5*sum(theta_g[1:2]^2)/f^2
    
    
    if(!is.na(log_acpt)&log(runif(1))<log_acpt)
    {
      theta <- newtheta
      theta_g <- newtheta_g
      coef_g <- newcoef_g
      u_g <- newu_g
      acpt_g <- acpt_g+1
    }
    
    
    ##########  update sigma2_g  ##########
    sigma2_g <- 1/rgamma(1,shape=ncolz/2+prior$alpha,rate=sum(u_g^2)/2+prior$beta)
    
    
    ##########  update theta_h  ##########
    newtheta_h <- mvrnorm(1,theta_h,var_h)
    newu_h <- newtheta_h[-c(1:2)]
    
    newtheta <- theta
    newtheta[-idx_g] <- newtheta_h
    
    newcoef_h <- as.numeric(data_h %*% newtheta_h)
    newloglik <- sum(data$y*log(exp(coef_g+newcoef_h)/(1+exp(coef_g+newcoef_h))) +
        (1-data$y)*log(1/(1+exp(coef_g+newcoef_h))))
    
    coef_h <- as.numeric(data_h %*% theta_h)
    loglik <- sum(data$y*log(exp(coef_g+coef_h)/(1+exp(coef_g+coef_h))) +
        (1-data$y)*log(1/(1+exp(coef_g+coef_h))))
    
    log_acpt <- newloglik-0.5*sum(newu_h^2)/sigma2_h-0.5*sum(newtheta_h[1:2]^2)/f^2-
      loglik+0.5*sum(u_h^2)/sigma2_h + 0.5*sum(theta_h[1:2]^2)/f^2
    
    if(!is.na(log_acpt)&log(runif(1))<log_acpt)
    {
      theta <- newtheta
      theta_h <- newtheta_h
      coef_h <- newcoef_h
      u_h <- newu_h
      acpt_h <- acpt_h+1
    }
    
    
    ##########  update sigma2_h  ##########
    sigma2_h <- 1/rgamma(1,shape=ncolz/2+prior$alpha,rate=sum(u_h^2)/2+prior$beta)
    
    
    ##########  save the values of the new iteration  ##########
    theta.all[l,] <- theta
    sigma2.all[l,] <- c(sigma2_g,sigma2_h)
    
    if(is.na(thin) == TRUE){
      if( (l>burnin)){
        
        ll <- ll+1
        theta.res[ll,] <- theta
        sigma2.res[ll,] <- c(sigma2_g,sigma2_h)
        
      }
    }else{
      
      if( (l>burnin) && ((l-burnin)%%thin==0) )
      {
        ll <- ll+1
        theta.res[ll,] <- theta
        sigma2.res[ll,] <- c(sigma2_g,sigma2_h)
      }
    }
  
    
    ##########  adaptively adjust the proposal variance based on the acceptance probability of the current MCMC chain  ##########
    if( (l>2000) && (l%%2000==1) )
    {
      all_g = rbind(theta_initial[idx_g],theta.all[1:(l-1),idx_g])[(l-2000):l,]
      all_h = rbind(theta_initial[-idx_g],theta.all[1:(l-1),-idx_g])[(l-2000):l,]
      var_g=2.4^2*(cov(all_g)+1e-10*diag(rep(1,length(idx_g))))/length(idx_g)      #a block of 2000 samples
      var_h=2.4^2*(cov(all_h)+1e-10*diag(rep(1,length(idx_g)+1)))/(length(idx_g)+1)     #the numbers of terms in g and h are the same

      # var_g=2.4^2*(cov(theta.all[1:(l-1),idx_g])+1e-10*diag(rep(1,length(idx_g))))/length(idx_g)      #a block of 2000 samples
      # var_h=2.4^2*(cov(theta.all[1:(l-1),-idx_g])+1e-10*diag(rep(1,length(idx_g))))/length(idx_g)     #the numbers of terms in g and h are the same
      # 
      if(acpt_g/l < 0.1)
      {
        var_g <- (multiplier3^2)*var_g
      }
      
      else if(acpt_g/l < 0.2)
      {
        var_g <- (multiplier2^2)*var_g
      }
      
      else if (acpt_g/l > 0.23)
      {
        var_g <- (multiplier4^2)*var_g
      }

      if(acpt_h/l < 0.1)
      {
        var_h <- (multiplier3^2)*var_h
      }
      
      else if(acpt_h/l < 0.2)
      {
        var_h <- (multiplier2^2)*var_h
      }
      
      else if(acpt_h/l > 0.23)
      {
        var_h <- (multiplier4^2)*var_h
      }
    }
    
  }
  
  
  heidel.theta=rep(NA,ncol(theta.res))
  
  for(t in 1:ncol(theta.res))
  {
    attemp <- try(heidel.diag(theta.res[,t])[,3])
    if (inherits(attemp, "try-error")) return(list(success=FALSE))
    else{
      heidel.theta[t] = heidel.diag(theta.res[,t])[,3]
    }
  }
  
  heidel.sigma2=rep(NA,2)
  
  for(t in 1:2)
  {
    attemp <- try(heidel.diag(sigma2.res[,t])[,3])
    if (inherits(attemp, "try-error")) return(list(success=FALSE))
    else{
      heidel.sigma2[t] = heidel.diag(sigma2.res[,t])[,3]
    }
  }
  
  
  heidel.conv <- !( min(abs(heidel.theta))<abs(0.05/((numIntKnots+2+2)*2+2)) | 
      min(abs(heidel.sigma2))<abs(0.05/((numIntKnots+2+2)*2+2)))
  
  if(!is.na(heidel.conv) & heidel.conv)
  {
    return(list(success=TRUE, theta=theta.res, sigma2=sigma2.res,				# MCMC samples 
      acpt_g=acpt_g/L, acpt_h=acpt_h/L,							# acceptance rates
      heidel_theta=heidel.theta, heidel_sigma2=heidel.sigma2 ))			#  test statistics 
  }
  
  else
  {
    # return(list(success=FALSE)) 
    return(list(success=FALSE, theta=theta.res, sigma2=sigma2.res,				# MCMC samples 
      acpt_g=acpt_g/L, acpt_h=acpt_h/L,							# acceptance rates
      heidel_theta=heidel.theta, heidel_sigma2=heidel.sigma2 ))			#  test statistics 
  }
  
  
}




