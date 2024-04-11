##################################################################################################################
### Implement adaptive randomization based on estimate the response rate for previously observed endpoint Y    ###
### as a function of continuous marker x and the next patient's marker profile.	                               ###
### Two arms (one control group and one experimental group) are assumed.                                   	   ###
### Parameters:                                                                                                ###
### x.eval = a vector of about to be evaluated continuous marker (vector)                                      ###
### previous.n = the number of enrolled patients (int)                                                         ###
### n = total sample size planned for the study (int)                                                          ###
### tuning = the hyper-parameter used to calculate probability of one arm is superior than the other (double)  ###
### random_e = lowest allowed randomized probability to experiment,                                            ###
###           and highest probability = 1- random_e (double)                                                   ###
### B_matrix = basis matrix from previous enrolled patients information (matrix)                               ###
### LZ_matrix = linear transformed basis matrix                                                                ###
### previous_ref = reference marker value used previously (double)                                             ###
### g_est = estimated g_matrix (prognostics association) previously (matrix)                                   ###
### h_est = estimated h_matrix (predictive association) previously (matrix)                                    ###
###                                                                                                            ###
### A list with the first object as the output (vector),                                                       ###
### which is the assigned treatment arms for each incoming patients                                            ### 
##################################################################################################################

adaptive_randomization <- function(
  x.eval,
  tuning = 0.5,          
  previous.n,
  n ,
  random_e = 0.1,         
  B_matrix,
  LZ_matrix,
  previous_ref,
  g_est,
  h_est

  ){
  
    new_B <- predict(B_matrix,x.eval)
    
    new_Z <- new_B%*%LZ_matrix
    
    new_Z <- scale(new_Z)
    
    next_pt <- as.data.frame(cbind(c(1), x.eval, new_Z))
    
    next_pt_g <- next_pt[,-c(1)]
    
    next_pt_g_byref <- sweep(next_pt_g,2,previous_ref)
    
    posterior_g_next <- g_est %*% t(next_pt_g_byref) 
    
    next_pt_h <- next_pt
    
    posterior_h_next <- h_est %*% t(next_pt_h)
    
    posterior_exp_next <- posterior_g_next + posterior_h_next
    
    expo_post_g_next <- exp(posterior_g_next)
    
    expo_post_exp_next <- exp(posterior_exp_next)
    
    next_pt_prob_ctl <- expo_post_g_next/(1 + expo_post_g_next)
    
    next_pt_prob_exp <- expo_post_exp_next/(1 + expo_post_exp_next)
    
    #NA value will be filled with 1
    next_pt_prob_ctl[which(is.na(next_pt_prob_ctl))] <- 1
    
    next_pt_prob_exp[which(is.na(next_pt_prob_exp))] <- 1
    
    
    prob_exp_superior <- apply(next_pt_prob_exp >= next_pt_prob_ctl,2,sum)/nrow(next_pt_prob_exp)
    
    if(tuning == 0.5){
      randomization_ratio <- prob_exp_superior^tuning/(prob_exp_superior^tuning + (1- prob_exp_superior)^tuning)
    
      }else{

      tuning <- (previous.n+c(1:length(x.eval)))/(2*n)
      
      randomization_ratio <- prob_exp_superior^tuning/(prob_exp_superior^tuning + (1- prob_exp_superior)^tuning)

    }
    
    randomization_ratio[randomization_ratio < random_e] <- random_e
    
    randomization_ratio[randomization_ratio > 1-random_e] <- 1-random_e
    
    
    output <- rbinom(length(x.eval),1,randomization_ratio)
    
    return(list(output = output))

}







