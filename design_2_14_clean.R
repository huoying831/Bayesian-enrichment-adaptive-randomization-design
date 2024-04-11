adaptive_enrich_design <- function(
    n,
  seed,
  MCMC.specs = list(burnin=20000,B=2000,thin=20,multiplier=0.5), 
  MCMC.specs.fut = list(burnin=20000,B=2000,thin=20,multiplier=0.5),
  marker.specs = list(dist="beta",param=c(1,1),range=c(0,1),ref="median"), 	                        
  percent.exp = 0.5,   
  FUN.ctl = function(x) {0*x-2.95},              # mapping function of the log hazard as a function of marker for the control group
  FUN.exp = function(x) {0*x-2.95},              # mapping function of the log hazard as a function of marker for the experimental group
  maxIntKnots = 6,
  adaptive_r = FALSE,   
  start_num = 300, 
  block_size = 50, 
  tuning = 0.5,  
  alpha = 0.05, 
  default_marker_seq = FALSE,
  lambda = 0.001, 
  p_eff_all = rep(0.95, 2),  			  # a vector of overall cohort posterior thresholds to claim efficacy (length=number of checks)
  p_eff_grp = rep(0.95, 2),  				# a vector of marker subgroup posterior thresholds to claim efficacy (length=number of checks)
  p_fut_all = rep(0.05, 1),				# a vector of overall cohort posterior thresholds for early stopping due to futility (length=number of checks - 1)
  p_fut_grp = rep(0.05, 1),				# a vector of marker subgroup posterior thresholds for early stopping due to futility (length=number of checks - 1)
  stop_eff = TRUE,
  random_e = 0.1,
  futility_loop = 100,
  p_fut = 0.1){
  
  # n = 500
  # seed = 3
  # MCMC.specs = list(burnin=2000,B=1000,thin=5,multiplier=0.5)
  # marker.specs = list(dist="beta",param=c(1,1),range=c(0,1),ref="median")
  # percent.exp = 0.5
  # FUN.ctl = function(x) {0*x-2.95}              # mapping function of the log hazard as a function of marker for the control group
  # FUN.exp = function(x) {0*x-2.95}             # mapping function of the log hazard as a function of marker for the experimental group
  # maxIntKnots = 6
  # adaptive_r = TRUE
  # start_num = 300
  # block_size = 100
  # tuning = 1
  # alpha = 0.05
  # default_marker_seq = FALSE
  # lambda = 0.001
  #  stop_eff = TRUE
  # random_e = 0.1
  # futility_loop = 10
  # p_fut = 0.1
  
  library(dplyr)  
  #source("bin_sim_all_data.R")
  source("bin_Bayes_penalized_splines.R")   
  source("adaptive_random.R")
  set.seed(seed)
  
  # #number of interim analysis
  # K <- floor((n-start_num)/block_size) + 1
  check_num <- seq(start_num,n,by = block_size)
  # check_num = ceiling(check*n)
  K = length(check_num)

  
  ##### Generate Patient Data #####
  ##### marker from 0.01 to 0.99 a fix sequence, sample to control/experimental arms
  
  sim_data <- function(n = 300, 								                                       # total trial size (both arms combined)
    percent.exp = 0.5,                                             # randomization ratio to the experimental group
    #marker.specs = list(dist="beta", param=c(1,1),range=c(0,1), ref="median"),        	# a list giving the distribution, range, and reference point of the continuous marker x                                     
    FUN.ctl = function(x) {0*x-2.95},                                   # mapping function of the probability of response as a function of marker x for the control group g(x)
    FUN.exp = function(x) {0*x-2.95} )                                  # mapping function of the probability of response as a function of marker x for the experimental group g(x)+h(x)
    {
    # Compute number of patients in each arm
    nsize.exp <- round(n*percent.exp)
    nsize.ctl <- n-nsize.exp
    
    marker_seq <- round(seq(0.01, 0.99, length.out = 99),6)
    
    #make sure reference value is in x.ctl
    ref <- median(marker_seq)
    #ref <- marker_seq[which.min(abs(marker_seq-ref))]
    x.ctl <- c(ref, sample(marker_seq,size = nsize.ctl-1, replace = TRUE))
    
    #make sure all marker values are covered

    marker_need <- marker_seq[!(marker_seq %in% x.ctl)]
    if(length(marker_need) == 0 ){
      
      x.exp <- sample(marker_seq,size = nsize.exp, replace = TRUE)
      
    }else{
      x.exp <- c(marker_need,sample(marker_seq,size = nsize.exp-length(marker_need), replace = TRUE))
    }
    
    
    

    
    # Calculate the probability of response, as a function of continuous marker x, relative to the reference point in the control group
    x.ctl_eff.ctl <- unlist(lapply(x.ctl,FUN.ctl))-FUN.ctl(ref)
    x.ctl_eff.exp <- unlist(lapply(x.ctl,FUN.exp))-FUN.ctl(ref)
    
    
    x.exp_eff.exp <- unlist(lapply(x.exp,FUN.exp))-FUN.ctl(ref)
    x.exp_eff.ctl <- unlist(lapply(x.exp,FUN.ctl))-FUN.ctl(ref)
    
    # Generate endpoint y based on the function of marker value
    x.clt_y.ctl <- rbinom(nsize.ctl,1,exp(x.ctl_eff.ctl)/(1+exp(x.ctl_eff.ctl)))	     # the control group
    x.exp_y.exp <- rbinom(nsize.exp,1,exp(x.exp_eff.exp)/(1+exp(x.exp_eff.exp)))	     # the experimental group
    
    x.clt_y.exp <- rbinom(nsize.ctl,1,exp(x.ctl_eff.exp)/(1+exp(x.ctl_eff.exp)))	     # the control group
    x.exp_y.ctl <- rbinom(nsize.exp,1,exp(x.exp_eff.ctl)/(1+exp(x.exp_eff.ctl)))	     # the experimental group
    
    
    # Set up treatment arm indicators
    trt <- c(rep(0, nsize.ctl), rep(1, nsize.exp),  rep(1, nsize.exp),rep(0, nsize.ctl))					# Treatment assignment 
    
    #label marker index(for merging purpose if marker values have duplication)
    #no duplication when initialization
    marker.order <- c(1:n,1:n)
    
    all_data <- data.frame(y=c(x.clt_y.ctl, x.exp_y.exp,x.clt_y.exp,x.exp_y.ctl), 
      x=c(x.ctl,x.exp,x.ctl,x.exp), trt= trt, marker.order = marker.order)
    
    return(list(sim_all_data = all_data, x = all_data$x[c(1:n)], ref = ref))
  }
  
  #elongate previous all_res_dat after change eligiblity creteria (mainly for generate response/marker.order)
  sim_data_rem <- function(
    x.rem,
    n.rem,
    ref,
    FUN.ctl = function(x) {0*x-2.95},                                   # mapping function of the probability of response as a function of marker x for the control group g(x)
    FUN.exp = function(x) {0*x-2.95} ){
    
    x.rem <- round(x.rem, 6)
    
    # Calculate the probability of response, as a function of continuous marker x, relative to the reference point in the control group
    x.rem_eff.ctl <- unlist(lapply(x.rem,FUN.ctl))-FUN.ctl(ref)
    x.rem_eff.exp <- unlist(lapply(x.rem,FUN.exp))-FUN.ctl(ref)
   
    
    # Generate endpoint y based on the function of marker value
    x.rem_y.ctl <- rbinom(n.rem,1,exp(x.rem_eff.ctl)/(1+exp(x.rem_eff.ctl)))	     # the control group
    x.rem_y.exp <- rbinom(n.rem,1,exp(x.rem_eff.exp)/(1+exp(x.rem_eff.exp)))	     # the experimental group
    
    trt <- c(rep(0, n.rem), rep(1, n.rem))	
    
    marker.order <- c(1:n.rem,1:n.rem)
    
    rem_data <- data.frame(y=c(x.rem_y.ctl, x.rem_y.exp), 
      x=c(x.rem,x.rem), trt= trt, marker.order = marker.order)
    
    return(rem_data)
  }
  
  futility_predict_response <- function(
    x.eval,
    trt,
    B_matrix,
    LZ_matrix,
    previous_ref,
    g_est,
    h_est
  ){
    
    new_B <- predict(B_matrix,x.eval)
    
    new_Z <- new_B%*%LZ_matrix
    
    next_pt <- as.data.frame(cbind(c(1), x.eval, new_Z))
    
    next_pt_g <- next_pt[,-c(1)]
    
    next_pt_g_byref <- sweep(next_pt_g,2,previous_ref)
    
    posterior_g_next <- g_est %*% t(next_pt_g_byref) 
    
    next_pt_h <- next_pt
    
    posterior_h_next <- h_est %*% t(next_pt_h)

    posterior_next <- posterior_g_next + trt * posterior_h_next 
    
    expo_post <- exp(posterior_next)
    
    next_pt_prob <- expo_post/(1 + expo_post)
    
    next_pt_prob[which(is.na(next_pt_prob))] <- 1
    
    reponse <- rbinom(length(x.eval),1,next_pt_prob)
    
    return(reponse)
  }
  
  
  
  #use sim_data
  sim_data <- sim_data(n=n, percent.exp=percent.exp, FUN.ctl=FUN.ctl, FUN.exp=FUN.exp)
  all_res_data <- sim_data$sim_all_data
  all_marker <- sim_data$x
  ref <- sim_data$ref  
  
  # temp <- all_marker[1]
  # switch_pos <- which(abs(all_marker-ref) == min(abs(all_marker - ref)))
  # all_marker[switch_pos] <- temp
  # all_marker[1] <- ref
  
  #for adpative == FALSE, sample for patients before and after interim
  # sample_stage <- sample(1:n, size = check_num[1])
  # if(! c(1) %in% sample_stage){
  #   sample_stage[1] <- 1
  # }
  
  # res_data <- data.frame(y = NA, x = NA, trt = NA)
  
  current_rr = percent.exp
  
  if(adaptive_r == TRUE){
    ##### generate dataset before the first interim
    
    #make sure equal randomization
    half_start <- ceiling(start_num/2)
    pre_x <- all_marker[1:(half_start*2)]
    pre_trt <- c(rep(0,half_start), rep(1,half_start))
    pre_marker_order <- 1:start_num
    
    # pre_marker_order[1] <- switch_pos
    # 
    # if(switch_pos <= start_num){
    #   pre_marker_order[switch_pos] <- 1
    # }
    
    
    pre_start <- data.frame(x = pre_x, trt = pre_trt, marker.order = pre_marker_order)
    
    res_data <- merge(pre_start, all_res_data, by = c("x", "marker.order", "trt"), all.x = TRUE, all.y = FALSE) %>% 
      dplyr::select(y, x, trt,marker.order) %>% arrange(marker.order)
  
    B_matrix = NA
    LZ_matrix = NA
    previous_ref = NA
    g_est = NA
    h_est = NA
   
    marker.stop = NULL     #stopped marker values
    marker.continuing = all_marker
    
    marker.eff = NA
    marker.noeff = NA
    marker.inf = NA
    marker.noinf = NA
    marker.fut = NA
    
    marker.enrolled = res_data$x
    marker.enrolled.status = NULL

    previous.j = 1 #index for check_num, to record number of newly enrolled pats
    
    ##### start interim + final analysis
    for(j in 1:K){
      
      if(nrow(res_data) %in% check_num[j]){
        
        ## for continuous all
         print(paste0("interimcheck",nrow(res_data)))      
        result.j <- Bayes_penalized_splines(data = res_data, xrange = marker.specs$range, numIntKnots=pmin(round(nrow(res_data)/20), maxIntKnots),    
          MCMC.specs=MCMC.specs, ref=ref, alpha_e=alpha, alpha_h = alpha, lambda = lambda,default_marker_seq = default_marker_seq)
        
        if(!result.j$success) {
          return(list(fail.MCMC=TRUE))
        } else{
          
          B_matrix = result.j$B_matrix
          
          LZ_matrix = result.j$LZ_matrix
          
          previous_ref = result.j$previous_ref
          
          g_est = result.j$g_est
          
          h_est = result.j$h_est
          
          ### save splines fit result#####
          fit.diff.all <- data.frame(x=result.j$x, est=result.j$est_diff, upr=result.j$upr_diff, lwr=result.j$lwr_diff,
            upr_pt=result.j$upr_diff_pt, lwr_pt=result.j$lwr_diff_pt)

          
          # check efficacy, inferiority, futility for j !=K
          if(j != K){
            
            #meaning no previous restriction on enrollment
            if(length(marker.stop) == 0){

              #marker value under evaluation/new enrolled  (enrolled = evaluated with no previous restriction)
              marker.j <- res_data$x

              fit.diff <- fit.diff.all %>% distinct(x, .keep_all = TRUE)
             
              #marker value not enrolled till the end of the trial
              marker.rest <- all_marker[(check_num[j]+1):length(all_marker)]
              
              # enrolled previous block + new enrolled + will be enrolled in future
              marker.all <- c(res_data$x,marker.rest)
              marker.all.status <- rep(0, length(marker.all))
              
            }else{
              
              #need to update marker.enroll.status using the new fitted model
              # what if new label is different from previous label (from -1 to 2 etc.)
              # currently, only change status if 0

              marker.j <- res_data$x[(check_num[previous.j]+1):check_num[j]]
              # marker.j.status <- rep(0, length(marker.j))  
              # # fit.diff <- fit.diff.all %>% filter(x %in% marker.j) %>% unique()
              fit.diff <- fit.diff.all %>% distinct(x, .keep_all = TRUE)
              # 
              # #all remaining marker value should be labelled as 0
              marker.rest <- all_marker[(check_num[j]+1):length(all_marker)]
              # marker.rest.status <- rep(0, length(marker.rest))
              # 
              # #no need to change marker.enroll/marker.enrolled.status
              # 
              
              marker.all <- c(marker.enrolled,marker.j,marker.rest)
              marker.all.status <- c(marker.enrolled.status,rep(0, length(marker.j)), rep(0, length(marker.rest)))
            }
            
            #if efficacy, labelled as 2
            #if inferiority, labelled as -1
            #if futility, labelled as 1
            #if continuing, labelled as 0
            
            #check for efficacy
            if(any(fit.diff$lwr_pt>0)) {
             
              marker.all.status_temp <- marker.all.status
              
              marker_eff <- fit.diff[fit.diff$lwr_pt>0,]$x
              marker_noeff <-fit.diff[fit.diff$lwr_pt<=0,]$x
              
              #if efficacy, labelled as 2
              
              if(length(marker_noeff) == 0){
                
                marker.all.status <- ifelse(marker.all.status == 0, 2, marker.all.status)
                
                marker.all.status <- marker.all.status + marker.all.status_temp
                
              }
              else{
                
                for(k in 1:length(marker.all)) {
                  
                  marker.all.status[k] <- ifelse(min(abs(marker.all[k]-marker_eff)) <= min(abs(marker.all[k]-marker_noeff)) & 
                      marker.all.status[k] == 0, 2,0)
                }
                
                #overlaid efficacy indicator and previous indicator
                marker.all.status <- marker.all.status + marker.all.status_temp
                
              }
              
            }
            
            #check for inferiority
            if(any(fit.diff$upr_pt<0)){
            
            
            marker_inf <- fit.diff[fit.diff$upr_pt<0,]$x
            marker_noinf <-fit.diff[(fit.diff$upr_pt>=0),]$x #including marker_eff
            
           
            #save marker.j.status/marker.rest.status from with efficacy indicator
            marker.all.status_temp <- marker.all.status
            
            if(length(marker_noinf) == 0){
              
              
              marker.all.status <- ifelse(marker.all.status == 0, -1, marker.all.status)
              
              marker.all.status <- marker.all.status + marker.all.status_temp
              
            }else{
              
              for(k in 1:length(marker.all)) {
                marker.all.status[k] <- ifelse(min(abs(marker.all[k]-marker_inf)) <= min(abs(marker.all[k]-marker_noinf)) &
                    marker.all.status[k] == 0 , -1, 0)
              }
              
              #overlaid efficacy indicator and inferiority indicator
              marker.all.status <- marker.all.status + marker.all.status_temp
              
              
            }
          }
            
            
            #check futility for marker status = 0
            
            marker_check_fut <- marker.all[marker.all.status == 0]
            
            if(length(marker_check_fut) > 0){
              
              fut.matrix <- matrix(NA, nrow = futility_loop, ncol = n) 
              
              print(paste0("futilitycheck",nrow(res_data))) 
              
              
              for(l in 1:futility_loop){
                
                result.fut.success <- FALSE
                while_loop_counter <- 1
                
                while(result.fut.success != TRUE & while_loop_counter < 5 ){
                  
                  tryCatch({
                    
                    print(paste0("while_loop_",while_loop_counter)) 
                    
                    n.rem <- n-length(res_data$x)
                    x.rem <- sample(x = marker_check_fut, size =n.rem, replace = TRUE)
                    
                    
                    trt.rem <- adaptive_randomization(x.eval = x.rem, n = n, previous.n = length(res_data$x),
                      tuning = tuning, 
                      random_e =random_e,
                      B_matrix = B_matrix, LZ_matrix = LZ_matrix, previous_ref = previous_ref, 
                      g_est= g_est, h_est = h_est)
                    
                    reponse.rem <-futility_predict_response( x.eval = x.rem,
                      trt = trt.rem$output,
                      B_matrix = B_matrix,
                      LZ_matrix = LZ_matrix,
                      previous_ref = previous_ref,
                      g_est = g_est[sample(1:nrow(g_est),1),],
                      h_est = h_est[sample(1:nrow(h_est),1),])
                    
                    dat_rem <- data.frame(y= reponse.rem, x = x.rem, trt = trt.rem$output)
                    
                    futility_dat <- rbind(res_data[,c("y","x","trt")], dat_rem)
                    
                    result.fut <- Bayes_penalized_splines(data = futility_dat, xrange = marker.specs$range, numIntKnots=pmin(round(nrow(res_data)/20), maxIntKnots),    
                      MCMC.specs=MCMC.specs.fut, ref=ref, alpha_e=alpha, alpha_h = alpha, lambda = lambda,default_marker_seq = default_marker_seq)
                    
                    result.fut.success <- result.fut$success
                    
                    while_loop_counter <- while_loop_counter+1
                    
                  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
                  
                 
                }
                
                
                if(!result.fut.success) {
                  return(list(fail.MCMC=TRUE))
                } else{
                  
                  fut.matrix[l,] <- as.numeric(result.fut$lwr_diff_pt > 0)
                  
                }
                
              }
              
              mean.fut <- apply(fut.matrix,2,mean)
              fut.res <- data.frame(x = result.fut$x, if_fut = as.numeric(mean.fut < p_fut)) #if_fut = 1 for marker value categorized as futility
              
              marker.fut <- fut.res$x[fut.res$if_fut == 1]
              
              #save marker.all.status from with efficacy/inferiority indicator
              marker.all.status_temp <- marker.all.status
              
              marker.all.status <- ifelse(marker.all %in% marker.fut & marker.all %in% marker_check_fut, 1, 0)
              
              #overlaid futility indicator with efficacy and inferiority indicator
              marker.all.status <- marker.all.status + marker.all.status_temp
              
            }
            
            if(j == 1){
              # at the first interim is the only time guarantee that marker.all has the marker values from the unselected population
              marker.all.dat <- data.frame(marker = marker.all, status = marker.all.status)
              marker.all.dat$end_time <-as.numeric(ifelse(marker.all.dat$status != 0, j, NA))
              marker.all.dat$status <-as.numeric(ifelse(marker.all.dat$status == 0, NA, marker.all.dat$status))
              
              
            }else{
              
            
              marker.all.dat.temp <- data.frame(marker = marker.all, status = marker.all.status) %>% distinct(marker, .keep_all = TRUE)
              marker.all.dat.temp$end_time <-as.numeric(ifelse(marker.all.dat.temp$status != 0, j, NA))
              
              marker.dat1 <- merge(marker.all.dat, marker.all.dat.temp, by = c("marker"), all.x = TRUE)%>% 
                mutate(status = coalesce(status.x,status.y)) %>% 
                mutate(end_time = coalesce(end_time.x,end_time.y)) %>% 
                dplyr::select(marker,status,end_time)%>%
                arrange(marker)
              
              marker.all.dat <-marker.dat1
              
              marker.all.dat$status <- as.numeric(ifelse(marker.all.dat$status == 0, NA, marker.all.dat$status))
              
            }

      
            
            # marker.rest <-marker.all[(length(res_data$x)+1):n]
            # marker.rest.status <-marker.all.status[(length(res_data$x)+1):n]
            # 
            # marker.rest.temp <- data.frame(marker = marker.rest,status = marker.rest.status) %>% 
            #   filter(status != 0)%>% 
            #   distinct(marker, .keep_all = TRUE)
            # 
            # marker.rest.dat <- rbind(marker.rest.dat,marker.rest.temp)
            # 
            # marker.stop <- unique(c(marker.stop,marker.all[marker.all.status != 0]))  #stopped marker values
            # marker.continuing <- unique(marker.all[marker.all.status == 0]) #continuing marker values
            
            marker.stop <- unique(marker.all.dat$marker[!is.na(marker.all.dat$status)])
            marker.continuing <- unique(marker.all.dat$marker[is.na(marker.all.dat$status)])
            
            # removed marker.rest
            marker.enrolled <- marker.all[1:length(res_data$x)]
            marker.enrolled.status <- marker.all.status[1:length(res_data$x)]
            marker.endtime <-marker.all.dat$end_time[1:length(res_data$x)]
            
            
            enroll.dat <- data.frame(x = marker.enrolled, status = marker.enrolled.status) #with duplicates and in the enrolled order
            
            enroll.res <- enroll.dat %>% group_by(x,status) %>% 
              count()%>% as.data.frame() #no dup, with count
            
           
            
           # when all markers are evaluated/ early stop
           if(length(marker.continuing) == 0) {
             
             samplesize <- sum(enroll.res$n)
             enroll.res$if_end <- rep(0,nrow(enroll.res)) #indicator if the early stopped
             
             marker.all.dat$if_end <-rep(0,nrow(marker.all.dat))
             
             marker.all.dat <- marker.all.dat %>% distinct(marker, .keep_all = TRUE)
             
             return(list(output = enroll.res, samplesize = samplesize, marker.stop= marker.stop, marker.continuing = marker.continuing,
               marker.all.dat = marker.all.dat, res_data = res_data ))
             
           } else{
             
             # adaptive randomization to next block
             
             previous.j <- j
             
             n.rem <- n-length(res_data$x)
             
             x.rem <- sample(x = marker.continuing, n.rem, replace = TRUE)
             
             all_marker <- c(all_marker[1:(check_num[j])], x.rem)
             
             x.next_block <- all_marker[(check_num[j]+1):check_num[j+1]]
             
             trt.next_block <- adaptive_randomization(x.eval = x.next_block, n = n, previous.n = length(res_data$x),
               tuning = tuning, 
               random_e =random_e,
               B_matrix = B_matrix, LZ_matrix = LZ_matrix, previous_ref = previous_ref, 
               g_est= g_est, h_est = h_est)
             
             marker.next_block <- 1:length(x.next_block)
             
             rem_data <- sim_data_rem(x.rem = x.rem, n.rem = n.rem, ref = ref, FUN.ctl =FUN.ctl, FUN.exp =FUN.exp)
            
             next_block_dat_temp <- data.frame(x = x.next_block, trt = trt.next_block$output, marker.order = marker.next_block)
             
             next_block_dat <- merge(next_block_dat_temp, rem_data, by = c("x",  "trt","marker.order"), all.x = TRUE, all.y = FALSE) %>% 
               dplyr::select(y, x, trt,marker.order) %>% arrange(marker.order)
             
             
             res_data <- rbind(res_data,next_block_dat)
             
             }

          }else{
            #when j == K only check efficacy
            
            # marker.j <- res_data$x[(check_num[previous.j]+1):check_num[j]]
            # marker.j.status <- rep(0, length(marker.j))  
            # fit.diff <- fit.diff.all %>% filter(x %in% marker.j) %>% unique()

            marker.j <- res_data$x[(check_num[previous.j]+1):check_num[j]]
            
            fit.diff <- fit.diff.all %>% distinct(x, .keep_all = TRUE)
              
            # enrolled previous block + new enrolled + will be enrolled in future
          
            
            marker.all <- res_data$x
            marker.all.status <- c(marker.enrolled.status,rep(0, length(marker.j)))
            if_end <- rep(0,length(marker.all))
            
            marker.notin <- marker.continuing[!(marker.continuing %in% marker.all)]
            
            if(length(marker.notin) != 0){
              
              marker.notin.status <- rep(0,length(marker.notin))
              
            }
            
            
            
            if(any(fit.diff$lwr_pt>0)){
              
              #if efficacy, labelled as 2
              #if not efficacy, labelled as 1
              #change if_end to 1 accordingly
              
              marker_eff <- fit.diff[fit.diff$lwr_pt>0,]$x
              marker_noeff <-fit.diff[fit.diff$lwr_pt<=0,]$x
              
              marker.all.status_temp <- marker.all.status
              
              if(length(marker_noeff) == 0){
                
                marker.all.status <- ifelse(marker.all.status == 0, 2, marker.all.status)
                
                marker.all.status <- marker.all.status + marker.all.status_temp
                
                if(length(marker.notin) != 0){
                  
                  marker.notin.status <- rep(2,length(marker.notin))
                  
                }
                
                
              }
              else{
              
                for(k in 1:length(marker.all)) {
                  
                  marker.all.status[k] <- ifelse(min(abs(marker.all[k]-marker_eff)) <= min(abs(marker.all[k]-marker_noeff)) & 
                      marker.all.status[k] == 0, 2, 1)
                  
                  #any marker with previous status 0 will be evaluate and labelled into 1 or 2 
                  if_end[k] <- ifelse(marker.all.status_temp[k] == 0, 1, 0)
                  
                }
              
               marker.all.status <- marker.all.status + marker.all.status_temp
               
               if(length(marker.notin) != 0){
                 
                 for(k in 1:length(marker.notin)) {
                   
                   marker.notin.status[k] <- ifelse(min(abs(marker.notin[k]-marker_eff)) <= min(abs(marker.notin[k]-marker_noeff)), 2, 1)
                   
                   
                 }
               }
              
              }

            }else{
              
              if_end <- ifelse(marker.all.status== 0, 1, 0)
              
              marker.all.status <- ifelse(marker.all.status == 0, 1, marker.all.status)
              
              if(length(marker.notin) != 0){
                
                marker.notin.status <- rep(1,length(marker.notin))
              }
              
            }
            
            
            enroll.dat <- data.frame(x = marker.all, status = marker.all.status, 
                    if_end = if_end) #with duplicates and in the enrolled order
            
            enroll.res <- enroll.dat %>% group_by(x,status,if_end) %>% 
              count()%>% as.data.frame()
            
            samplesize <- sum(enroll.res$n)
          
            
            marker.all.dat.temp <- data.frame(marker = marker.all, status = marker.all.status) %>% distinct(marker, .keep_all = TRUE)
            
            marker.dat1 <- merge(marker.all.dat, marker.all.dat.temp, by = c("marker"), all.x = TRUE)%>% 
              mutate(status = coalesce(status.x,status.y)) %>% 
              dplyr::select(marker,status,end_time)%>%
              arrange(marker)
            
            marker.dat1$if_end <- as.numeric(ifelse(is.na(marker.all.dat$status), 1, 0))
            marker.dat1$end_time <- as.numeric(ifelse(is.na(marker.dat1$end_time), j, marker.dat1$end_time))
            
            if(length(marker.notin) != 0){
              
              marker.notin.dat <- data.frame(marker =  marker.notin, status =  marker.notin.status, if_end = rep(1,length(marker.notin)),
                end_time = rep(as.numeric(j),length(marker.notin))) %>%
                distinct(marker, .keep_all = TRUE)
              
              marker.dat2 <- merge(marker.dat1, marker.notin.dat, by = c("marker"), all.x = TRUE)%>% 
                mutate(status = coalesce(status.x,status.y)) %>% 
                mutate(if_end = coalesce(if_end.x,if_end.y)) %>% 
                mutate(end_time = coalesce(end_time.x,end_time.y)) %>% 
                dplyr::select(marker,status,if_end,end_time)%>%
                arrange(marker)
              
              marker.all.dat <-marker.dat2
              
            }else{
              
              marker.all.dat <- marker.dat1
              
            }
            
            marker.all.dat <- marker.all.dat %>% distinct(marker, .keep_all = TRUE)
         
            return(list(output = enroll.res, samplesize = samplesize, marker.stop= marker.stop, marker.continuing = marker.continuing,
              marker.all.dat = marker.all.dat, res_data = res_data ))
            
          }
        }
    }
  }
 }
}
