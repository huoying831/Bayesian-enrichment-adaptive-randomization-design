
Freq_adaptive_design <- function(
  n,
  seed,
  percent.exp = 0.5,   
  FUN.ctl = function(x) {0*x-2.95},              # mapping function of the log hazard as a function of marker for the control group
  FUN.exp = function(x) {0*x-2.95},   
  #####################################################################################################	 	
  ### surdata = a data frame with the observations arranged by row, and including the columns below (in the following order) ###
  ### treatment arm indicator trt, entry dates enter, continuous marker values x, event times T, and observed TTE dates T_date. ###  
  ##########		     							   		   
  # n_events = nrow(surdata),						             # the number of events expected to be observed at the end of the trial
  # m.prev.lwr = 0.2,						   	                 # the lower bound of marker prevalence allowed for the grid search of cutpoint
  check_num = c(250,500), 			   	     # interim analysis times based on % of observed events for T
  m.increment = 0.01, 							                 # the increment of the sequence of marker cutpoint searched 
  p_int = rep(0.01, length(check_num)),							   # a vector of minimum p-value cutoffs for the treatment-by-marker interaction to decide whether to dichotomize (length=number of checks)
  p_eff_all = rep(0.01, length(check_num)),  					 # a vector of overall cohort one-sided log-rank test p-value cutoffs to claim efficacy (length=number of checks)
  p_eff_grp = rep(0.1, length(check_num)),						   # a vector of marker subgroup one-sided log-rank test p-value cutoffs to claim efficacy (length=number of checks)
  p_fut_all = rep(0.1, length(check_num)-1),						 # a vector of overall cohort thresholds of HR above which to claim futility (length=number of checks - 1)
  p_fut_grp = rep(0.1, length(check_num)-1),						 # a vector of marker subgroup thresholds of HR above which to claim futility (length=number of checks - 1)
  stop_eff = TRUE,                                      # a logical indictor of whether to allow early stopping due to efficacy		 	                                              		
  theta_alter = 0.05,
  m.prev.lwr = 0.2
  ) { 							                 
  
  
  
  ### Code here to load in the functions that will be called later
  source("/home1/yuetu/freq_design/freq_marker_dichotomize.R")  
  
  set.seed(seed)
  ### Number of analyses (including interim and final analysis)
  K <- length(check_num)
  
  
  ### Create vectors to save results
  # Initialize indicator vectors for stopping (overall, cohort), endpoint, and reasons
  stop.eff.pos <- rep(0, K)
  stop.eff.neg <- rep(0, K)
  stop.eff.all <- rep(0, K)
  
  stop.fut.pos <- rep(0, K)
  stop.fut.neg <- rep(0, K)
  stop.fut.all <- rep(0, K)
  
  stop.inf.pos <- rep(0, K)
  stop.inf.neg <- rep(0, K)
  stop.inf.all <- rep(0, K)
  
  # Initialize indicator vectors for continuing (rather than stopping early)
  continue.pos <- rep(0, K-1)
  continue.neg <- rep(0, K-1)
  continue.all <- rep(0, K-1)
  
  
  # Initialize scalar to store final Bayesian sample sizes 
  sampsize <- NA
  sampsize.pos <- NA
  sampsize.neg <- NA
  
  # Initialize scalar to store marker prevalence and marker status for each patient if dichotomized
  marker.prevalence <- NA
  marker.status <- rep(NA, n)
  
  
  
  ### Initialize the trial stopping indicator (of either cohort if marker is dichotomized), which is updated at each interim check 
  stop.any <- FALSE
  stop.pos <- FALSE
  stop.neg <- FALSE
  
  
  
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
      x=c(x.ctl,x.exp,x.ctl,x.exp), trt= trt)
    
    return(list(sim_all_data = all_data, x = all_data$x[c(1:n)], ref = ref, marker_seq = marker_seq))
  }
  
  
  sim_rem_data <- function(n = 300,         								      #remaining sample (both arm combined)
    percent.exp = 0.5,                                             # randomization ratio to the experimental group
    #marker.specs = list(dist="beta", param=c(1,1),range=c(0,1), ref="median"),        	# a list giving the distribution, range, and reference point of the continuous marker x                                     
    FUN.ctl = function(x) {0*x-2.95},                                   # mapping function of the probability of response as a function of marker x for the control group g(x)
    FUN.exp = function(x) {0*x-2.95},                                 # mapping function of the probability of response as a function of marker x for the experimental group g(x)+h(x)
    rem.marker.range = c(0,1),
    ref
    )                                  
  {
    # Compute number of patients in each arm
    nsize.exp <- round(n*percent.exp)
    nsize.ctl <- n-nsize.exp
    
    marker_seq <- round(seq(rem.marker.range[1], rem.marker.range[1], by = 0.01),6)
    x.ctl <- sample(marker_seq,nsize.ctl,replace = TRUE)
    x.exp <- sample(marker_seq,nsize.exp,replace = TRUE)
    
    # Calculate the probability of response, as a function of continuous marker x, relative to the reference point in the control group
    x.ctl_eff.ctl <- unlist(lapply(x.ctl,FUN.ctl))-FUN.ctl(ref)
    x.exp_eff.exp <- unlist(lapply(x.exp,FUN.exp))-FUN.ctl(ref)
    
    # Generate endpoint y based on the function of marker value
    x.clt_y.ctl <- rbinom(nsize.ctl,1,exp(x.ctl_eff.ctl)/(1+exp(x.ctl_eff.ctl)))	     # the control group
    x.exp_y.exp <- rbinom(nsize.exp,1,exp(x.exp_eff.exp)/(1+exp(x.exp_eff.exp)))	     # the experimental group
    
    # Set up treatment arm indicators
    trt <- c(rep(0, nsize.ctl), rep(1, nsize.exp))					# Treatment assignment 
    
    
    all_data <- data.frame(y=c(x.clt_y.ctl, x.exp_y.exp), 
      x=c(x.ctl,x.exp), trt= trt)
    
    return(list(sim_all_data = all_data, x = all_data$x, marker_seq = marker_seq))
  }
  
  #use sim_data
  sim_data <- sim_data(n=n, percent.exp=percent.exp, FUN.ctl=FUN.ctl, FUN.exp=FUN.exp)
  all_res_data <- sim_data$sim_all_data
  all_marker <- sim_data$x
  ref <- sim_data$ref
  marker_seq <- sim_data$marker_seq
  
  no_rand_data <- all_res_data[1:n,]
  no_rand_data <- no_rand_data[sample(1:n),]
 

  #############################################
  #############################################
  
  #############################################
  ###### Loop over interim analysis times #####
  #############################################
  
  
  for(j in 1:K){
    
    #####################################################
    #####################################################
    ################# Interim Analyses ##################
    #####################################################
    #####################################################
    
    
    #####################################################
    ###### Add estimation / interim decision code here ##### 
    
    if(!stop.any) {
      
      data.j = data.frame(no_rand_data[1:(check_num[j]),])
      
      
      ##### Call the function to decide if and how to dichotomize the marker based on data available at this time
      result.j = Freq_detect_marker(data=data.j, x.all= marker_seq, m.prev.lwr = m.prev.lwr, m.increment = m.increment, p_int = p_int[j])
      
      
      ##### If dichotomize the marker #####
      if( result.j$dicho ) {
        
        # dichotomize the continuous marker for individuals enrolled at this time
        data.j$marker <- result.j$marker
        
        marker.pos <- marker_seq[result.j$marker.all == 1] 
        marker.neg <- marker_seq[result.j$marker.all == 0] 
        
        samp.pos <- sum(data.j$marker==1)
        samp.neg <- sum(data.j$marker==0)
        
        pos.marker.range <-result.j$pos.marker.range
        neg.marker.range <-result.j$neg.marker.range
        
        # subset the data available at this time based on marker positivity and arm
        data.pos <- data.j[data.j$marker==1, ]
        data.neg <- data.j[data.j$marker==0, ]
        
        test.pos <- glm(data = data.pos, y ~ trt, family = binomial())			# the marker positive cohort
        pvalue.pos <- coef(summary(test.pos))[2,4]
        z.pos <- coef(summary(test.pos))[2,3]
        se.pos <- coef(summary(test.pos))[2,2]
        est.pos <- coef(summary(test.pos))[2,1]

        test.neg <- glm(data = data.neg, y ~ trt, family = binomial())			# the marker negative cohort
        pvalue.neg <- coef(summary(test.neg))[2,4]
        z.neg <- coef(summary(test.neg))[2,3]
        se.neg <-coef(summary(test.neg))[2,2]
        est.neg <- coef(summary(test.neg))[2,1]
        
        
        ###### If this is the final check, skip the futility analysis ######
        if(j==K) {
          stop.eff.pos[j] <- as.numeric(pvalue.pos < p_eff_grp[j] & est.pos > 0)
          stop.eff.neg[j] <- as.numeric(pvalue.neg < p_eff_grp[j] & est.neg > 0)
          stop.inf.pos[j] <- as.numeric(pvalue.pos < p_eff_grp[j] & est.pos < 0)
          stop.inf.neg[j] <- as.numeric(pvalue.neg < p_eff_grp[j] & est.neg < 0)
          stop.fut.pos[j] <- as.numeric(!(as.logical(stop.eff.pos[j])|as.logical(stop.inf.pos[j])))
          stop.fut.neg[j] <- as.numeric(!(as.logical(stop.eff.neg[j])|as.logical(stop.inf.neg[j])))  
          sampsize <- nrow(data.j)
          sampsize.pos <- samp.pos
          sampsize.neg <- samp.neg
          marker.prevalence <- mean(result.j$marker.all==1)
          marker.status <- result.j$marker.all
        }
        
        ##### If this is not the final check, perform efficacy, inferiority and futility analysis #####	
        else {
          
          if(stop_eff & (pvalue.pos < p_eff_grp[j] & est.pos > 0)) {					# stop the marker positive cohort due to efficacy
            stop.pos <- TRUE 
            stop.eff.pos[j] <- 1
          }
          
          else if(stop_eff & (pvalue.pos < p_eff_grp[j] & est.pos < 0)) {			# stop the marker positive cohort due to inferiority
            stop.pos <- TRUE
            stop.inf.pos[j] <- 1
          }
          
          else {
            # test for marker positive cohort for futility using conditional power
            #https://www.ncss.com/wp-content/themes/ncss/pdf/Procedures/PASS/Conditional_Power_and_Sample_Size_Reestimation_of_Tests_for_the_Difference_Between_Two_Proportions.pdf
            #set expected difference under the alternative hypothesis as 0.05
            
            data.pos.j <- data.j[data.j$x %in% marker.pos,]
            pos.j <- nrow(data.pos.j)
            pos.ctl.j <- length(data.pos.j$trt == 0)
            pos.exp.j <- length(data.pos.j$trt == 1)
            pos.K <- n-nrow(data.j) + pos.j  #assume only enroll positive subgroup till end
            info.pos.j <- 1/se.pos^2*1/(1/pos.ctl.j + 1/pos.exp.j) 
            info.pos.K <-1/se.pos^2*1/(4/pos.K) #assume equal randomization
            fut.pos.z <- z.pos* sqrt(info.pos.j)-qnorm(1-p_eff_grp[j])*sqrt(info.pos.K) + theta_alter*(info.pos.K-info.pos.j)
            fut.pos <- pnorm(fut.pos.z)
            
            if( stop_eff & fut.pos < p_fut_grp[j]){
              stop.pos <- TRUE
              stop.fut.pos[j] <- 1
            }
            else {
              continue.pos[j] <- 1
              }
            }
      
          
          
          if(stop_eff & (pvalue.neg < p_eff_grp[j] & est.neg > 0)) {					# stop the marker negative cohort due to efficacy
            stop.neg <- TRUE 
            stop.eff.neg[j] <- 1
          }
          
          else if(stop_eff & (pvalue.neg < p_eff_grp[j] & est.neg < 0)){			# if the marker negative cohort fails to stop early for efficacy, test futility
            stop.neg <- TRUE
            stop.inf.neg[j] <- 1
          }
            
            else {
              
              data.neg.j <- data.j[data.j$x %in% marker.neg,]
              neg.j <- nrow(data.neg.j)
              neg.ctl.j <- length(data.neg.j$trt == 0)
              neg.exp.j <- length(data.neg.j$trt == 1)
              neg.K <- n - nrow(data.j)+neg.j #assume only enroll negative subgroup till end
              info.neg.j <- 1/se.neg^2*1/(1/neg.ctl.j + 1/neg.exp.j) #assume equal randomization
              info.neg.K <-1/se.neg^2*1/(4/neg.K)
              fut.neg.z <- z.neg* sqrt(info.neg.j)-qnorm(1-p_eff_grp[j])*sqrt(info.neg.K) + theta_alter*(info.neg.K-info.neg.j)
              fut.neg <- pnorm(fut.neg.z)
              
              if( stop_eff & fut.neg < p_fut_grp[j]){
                stop.neg <- TRUE
                stop.fut.neg[j] <- 1
              } else {
              continue.neg[j] <- 1
            }
          }
          
          
          ##### Determine if to stop (any or both) marker cohorts at this check
          stop.any <- (stop.pos | stop.neg)
          
          
          
          if(stop.any) {
            
            marker <- result.j$marker.all
            marker.prevalence <- mean(marker==1)
            marker.status <- marker

            if(stop.pos & stop.neg) {
              sampsize <- nrow(data.j)
              sampsize.pos <- samp.pos
              sampsize.neg <- samp.neg							
              break
            } else if(stop.neg){ #continuous only positive subgroup
              
              rem.n <- n-nrow(data.j)
              sim_rem_data <- sim_rem_data(n = rem.n, percent.exp = percent.exp, FUN.ctl = FUN.ctl,
                FUN.exp = FUN.exp, ref=ref, rem.marker.range = pos.marker.range)
              
              rem.data <- sim_rem_data$sim_all_data
              rem.data <- rem.data[sample(1:nrow(rem.data)),]
              rem.data$marker <- rep(as.numeric(1),nrow(rem.data))
              no_rand_data <- rbind(data.j,rem.data)
              marker.status <- no_rand_data$marker
              
            } else if(stop.pos) { #continuous only negative subgroup
              rem.n <- n-nrow(data.j)
              sim_rem_data <- sim_rem_data(n = rem.n, percent.exp = percent.exp, FUN.ctl = FUN.ctl,
                FUN.exp = FUN.exp, ref=ref, rem.marker.range = neg.marker.range)
              
              rem.data <- sim_rem_data$sim_all_data
              rem.data <- rem.data[sample(1:nrow(rem.data)),]
              rem.data$marker <- rep(as.numeric(0),nrow(rem.data))
              no_rand_data <- rbind(data.j,rem.data)
              marker.status <- no_rand_data$marker
              
            }
          }
        }
      }
      
      
      ##### If not to dichotomize the marker #####
      else {
        
        test.all <- glm(data = data.j, y ~ trt, family = binomial())
        pvalue.all <- coef(summary(test.all))[2,4]
        z.all <-coef(summary(test.all))[2,3]
        se.all <-coef(summary(test.all))[2,2]
        est.all <- coef(summary(test.all))[2,1]
        
        
        ###### do log-rank test and fit proportional hazard model in the overall cohort ######
        # test.all <- survdiff(Surv(futime, status) ~ trt, data=data.j)		
        # idx.all <- which(names(test.all$n) == "trt=1")
        # z.all <- (test.all$obs[idx.all]-test.all$exp[idx.all])/sqrt(test.all$var[idx.all,idx.all])
        # pvalue.all <- pnorm(z.all, lower.tail = TRUE)
        # ph.all <- coxph(Surv(futime, status) ~ trt, data=data.j)
        # 
        
        ###### If this is the final check, skip the futility analysis ######
        if(j==K) {
          stop.eff.all[j] <- as.numeric(pvalue.all < p_eff_all[j] & est.all > 0)
          stop.inf.all[j] <- as.numeric(pvalue.all < p_eff_all[j] & est.all < 0)
          stop.fut.all[j] <- as.numeric(!(as.logical(stop.eff.all[j])|as.logical(stop.inf.all[j])))  
          sampsize <- nrow(data.j)
          sampsize.pos <- NA
          sampsize.neg <- NA	
          marker.prevalence <- NA
          marker.status <- rep(NA, n)
        }
        
        ##### If this is not the final check, perform efficacy inferiority and futility analysis #####	
        else {
          if(stop_eff & (pvalue.all < p_eff_all[j]& est.all > 0)) {					# stop the overall cohort due to efficacy
            stop.any <- TRUE 
            stop.eff.all[j] <- 1
            sampsize <- nrow(data.j)
            sampsize.pos <- NA
            sampsize.neg <- NA	
            marker.prevalence <- NA	
            marker.status <- rep(NA, n)						
            break
          }
          
          else if(stop_eff & (pvalue.all < p_eff_all[j]& est.all < 0)) {					# if not, test futility for the overall cohort
            stop.any <- TRUE 
            stop.inf.all[j] <- 1
            sampsize <- nrow(data.j)
            sampsize.pos <- NA
            sampsize.neg <- NA
            marker.prevalence <- NA	
            marker.status <- rep(NA, n)							
            break
          }
          else {
            
            all.j <- nrow(data.j)
            all.K <- n
            all.ctl.j <- length(data.j$trt == 0)
            all.exp.j <- length(data.j$trt == 1)
            info.all.j <- 1/se.all^2*1/(1/all.ctl.j + 1/all.exp.j) #assume equal randomization
            info.all.K <-1/se.all^2*1/(4/all.K)
            fut.all.z <- z.all* sqrt(info.all.j)-qnorm(1-p_eff_all[j])*sqrt(info.all.K) + theta_alter*(info.all.K-info.all.j)
            fut.all <- pnorm(fut.all.z)
            
            if( stop_eff & fut.all < p_fut_all[j]){
              stop.any <- TRUE
              stop.fut.all[j] <- 1
              sampsize <- nrow(data.j)
              sampsize.pos <- NA
              sampsize.neg <- NA	
              marker.prevalence <- NA	
              marker.status <- rep(NA, n)						
              break
              
              }
            
              else {
              continue.all[j] <- 1
            }
          }
          
        }
      }
      
    
    
    }else if (stop.pos & (!stop.neg)) {							# only the marker negative cohort continues
      
      # extract the marker negative cohort based on the dichotomization rule
      data.j <- data.frame(no_rand_data[1:(check_num[j]),])
      data.neg.j <- data.j[data.j$marker == 0,]
  
      samp.neg <- nrow(data.neg.j)
      
      
      ###### Do log-rank test and fit proportional hazard model in the marker-negative cohort ######
      # test.neg <- survdiff(Surv(futime, status) ~ trt, data=data.j)		
      # idx.neg <- which(names(test.neg$n) == "trt=1")
      # z.neg <- (test.neg$obs[idx.neg]-test.neg$exp[idx.neg])/sqrt(test.neg$var[idx.neg,idx.neg])
      # pvalue.neg <- pnorm(z.neg, lower.tail = TRUE)
      # ph.neg <- coxph(Surv(futime, status) ~ trt, data=data.j)
      
      test.neg <- glm(data = data.neg.j, y ~ trt, family = binomial())			# the marker negative cohort
      pvalue.neg <- coef(summary(test.neg))[2,4]
      z.neg <- coef(summary(test.neg))[2,3]
      se.neg <-coef(summary(test.neg))[2,2]
      est.neg <- coef(summary(test.neg))[2,1]
      
      
      ###### If this is the final check, skip the futility analysis ######
      if(j==K) {
        stop.eff.neg[j] <- as.numeric(pvalue.neg < p_eff_grp[j] & est.neg > 0)
        stop.inf.neg[j] <- as.numeric(pvalue.neg < p_eff_grp[j] & est.neg < 0)
        stop.fut.neg[j] <- as.numeric(!(as.logical(stop.eff.neg[j])|as.logical(stop.inf.neg[j])))  
        sampsize <- samp.pos + samp.neg
        sampsize.pos <- samp.pos
        sampsize.neg <- samp.neg
      }
      
      ##### If this is not the final check, perform efficacy and futility analysis #####	
      else {
        if(stop_eff & (pvalue.neg < p_eff_grp[j] & est.neg > 0)) {				# stop the marker negative cohort due to efficacy
          stop.neg <- TRUE 
          stop.eff.neg[j] <- 1
          sampsize <- samp.pos + samp.neg
          sampsize.pos <- samp.pos
          sampsize.neg <- samp.neg
          break
        }
        
        else if (stop_eff &(pvalue.neg < p_eff_grp[j] & est.neg < 0)) {				# if not, test futility for the marker negative cohort
          stop.neg <- TRUE
          stop.inf.neg[j] <- 1
          sampsize <- samp.pos + samp.neg
          sampsize.pos <- samp.pos
          sampsize.neg <- samp.neg
          break
        }
        else {
          
          neg.j <- nrow(data.neg.j)
          neg.ctl.j <- length(data.neg.j$trt == 0)
          neg.exp.j <- length(data.neg.j$trt == 1)
          neg.K <- n-nrow(data.j)+ neg.j
          info.neg.j <- 1/se.neg^2*1/(1/neg.ctl.j + 1/neg.pos.j) #assume equal randomization
          info.neg.K <-1/se.neg^2*1/(4/neg.K)
          fut.neg.z <- z.neg* sqrt(info.neg.j)-qnorm(1-p_eff_grp[j])*sqrt(info.neg.K) + theta_alter*(info.neg.K-info.neg.j)
          fut.neg <- pnorm(fut.neg.z)
          
          if( stop_eff & fut.neg < p_fut_grp[j]){
            stop.neg <- TRUE
            stop.fut.neg[j] <- 1
            sampsize <- samp.pos + samp.neg
            sampsize.pos <- samp.pos
            sampsize.neg <- samp.neg
            break
            
          }
          
          else {
            continue.neg[j] <- 1
          }
        }
      }
    }else if (stop.neg & (!stop.pos)) {							# only the marker positive cohort continues
      
      
      # extract the marker positive cohort based on the dichotomization rule
      data.j <- data.frame(no_rand_data[1:(check_num[j]),])
      data.pos.j <- data.j[data.j$marker == 1,]
      
      samp.pos <- nrow(data.pos.j)
      
      test.pos <- glm(data = data.pos.j, y ~ trt, family = binomial())			# the marker positive cohort
      pvalue.pos <- coef(summary(test.pos))[2,4]
      z.pos <- coef(summary(test.pos))[2,3]
      se.pos <- coef(summary(test.pos))[2,2]
      est.pos <- coef(summary(test.pos))[2,1]
      
      ###### If this is the final check, skip the futility analysis ######
      if(j==K) {
        stop.eff.pos[j] <- as.numeric(pvalue.pos < p_eff_grp[j] & est.pos > 0)
        stop.inf.pos[j] <- as.numeric(pvalue.pos < p_eff_grp[j] & est.pos < 0)
        stop.fut.pos[j] <- as.numeric(!(as.logical(stop.eff.pos[j])|as.logical(stop.inf.pos[j])))
        sampsize <- samp.pos + samp.neg
        sampsize.pos <- samp.pos
        sampsize.neg <- samp.neg
      }
      
      ##### If this is not the final check, perform efficacy and futility analysis #####	
      else {
        if(stop_eff & (pvalue.pos < p_eff_grp[j] & est.pos > 0)) {				# stop the marker positive cohort due to efficacy
          stop.pos <- TRUE 
          stop.eff.pos[j] <- 1
          sampsize <- samp.pos + samp.neg
          sampsize.pos <- samp.pos
          sampsize.neg <- samp.neg
          break
        }
        
        else if (stop_eff & (pvalue.pos < p_eff_grp[j] & est.pos < 0)) {				# if not, test futility for the marker negative cohort
          stop.pos <- TRUE
          stop.inf.pos[j] <- 1
          sampsize <- samp.pos + samp.neg
          sampsize.pos <- samp.pos
          sampsize.neg <- samp.neg
          break
        }
        else {
          
          pos.j <- nrow(data.pos.j)
          pos.ctl.j <- length(data.pos.j$trt == 0)
          pos.exp.j <- length(data.pos.j$trt == 1)
          pos.K <- n - nrow(data.j) + pos.j
          info.pos.j <- 1/se.pos^2*1/(1/pos.ctl.j + 1/pos.exp.j) 
          info.pos.K <-1/se.pos^2*1/(4/pos.K) #assume equal randomization
          fut.pos.z <- z.pos* sqrt(info.pos.j)-qnorm(1-p_eff_grp[j])*sqrt(info.pos.K) + theta_alter*(info.pos.K-info.pos.j)
          fut.pos <- pnorm(fut.pos.z)
          
          if( stop_eff & fut.pos < p_fut_grp[j]){
            stop.pos <- TRUE
            stop.fut.pos[j] <- 1
            sampsize <- samp.pos + samp.neg
            sampsize.pos <- samp.pos
            sampsize.neg <- samp.neg
            break
          }
          else {
            continue.pos[j] <- 1
          }
          
        }
      }
    }
    
    
    ######################################################                      
    
  }   #  end loop over interim checks
  
  
  ### Determine the final result for this trial (overall efficacy, subgroup efficacy, or overall futility)
  if(is.na(marker.prevalence)) {										# the marker is not dichotomized
    if(sum(stop.eff.all)>0) {										# overall efficacy
      res.final <- "eff_all"
    }
    else if(sum(stop.inf.all)>0){
      res.final <- "inf_all"										
    }else if(sum(stop.fut.all)>0){
      res.final <- "fut_all"										
    }
  }else {														# the marker is dichotomized
    if(sum(stop.eff.pos)>0 & sum(stop.eff.neg)>0) {							# efficacy in both marker subgroups
      res.final <- "eff_all"
    }
    else if (sum(stop.fut.pos)>0 & sum(stop.fut.neg)>0) {					# futility in either marker subgroup
      res.final <- "fut_all"
    }	
    else if (sum(stop.inf.pos)>0 & sum(stop.inf.neg)>0) {					# inferiority in either marker subgroup
      res.final <- "inf_all"
    }	
    
    else if (sum(stop.eff.pos)>0 & sum(stop.fut.neg)>0) {				
      res.final <- "eff_pos,fut_neg"
    }	
    else if (sum(stop.eff.pos)>0 & sum(stop.inf.neg)>0) {				
      res.final <- "eff_pos,inf_neg"
    }	
    else if (sum(stop.inf.pos)>0 & sum(stop.eff.neg)>0) {				
      res.final <- "inf_pos,eff_neg"
    }	
    else if (sum(stop.inf.pos)>0 & sum(stop.fut.neg)>0) {				
      res.final <- "inf_pos,fut_neg"
    }	
    else if (sum(stop.fut.pos)>0 & sum(stop.eff.neg)>0) {				
      res.final <- "fut_pos,eff_neg"
    }	
    else if (sum(stop.fut.pos)>0 & sum(stop.inf.neg)>0) {				
      res.final <- "fut_pos,inf_neg"
    }	
  }	
  
  
  return(list(stop.eff.pos=stop.eff.pos, stop.eff.neg=stop.eff.neg, stop.eff.all=stop.eff.all, 
    stop.fut.pos=stop.fut.pos, stop.fut.neg=stop.fut.neg, stop.fut.all=stop.fut.all,
    continue.pos=continue.pos, continue.neg=continue.neg, continue.all=continue.all, 
    res.final=res.final, sampsize=sampsize, sampsize.pos=sampsize.pos, sampsize.neg=sampsize.neg, 
    marker.prevalence=marker.prevalence, x=no_rand_data$x, marker.status=marker.status))   
  
} #end function




