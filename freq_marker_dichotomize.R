
Freq_detect_marker <- function(data, 				
  #####################################################################################################	 	
  ### data = a data frame with the observations arranged by row, and including the columns below (in the following order):
  ###        1) futime: time to event or right censoring, whichever happens first
  ###        2) status: the censoring status indicator (=1 for event; =0 for right censored)
  ###        3) x: the value of the continuous marker
  ###        4) trt: the group indicator (=1 for experimental group; =0 for control group)
  ##########
  x.all, 								# the continuous marker values of the entire dataset
  m.prev.lwr = 0.05,						# the lower bound of marker prevalence allowed for the grid search of cutpoint
  m.increment = 0.01, 						# the increment of the sequence of marker cutpoint searched 
  p_int = 0.01){							# minimum p-value cutoff for the interaction term to decide whether to dichotomize the marker
  #####################################################################################################
  
  
  m.grid = quantile(data$x, seq(m.prev.lwr, 1-m.prev.lwr, m.increment))
  p.value = rep(NA, length(m.grid))
  
  for (l in 1:length(m.grid)) {
    
    m.cutpoint = m.grid[l]
    
    data.l = data
    data.l$xd = as.numeric(data.l$x >= m.cutpoint)
    
    fit.l <- glm(data = data.l, y ~ trt +  xd + trt*xd, family = binomial())
    coef.fit.l <- coef(summary(fit.l))
    
    #check if number of patients in each treatment is too imbalanced
    if(nrow(coef.fit.l) < 4) {
      
      p.value[l] <- 1
      
    }else{
      p.value[l] <- coef.fit.l[4,4]
    }
    
    #p.value[l] <- coef.fit.l[3,5]
 
    
  }
  
  if(min(p.value) < p_int) {
    
    m.cutpoint = m.grid[which.min(p.value)]
    
    data.l = data
    data.l$xd = as.numeric(data.l$x >= m.cutpoint)
    
    # fit.l <- coxph(Surv(futime, status) ~ trt + xd + trt*xd, data=data.l)
    # coef.fit.l <- coef(summary(fit.l))
    fit.l <- glm(data= data.l, y ~ trt +  xd + trt*xd, family = binomial())
    coef.fit.l <- coef(summary(fit.l))
    
    # if(coef.fit.l[3,1] < 0) {
    if(coef.fit.l[4,1] > 0) {
      marker = as.numeric(data.l$x >= m.cutpoint)
      marker.all = as.numeric(x.all >= m.cutpoint)
      pos.marker.range = c(m.cutpoint, max(x.all))
      neg.marker.range = c(min(x.all), m.cutpoint-m.increment)
    }
    
    else {
      marker = as.numeric(data.l$x <= m.cutpoint)
      marker.all = as.numeric(x.all <= m.cutpoint)
      pos.marker.range = c(min(x.all),m.cutpoint)
      neg.marker.range = c(m.cutpoint-m.increment, max(x.all))
    }
    
    return(list(dicho=TRUE, cutpoint=m.cutpoint, pvalue=min(p.value), marker=marker, marker.all=marker.all, 
      pos.marker.range = pos.marker.range, neg.marker.range = neg.marker.range))
  }
  
  else {
    return(list(dicho=FALSE, pvalue=min(p.value)))
  }
  
}

