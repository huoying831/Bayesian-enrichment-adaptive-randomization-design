

setwd("/home1/yuetu/freq_design")
# source("/Users/ty/Documents/paper/biomarker/randomization/continuous_marker_binary/local_run_Jul_2023/freq_design.R")
# source("/Users/ty/Documents/paper/biomarker/randomization/continuous_marker_binary/local_run_Jul_2023/freq_marker_dichotomize.R")
args <- commandArgs(TRUE)
phen <- args[1]


source("/home1/yuetu/freq_design/freq_design.R")
# 
# case = 1
# phen = 1
setting = as.numeric(phen)
#setting = as.numeric(case)

# Nrep = 100
# seed = 12342423		            # Number of iterations and seed to reproduce results
n = 500 							   	 	              # total number of patients expected to enroll (both arms combined)
percent.exp = 0.5                        # randomization ratio to the experimental group
check_num = c(250,500)
m.increment = 0.01
p_int = rep(0.05, length(check_num))						   # a vector of minimum p-value cutoffs for the treatment-by-marker interaction to decide whether to dichotomize (length=number of checks)
p_eff_all = rep(0.025, length(check_num))  					 # a vector of overall cohort one-sided log-rank test p-value cutoffs to claim efficacy (length=number of checks)
p_eff_grp = rep(0.025, length(check_num))						   # a vector of marker subgroup one-sided log-rank test p-value cutoffs to claim efficacy (length=number of checks)
p_fut_all = rep(0.1, length(check_num)-1)						 # a vector of overall cohort thresholds of HR above which to claim futility (length=number of checks - 1)
p_fut_grp = rep(0.1, length(check_num)-1)						 # a vector of marker subgroup thresholds of HR above which to claim futility (length=number of checks - 1)
stop_eff = TRUE                                      # a logical indictor of whether to allow early stopping due to efficacy		 	                                              		
theta_alter = 0.01
m.prev.lwr = 0.2



fun_control = list(function(x) {0*x-2.95}, function(x) {0*x-2.95}, function(x) {2.1*x-2.95}, function(x) {0*x-2.95}, function(x) {0*x-2.95},
  function(x) {0*x-2.95},function(x) {0*x-2.95},function(x) {0*x-2.95},function(x) {0*x-2.95},function(x) {0*x-2.95},
  function(x) {0*x -0.85},
 function(x) {1.22*x-2.95}
)
fun_experiment = list(function(x) {0*x-2.95}, function(x) {0*x -0.85},function(x) {2.1*x-2.95},function(x) {2.1*as.numeric(x>=0.5)-2.95},
  function(x) {
    x2 <- exp(30*(x-0.5))
    y<- 2.1*x2/(1+x2)-2.95
    return(y)
  },
  function(x) {2.1*x-2.95},
  function(x) {-0.85-exp(-7*(x-0.106))},
  # function(x) {2.1*sqrt(x)},
  function(x) {exp(1.13*x)-3.95},
  function(x) {
    if(x<=0.5) {
      -0.85-2.1*exp(30*(x-0.2))/(1+exp(30*(x-0.2)))
    }
    else {
      -2.95+2.1*exp(30*(x-0.8))/(1+exp(30*(x-0.8)))
    }
  },
  function(x) {
    if(x<=0.5) {
      -2.95+2.1*exp(30*(x-0.2))/(1+exp(30*(x-0.2)))
    }
    else {
      -0.85-2.1*exp(30*(x-0.8))/(1+exp(30*(x-0.8)))
    }
  },
  function(x) {0*x-2.95},
 function(x) {
    x2 <- exp(30*(x-0.5))
    y<- 1.33*x2/(1+x2)-2.95 + 1.22*x
    return(y)
  }
)

for(seed in 1:1000){
  if(file.exists(paste0("/home1/yuetu/freq_design/output/case",setting,"/temp",seed,".RData"))){
    quit(save="no")}else{ 
  FUN.ctl <-fun_control[[setting]]
  FUN.exp <- fun_experiment[[setting]]
  
  result <- Freq_adaptive_design(
    FUN.ctl = FUN.ctl,
    FUN.exp = FUN.exp,
    seed = seed,		            # Number of iterations and seed to reproduce results
        n = n, 							   	 	              # total number of patients expected to enroll (both arms combined)
    check_num = check_num,
    m.increment = m.increment,
    p_int = p_int,
    p_eff_all = p_eff_all,
    p_eff_grp = p_eff_grp,
    p_fut_all = p_fut_all,
    p_fut_grp = p_fut_grp,
    stop_eff = stop_eff,
    theta_alter = theta_alter,
    m.prev.lwr = m.prev.lwr
      )
      
  # setwd("/Users/ty/Documents/paper/biomarker/randomization/continuous_marker_binary/local_run_Jul_2023/")
  # save(result, file=paste0("temp",seed,".RData"))
  setwd(paste0("/home1/yuetu/freq_design/output"))
  save(result, file=paste0("case",setting,"seed",seed,".RData"))
  
    }
  }
