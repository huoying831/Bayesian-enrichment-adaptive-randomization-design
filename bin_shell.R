

setwd("/home1/yuetu/attempt_design")
source("design_2_14_clean.R")

args <- commandArgs(TRUE)
phen <- args[1]
case <- args[2]


seed = as.numeric(phen)
setting = as.numeric(case)

# Nrep = 100
# seed = 12342423		            # Number of iterations and seed to reproduce results
n = 500 							   	 	              # total number of patients expected to enroll (both arms combined)
percent.exp = 0.5                        # randomization ratio to the experimental group
marker.specs = list(dist="beta",param=c(1,1),range=c(0,1),ref="median") 		# a list giving the distribution, range, reference point of continuous marker, and whether to include right censored subjects in calculation of interior knots location                                   
# mapping function of the log hazard as a function of marker for the experimental group
 MCMC.specs = list(burnin=25000,B=10000,thin=50,multiplier=0.5)         		# Bayes: a list giving the MCMC specs for spline modeling, including burn-in length, number of posterior samples, thinning parameter and scaling factor for the size of random walk proposal
# MCMC.specs = list(burnin=20000,B=2000,thin=20,multiplier=0.5)

MCMC.specs.fut = list(burnin=20000,B=2000,thin=20,multiplier=0.5)
# MCMC.specs = list(burnin=2500,B=1000,thin=5,multiplier=0.5)         		# Bayes: a list giving the MCMC specs for spline modeling, including burn-in length, number of posterior samples, thinning parameter and scaling factor for the size of random walk proposal
maxIntKnots = 6								          # Bayes: the maximum number of interior knots allowed for penalized splines
adaptive_r = TRUE                       #whether or not assign patients adaptively
start_num = 100        ## after this the number of total patients, we start on adaptive randomization;
## before this number, equal randomization is applied
block_size = 50        ## the frequency of updating on posterior distribution for the outcome
tuning = 1          ## the hyperparameter on probability of one arm is superior than the other
lambda = 0.001 
alpha = 0.05
default_marker_seq = FALSE #for plotting purpose
random_e = 0.1
futility_loop = 100
# futility_loop=30
p_fut = 0.1

fun_control = list(function(x) {0*x-2.95}, function(x) {0*x-2.95}, function(x) {2.1*x-2.95}, function(x) {0*x-2.95}, function(x) {0*x-2.95},
  function(x) {0*x-2.95},function(x) {0*x-2.95},function(x) {0*x-2.95},function(x) {0*x-2.95},function(x) {0*x-2.95},
function(x) {0*x -0.85},
#   function(x) {
#     x2 <- exp(30*(x-0.5))
#     y<- 2.1*x2/(1+x2)-2.95
#     return(y)
#   })

 # function(x) {0*x -0.85})
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
  # function(x) {0*x-2.95}
#  function(x) {
#    x2 <- exp(30*(x-0.5))
#    y<- 2.1*x2/(1+x2)-2.95
#    return(y)
#  }

 function(x) {
    x2 <- exp(30*(x-0.5))
    y<- 1.33*x2/(1+x2)-2.95 + 1.22*x
    return(y)
  }
)

FUN.ctl <-fun_control[[setting]]
FUN.exp <- fun_experiment[[setting]]

if(file.exists(paste0("/home1/yuetu/attempt_design/output/case",case,"/temp",seed,".RData"))){
  quit(save="no")}else{    
result <- adaptive_enrich_design( seed = seed,		            # Number of iterations and seed to reproduce results
    n = n, 							   	 	              # total number of patients expected to enroll (both arms combined)
    percent.exp = percent.exp,                        # randomization ratio to the experimental group
    marker.specs = marker.specs, 		# a list giving the distribution, range, reference point of continuous marker, and whether to include right censored subjects in calculation of interior knots location
    FUN.ctl = FUN.ctl,              # mapping function of the log hazard as a function of marker for the control group
    FUN.exp = FUN.exp,              # mapping function of the log hazard as a function of marker for the experimental group
    MCMC.specs = MCMC.specs,         		# Bayes: a list giving the MCMC specs for spline modeling, including burn-in length, number of posterior samples, thinning parameter and scaling factor for the size of random walk proposal
    MCMC.specs.fut = MCMC.specs.fut,
    maxIntKnots = maxIntKnots,								          # Bayes: the maximum number of interior knots allowed for penalized splines
    adaptive_r = adaptive_r,                       #whether or not assign patients adaptively
    start_num = start_num,        ## after this the number of total patients, we start on adaptive randomization;
    ## before this number, equal randomization is applied
    block_size = block_size,        ## the frequency of updating on posterior distribution for the outcome
    tuning = tuning,          ## the hyperparameter on probability of one arm is superior than the other
    lambda = lambda,
    alpha = alpha,
    random_e = random_e,
    futility_loop = futility_loop,
    p_fut = p_fut
    )
  
setwd(paste0("/home1/yuetu/attempt_design/output/case",setting))
save(result, file=paste0("temp",seed,".RData"))
}


