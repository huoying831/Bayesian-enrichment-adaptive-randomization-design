# setwd("/Users/ty/USCHPC/freq_design/output")
setwd("/home1/yuetu/freq_design/output")

## only take res.final,marker.prevalence,samplesize

sim_result <- matrix(NA, nrow = 12*1000, ncol = 5)

for(setting in 1:12){
  for(seed in 1:1000){
    
    load(paste0("case",setting,"seed",seed,".RData"))
    sim_result[(seed+(setting-1)*1000),1] <- setting
    sim_result[(seed+(setting-1)*1000),2] <- seed
    sim_result[(seed+(setting-1)*1000),3] <- result$res.final
    sim_result[(seed+(setting-1)*1000),4] <- result$marker.prevalence
    sim_result[(seed+(setting-1)*1000),5] <- result$sampsize
    
  }
}

save(sim_result, file=paste0("sim_all.RData"))
