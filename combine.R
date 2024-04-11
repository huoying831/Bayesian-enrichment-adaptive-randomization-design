library(dplyr)

args <- commandArgs(TRUE)
case <- args[1]

empty <- NULL
samplesize_vet <- rep(NA,1000)
status.matrix <- matrix(NA, nrow = 1000, ncol = 99) 
if_end.matrix <- matrix(NA, nrow = 1000, ncol = 99)
end_time.matrix <-matrix(NA, nrow = 1000, ncol = 99)
marker_n.matrix <-matrix(NA, nrow = 1000, ncol = 99)
marker_exp.matrix <-matrix(NA, nrow = 1000, ncol = 99)
marker_seq <- data.frame(x = round(seq(0.01, 0.99, length.out = 99),6))

for(seed in 1:1000){
  if(file.exists(paste0("/home1/yuetu/attempt_design/output/case",case,"/temp",seed,".RData"))){
    
    load(paste0("/home1/yuetu/attempt_design/output/case",case,"/temp",seed,".RData"))
    
   # if(seed == 45){
   #   next
   # }
    result$marker.all.dat <- result$marker.all.dat %>% distinct(marker, .keep_all = TRUE)
    


    samplesize_vet[seed] <- result$samplesize
    status.matrix[seed,] <- result$marker.all.dat$status
    if_end.matrix[seed,] <- result$marker.all.dat$if_end
    end_time.matrix[seed,] <- result$marker.all.dat$end_time
    
    temp_n <- result$res_data %>% group_by(x) %>% summarize(n = n()) %>% 
      left_join(marker_seq,., by = "x") %>% mutate_at(vars(n), ~replace(., is.na(.), 0))
    marker_n.matrix[seed,] <- temp_n$n
    
    temp_exp_n <- result$res_data %>% group_by(x,trt) %>% summarize(n = n(),.groups = "keep") %>%
      ungroup() %>%
      filter(trt == 1) %>% select(-trt) %>% 
      left_join(marker_seq,., by = "x") %>% 
      mutate(exp_n = replace(n, is.na(n), 0)) %>% select(-n)
    marker_exp.matrix[seed,] <- temp_exp_n$exp_n
    
    
    if(seed == 1){
      marker.value <-result$marker.all.dat$marker
    } 
    
    
  }else{
    empty = c(empty,seed)
    next
  }
}

save.image(file = paste0("/home1/yuetu/attempt_design/output/case",case,"/summary_new.RData"))
