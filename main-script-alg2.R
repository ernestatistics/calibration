#!/usr/local/bin/Rscript

# Packages we need
if(!require(plyr)){install.packages("plyr")}
if(!require(dplyr)){install.packages("dplyr")}
if(!require(np)){install.packages("np")}
if(!require(remotes)){install.packages("remotes")}
if(!require(sl3)) remotes::install_github("tlverse/sl3")
if(!require(SuperLearner)){install.packages("SuperLearner")}

args <- commandArgs(TRUE)

# local pc flag
if(getwd() == "/Users/Ernesto/Dropbox/UW/Research/Dissertation/Project 3"){
  localpc <- TRUE
}else{
  localpc <- FALSE
}

## Parse the arguments 
if(length(args)==0){
  print("No arguments supplied.")
} else {
  for (i in 1:length(args)){
    eval(parse(text = args[i]))
  }
}

if(localpc){
  scen <- 3
  source('Code/functions alg2.R')
  source('Code/user-inputs.R')
  n.all <- 100
  nsim <- 2
  nk <- 10
  n.test <- 200
  n.gb <- 1000
  # meta parameters for cluster 
  # nall, ntest, nker
}else{

  # loads from shell file: 
  n.all <- nall
  n.test <- ntest
  # scen also loads
  
  source('user-inputs.R')
  source('functions alg2.R')
  nsim <- 500
  n.gb <- 100000
  nk <- 10
  
}

start_time <- Sys.time()
all_res <- do.many.sim(nsim, n.all = n.all, nk = nk, n.test = n.test, n.gb = n.gb,
                      sim_data = sim_data, sim_data_cal = sim_data_cal,
                      Qbar = Qbar, gbar = gbar, outcome = 'Y', 
                      covars.outcome = covars.outcome, covars.trt = covars.trt,
                      cate.estimators = cate.estimators, cate.names = cate.names,
                      learners.outcome = learners.outcome, learners.trt = learners.trt ,
                      outcome.type =  outcome.type, thresh = 0.01)
  
# wrapper functions to take the mean and store results
all.res.mse <- wrapper_res(all_res$all.pooled.mse, all_res$all.unpooled.mse,
                           all_res$all.uncal.tau.mse, 'mse', cate.names, n.all, nk)

all.res.cal <- wrapper_res(all_res$all.iso.cal, all_res$all.iso.cal.unpooled,
                           all_res$all.uncal.tau.cal, 'cal', cate.names, n.all, nk)

all.res.cal.2 <- wrapper_res(all_res$all.iso.cal.2, all_res$all.iso.cal.2.unpooled,
                             all_res$all.uncal.tau.cal.2, 'cal2', cate.names, n.all, nk)

all.res.quant.upper <- wrapper_res(all_res$all.quantile.upper.cal, 
                                   all_res$all.quantile.upper.cal.unpooled,
                                   all_res$all.quantile.upper.uncal, 'upperq', 
                                   cate.names, n.all, nk)

all.res.quant.lower <- wrapper_res(all_res$all.quantile.lower.cal, 
                                   all_res$all.quantile.lower.cal.unpooled,
                                   all_res$all.quantile.lower.uncal, 'lowerq', 
                                   cate.names, n.all, nk)

all.res.quant.upper.tau0 <- wrapper_res(all_res$all.quantile.upper.tau0.cal, 
                                        all_res$all.quantile.upper.tau0.cal.unpooled,
                                        all_res$all.quantile.upper.tau0.uncal, 'upperqtau0', 
                                        cate.names, n.all, nk)

all.res.quant.lower.tau0 <- wrapper_res(all_res$all.quantile.lower.tau0.cal, 
                                        all_res$all.quantile.lower.tau0.cal.unpooled,
                                        all_res$all.quantile.lower.tau0.uncal, 'lowerqtau0', 
                                        cate.names, n.all, nk)

end_time <- Sys.time()
total_time <- end_time - start_time
cat(total_time)

if(!localpc){
  save(all_res, all.res.mse, all.res.cal, all.res.cal.2, all.res.quant.upper,
       all.res.quant.lower, all.res.quant.upper.tau0, all.res.quant.lower.tau0,
       file = paste(paste(n.all, scen, 'alg2', sep = '-'),'.RData', sep = ''))
}else{
  save(all_res, all.res.mse, all.res.cal, all.res.cal.2, all.res.quant.upper,
       all.res.quant.lower, all.res.quant.upper.tau0, all.res.quant.lower.tau0,
       file = paste(paste(n.all, scen, 'alg2', sep = '-'),'.RData', sep = ''))
}

