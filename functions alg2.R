#!/usr/local/bin/Rscript
pred_cal_fun <- function(dat, tau.fit, covars.outcome, outcome, outcome.type, 
                         iso.cal.fit.unpooled, iso.cal.fit.pooled, 
                         nk, cate.names){
  
  est.cate.test.mat <- list()
  iso.preds.unpooled <- list()
  
  for(s in 1:nk){
    # First we extract predictions across the nk learners
    temp <- wrapper_est_cate(dat, tau.fit[[s]], covars.outcome, 
                             outcome, outcome.type)
    temp <- do.call(cbind, temp)
    est.cate.test.mat[[s]] <- data.frame(temp)
    
    # Now we obtain the predictions from the nk unpooled calibrated learners
    temp <- lapply(1:ncol(est.cate.test.mat[[s]]), function(x) iso.cal.fit.unpooled[[s]][[x]]$iso.fit(est.cate.test.mat[[s]][, x]))
    temp <- do.call(cbind, temp)
    iso.preds.unpooled[[s]] <- data.frame(temp)
  }
  
  # Then we obtain median across predictions of nk calibrated learners
  iso.test.unpooled.preds <- matrix(0, ncol = length(cate.names), nrow = nrow(dat))
  for(j in 1:length(cate.names)){
    temp <- lapply(iso.preds.unpooled, function(x){x[, j]})
    temp <- do.call(cbind, temp)
    iso.test.unpooled.preds[, j] <- apply(temp, 1, function(x) median(x))  
  }
  
  # Compute pooled predictions
  iso.preds.pooled <- list()
  for(s in 1:nk){
    # compute predictions of calibrated learners 
    temp <- lapply(1:ncol(est.cate.test.mat[[s]]), function(x) iso.cal.fit.pooled[[x]]$iso.fit(est.cate.test.mat[[s]][, x]))
    temp <- do.call(cbind, temp)
    iso.preds.pooled[[s]] <- data.frame(temp)
  }
  
  # Now we take median across predictions of nk calibrated learners
  iso.test.pooled.preds <- matrix(0, ncol = length(cate.names), nrow = nrow(dat))
  for(j in 1:length(cate.names)){
    temp <- lapply(iso.preds.pooled, function(x){ x[, j]})
    temp <- do.call(cbind, temp)
    iso.test.pooled.preds[, j] <- apply(temp, 1, function(x) median(x))  
  }
  
  return(list(iso.test.unpooled.preds = iso.test.unpooled.preds, 
              iso.test.pooled.preds = iso.test.pooled.preds,
              est.cate.test.mat = est.cate.test.mat))
}

compute_quantiles <- function(test.mat.preds, tau.0.test, cate.names){
  
  quantiles.cate <- list()
  for(i in 1:length(cate.names)){
    quantiles.cate[[i]]  <- quantile(test.mat.preds[, i], 
                                     probs = c(0, 0.1, 0.9, 1))
    names(quantiles.cate[[i]])[1:4] <- c('q0','q10','q90','q100')
    quantiles.cate[[i]] <- data.frame(t(quantiles.cate[[i]]))
  }
  
  quantile.res.lower <- list()
  quantile.res.upper <- list()
  true.cate.lower <- list()
  true.cate.upper <- list()
  
  for(i in 1:ncol(test.mat.preds)){
    lower.idx.test.mat.preds <-  quantiles.cate[[i]]$q0 <= test.mat.preds[, i] & test.mat.preds[, i] <= quantiles.cate[[i]]$q10
    upper.idx.test.mat.preds <-  quantiles.cate[[i]]$q90 <= test.mat.preds[, i] & test.mat.preds[, i] <= quantiles.cate[[i]]$q100
    
    quantile.res.lower[[i]] <- mean(test.mat.preds[,i][lower.idx.test.mat.preds])
    quantile.res.upper[[i]] <- mean(test.mat.preds[,i][upper.idx.test.mat.preds])
    
    true.cate.lower[[i]] <- mean(tau.0.test[lower.idx.test.mat.preds])
    true.cate.upper[[i]] <- mean(tau.0.test[upper.idx.test.mat.preds])
  }
  return(list(quantile.res.lower = quantile.res.lower, 
              quantile.res.upper = quantile.res.upper,
              true.cate.lower = true.cate.lower, 
              true.cate.upper =  true.cate.upper))
}

wrapper_xgboost <- function(pred_col, tau0.gb, depth, test.col, tau.0.test){
  # fit xgboost
  xgboost.lrnr <- Lrnr_xgboost$new(max_depth = depth)
  temp.dat <- data.frame(tau = tau0.gb, preds = pred_col)
  task.cal <- make_sl3_Task(data = temp.dat, 
                            covariates = 'preds',
                            outcome = 'tau', 
                            outcome_type = "continuous")
  trained.model <- xgboost.lrnr$train(task.cal)
  
  # obtain predictions of E[tau0.gb|theta(tau)(W)] on test data 
  pred.dat <- data.frame(preds = test.col, tau = tau.0.test)
  pred.task <- make_sl3_Task(data = pred.dat, 
                             covariates = 'preds', 
                             outcome_type = "continuous")
  pred.gb <-  trained.model$predict(pred.task)
  
  # compute cal measure 
  res <- mean((pred.gb - test.col)^2)
  res.2 <- mean((tau.0.test - test.col) * (pred.gb - test.col) )
  return(list(res = res, res.2 = res.2))
}

xgboost_cal_measure <- function(tau.hat, 
                                tau.test.mat = NULL, 
                                tau.0.test,
                                sim_data_cal, 
                                n.gb = 100000, 
                                nk = NULL,
                                covars.outcome, outcome, 
                                outcome.type, 
                                iso.flag = FALSE, 
                                iso.cal.fit.unpooled = NULL, 
                                iso.cal.fit.pooled = NULL, 
                                iso.test.mat.unpooled = NULL, 
                                iso.test.mat.pooled = NULL, 
                                cate.names, 
                                uncal.flag = FALSE){

  # generate large data set to obtain predictions for fitting the gradient boosted trees
  dat.gb <- sim_data_cal(n.gb)
  tau0.gb <- dat.gb$Q1 - dat.gb$Q0
  dat.gb.cate <- cbind(A = dat.gb$A, dat.gb$W, Y = dat.gb$Y)
  
  # compute cal(tau) for calibrated pooled and unpooled estimator 
  if(iso.flag){

    iso.res.pooled <- list()
    iso.res.2.pooled <- list()
    iso.res.unpooled <- list()
    iso.res.2.unpooled <- list()
    
    preds_gbcal <- pred_cal_fun(dat.gb.cate, tau.hat, 
                                covars.outcome, outcome, outcome.type, 
                                iso.cal.fit.unpooled = iso.cal.fit.unpooled, 
                                iso.cal.fit.pooled = iso.cal.fit.pooled, nk, cate.names)
    
    iso.gb.preds.unpooled <- preds_gbcal$iso.test.unpooled.preds
    iso.gb.preds.pooled <- preds_gbcal$iso.test.pooled.preds
    
    #iso.gb.preds <- lapply(1:length(iso.cal.fit), function(x) iso.cal.fit[[x]]$iso.fit(tau.preds[, x]))
    #iso.gb.mat <- do.call(cbind, iso.gb.preds)
    
    for(i in 1:length(cate.names)){
      # obtain maximum depth
      #cal.params.temp <- cal.params %>% filter(method == cate.names[i])
      #mdepth <- as.numeric(cal.params.temp$cv.depth.xgboost)
      
      temp.unpooled <- wrapper_xgboost(iso.gb.preds.unpooled[,i], tau0.gb, 
                                       depth = 5, iso.test.mat.unpooled[,i], tau.0.test)
      
      iso.res.unpooled[[i]] <- temp.unpooled$res
      iso.res.2.unpooled[[i]] <- temp.unpooled$res.2
      
      temp.pooled <- wrapper_xgboost(iso.gb.preds.pooled[,i], tau0.gb, 
                                     depth = 5, iso.test.mat.pooled[,i], tau.0.test)
      
      iso.res.pooled[[i]] <- temp.pooled$res
      iso.res.2.pooled[[i]] <- temp.pooled$res.2
      
    }

  return(list(iso.res.pooled = iso.res.pooled, 
              iso.res.2.pooled = iso.res.2.pooled,
              iso.res.unpooled = iso.res.unpooled, 
              iso.res.2.unpooled = iso.res.2.unpooled))

  }
  if(uncal.flag){ 
    
    tau.res <- list()
    tau.res.2 <- list() 
    
    # predictions from uncalibrated CATE estimators 
    est.tau.cate <- wrapper_est_cate(dat.gb.cate, tau.hat, 
                                     covars.outcome, 
                                     outcome, outcome.type)
    
    # transform predictions to run xgboost regression 
    tau.gb.preds <- data.frame(est.tau.cate$pred.cate) 
    
    # compute cal(tau) for uncalibrated CATE estimators
    for(i in 1:length(cate.names)){
      
      # obtain maximum depth for cate estimator 
      # cal.params.temp <- cal.params %>% filter(method == cate.names[i])
      # mdepth <- as.numeric(cal.params.temp$cv.depth.xgboost.cate.uncal)
      
      temp <- wrapper_xgboost(tau.gb.preds[, i], tau0.gb, 
                              depth = 5, tau.test.mat[, i], tau.0.test)
      
      tau.res[[i]] <- temp$res
      tau.res.2[[i]] <- temp$res.2
    }
    return(list(tau.res = tau.res, tau.res.2 = tau.res.2))
  }
}


isotonic_reg <- function(pseudo.out, cate){
  iso.reg <- isoreg(cate, pseudo.out)
  iso.fit <- as.stepfun(iso.reg)
  return(list(iso.reg = iso.reg, iso.fit = iso.fit))
}

calcate_pseud_out <- function(fit.Q1, fit.Q0,  pa.1, A, Y){
  z <- (A - pa.1)/ ((pa.1)*(1 - pa.1)) * (Y - A*(fit.Q1) - (1-A) * fit.Q0) + (fit.Q1 - fit.Q0)
}

wrapper_est_cate <- function(dat, tau.fit, covars.outcome, outcome, outcome.type){
  # change treatment values to generate predictions with A = 1, and A = 0
  dat.1 <- dat %>% mutate(A=1)
  dat.0 <- dat %>% mutate(A=0)
  
  # make the task
  task.1 <- make_sl3_Task(data = dat.1, covariates = covars.outcome,
                          outcome = outcome, outcome_type = outcome.type)
  task.0 <- make_sl3_Task(data = dat.0, covariates = covars.outcome,
                          outcome = outcome, outcome_type = outcome.type)
  
  # get fitted values
  pred.cate <- tau.fit$predict(task.1) - tau.fit$predict(task.0)
  return(list(pred.cate = pred.cate))
}

pseudo_est <- function(dat,covars.outcome, covars.trt, outcome.type = 'binomial',
                       learners.outcome = list(Lrnr_glm$new(), Lrnr_gam$new(), 
                                               Lrnr_xgboost$new(max_depth =4), 
                                               Lrnr_glmnet$new()),
                       learners.trt = list(Lrnr_glm$new(), Lrnr_gam$new(), 
                                           Lrnr_xgboost$new(max_depth =4), 
                                           Lrnr_glmnet$new()), folds.dat, thresh = 0.01){
  
  Y <- dat$Y
  A <- dat$A
  
  # Create new task to fit model for E[A|W]
  task.A <- sl3_Task$new(dat, covariates = covars.trt, outcome = "A", folds = folds.dat)
  
  # Create new task to fit model for E[Y| W, A]
  # We specify same folds as in previous model
  task.Y <- sl3_Task$new(dat, covariates = covars.outcome, outcome = "Y", folds = folds.dat)
  
  # list of learners for both tasks
  lrnr.list.outcome <- learners.outcome
  lrnr.list.trt <- learners.trt
  
  stack.outcome <- Stack$new(lrnr.list.outcome)
  lrnr.cv.outcome <- Lrnr_cv$new(stack.outcome)
  
  stack.trt <- Stack$new(lrnr.list.trt)
  lrnr.cv.trt <- Lrnr_cv$new(stack.trt)
  
  lrnr.crossfit.sl.outcome <- make_learner(Pipeline, lrnr.cv.outcome, 
                                           Lrnr_cv_selector$new(loss_squared_error))
  
  lrnr.crossfit.sl.trt <- make_learner(Pipeline, lrnr.cv.trt, 
                                       Lrnr_cv_selector$new(loss_squared_error))
  
  lrnr.A <- lrnr.crossfit.sl.trt$train(task.A)
  lrnr.Y <- lrnr.crossfit.sl.outcome$train(task.Y)
  
  # obtain predictions for E[Y|A = 1, W] and E[Y|A = 0, W]
  dat.1 <- dat %>% mutate(A=1)
  dat.0 <- dat %>% mutate(A=0)
  
  task.1 <- make_sl3_Task(data = dat.1, covariates = covars.outcome,
                          outcome = "Y", outcome_type = outcome.type, 
                          folds = folds.dat)
  
  task.0 <- make_sl3_Task(data = dat.0, covariates = covars.outcome,
                          outcome = "Y", outcome_type = outcome.type, 
                          folds = folds.dat)
  
  Q1 <- lrnr.Y$predict(task.1) 
  Q0 <- lrnr.Y$predict(task.0)
  
  task.A <- sl3_Task$new(dat, covariates = covars.trt, outcome = "A", folds = folds.dat)
  pred.trt <- lrnr.A$predict(task.A)
  
  # bound estimated treatment probabilities with threshold
  pred.trt <- ifelse(pred.trt < thresh, thresh, pred.trt)
  pred.trt <- ifelse(pred.trt > 1-thresh, 1-thresh, pred.trt)
  
  pseudo.out <- calcate_pseud_out(Q1, Q0, pred.trt, A, Y)
  
  return(list(pseudo.out = pseudo.out))
}


# this function does a simulation iteration using algorithm 2 from paper
do.one.sim <- function(n.all, nk = 10, n.test, n.gb = 10000,
                       sim_data, sim_data_cal, Qbar, gbar,
                       outcome = 'Y', covars.outcome, covars.trt,
                       cate.estimators = list(xgboost = lrnr.xgboost, gam = lrnr.gam, 
                                           ranger = lrnr.ranger, glm = lrnr.glm, 
                                           earth = lrnr.earth, glmnet = lrnr.glmnet),
                       cate.names = c('xboost', 'gam', 'ranger', 'glm', 'earth', 'glmnet'),
                       learners.outcome = list(Lrnr_glm$new(), Lrnr_gam$new(), 
                                                  Lrnr_xgboost$new(max_depth =4), 
                                                  Lrnr_glmnet$new()),
                       learners.trt = list(Lrnr_glm$new(), Lrnr_gam$new(), 
                                           Lrnr_xgboost$new(max_depth =4), 
                                           Lrnr_glmnet$new()),
                       outcome.type = 'binomial',
                       thresh = 0.01){

  # dat.all is used for fitting calibrators and uncalibrated CATE estimator 
  dat.all<- sim_data(n.all)
  dat.all <- cbind(A = dat.all$A, dat.all$W, Y = dat.all$Y, 
                   Q1 = dat.all$Q1, Q0 = dat.all$Q0)
  
  # define folds for cross validation 
  # nk = number of splits for algorithm 2 (10 is default)
  folds.dat <-  origami::folds_vfold(n.all, V = nk)
   
  ## use cross fitting with dat.all to get fitted pseudo-outcomes
  cross.fit.out <- pseudo_est(dat.all, covars.outcome, covars.trt,
                              outcome.type, learners.outcome,
                              learners.trt, folds.dat, thresh)
  
  pseudo.outcomes <- cross.fit.out$pseudo.out
  
  # define list where we store the k-th fold CATE learners
  tau.fit <- list()
  init.tau <- Stack$new(cate.estimators)
  
  # data where we store predictions for each fold
  pred.cate.init.tau <- list()
  all.preds.cate <- matrix(0, nrow = n.all, ncol = length(cate.estimators))
  
  # Lines 2-6 in pseudo-code
  for(s in 1:nk){  
    # data to fit CATE learner at each fold
    dat.s <- dat.all[folds.dat[[s]]$training_set, ]
    # data to obtain predictions outside the sample 
    dat.pred <- dat.all[folds.dat[[s]]$validation_set, ]
    # fit CATE learners 
    task.s <- make_sl3_Task(data = dat.s, covariates = covars.outcome,
                            outcome = outcome, outcome_type = outcome.type)  
    tau.fit[[s]] <- init.tau$train(task.s)
    # Obtain predicted CATE from tau learners at fold s  
    est.cate <- wrapper_est_cate(dat.pred, tau.fit[[s]], covars.outcome, 
                                     outcome, outcome.type)
    pred.cate.init.tau[[s]] <- est.cate$pred.cate
    # Store all predictions for pooled calibrator
    all.preds.cate[folds.dat[[s]]$validation_set, ] <- as.matrix(pred.cate.init.tau[[s]])
  }

  # Define unpooled calibrators
  iso.cal.fit.unpooled <- list()
  # Line 8 in pseudocode for unpooled calibrator
  for(s in 1:nk){
    pseudo.outcomes.s <- pseudo.outcomes[folds.dat[[s]]$validation_set]
    iso.cal.fit.unpooled[[s]] <- apply(pred.cate.init.tau[[s]], 2, function(x){isotonic_reg(pseudo.outcomes.s, x)})
    # a prediction from unpooled is median across s of iso.cal.fit.unpooled[[s]](tau.fit[[s]](newobs))
  }
  
  # Construct pooled calibrators
  iso.cal.fit.pooled <- apply(all.preds.cate, 2, function(x){isotonic_reg(pseudo.outcomes, x)})
  # a prediction from iso.cal.fit.pooled is median across s of iso.cal.fit.pooled[[index learner]]]$iso.fit(tau.fit[[s]][[index learner]](newobs))
  
  # Fit learners that are uncalibrated with all available data
  task.all <- make_sl3_Task(data = dat.all, covariates = covars.outcome,
                            outcome = outcome, outcome_type = outcome.type)
  tau.uncal.fit <- init.tau$train(task.all)
  
  ################
  ### Test methods
  ################
  
  ## Generate test data 
  dat.test.all <- sim_data_cal(n.test)
  dat.test.Q1 <- dat.test.all$Q1
  dat.test.Q0 <- dat.test.all$Q0
  tau.0.test <- dat.test.Q1 - dat.test.Q0
  dat.test <- cbind(A = dat.test.all$A, dat.test.all$W, Y= dat.test.all$Y)
  
  #####################################################
  ## Predictions on test data from uncalibrated learner
  #####################################################
  est.uncal.cate.test <- wrapper_est_cate(dat.test, tau.uncal.fit, covars.outcome, 
                                          outcome, outcome.type)
  uncal.tau.test.preds <- data.frame(est.uncal.cate.test$pred.cate) 
  colnames(uncal.tau.test.preds) <- cate.names
  
  ####################################################
  ## Predictions on test data from calibrated learners
  ####################################################

  pred_cal <- pred_cal_fun(dat.test, tau.fit, 
                       covars.outcome, outcome, outcome.type, 
                       iso.cal.fit.unpooled, iso.cal.fit.pooled, nk, cate.names)
  
  ## Store predictions on test data from unpooled learner
  iso.test.unpooled.preds <- pred_cal$iso.test.unpooled.preds 
  colnames(iso.test.unpooled.preds) <- cate.names
  
  ## Store predictions on test data from pooled learner
  iso.test.pooled.preds <- pred_cal$iso.test.pooled.preds 
  colnames(iso.test.pooled.preds) <- cate.names
  
  ## Store predictions from nk CATE learners
  est.cate.test.mat <- pred_cal$est.cate.test.mat
  
  ###############################
  ## Compute performance measures
  #################################
  
  ########
  ## MSE
  #######
  
  ## Compute MSE for uncalibrated learner 
  uncal.tau.mse <- apply(uncal.tau.test.preds, 2, function(x) mean((dat.test.Q1 - dat.test.Q0 - x)^2))
 
  ## Compute MSE for unpooled learner 
  unpooled.mse <- apply(iso.test.unpooled.preds, 2, function(x) mean((dat.test.Q1 - dat.test.Q0 - x)^2))

  ## Compute MSE for pooled learner 
  pooled.mse <- apply(iso.test.pooled.preds, 2, function(x) mean((dat.test.Q1 - dat.test.Q0 - x)^2))

  ########################################
  ## Calibration measure and quantiles
  ########################################
  
  ## Get calibration measure and compute quantiles for uncalibrated estimator
  cal.measure.uncal <- xgboost_cal_measure(tau.hat = tau.uncal.fit, 
                                       tau.test.mat = uncal.tau.test.preds,
                                       tau.0.test = tau.0.test, 
                                       sim_data_cal = sim_data_cal, 
                                       n.gb = n.gb, 
                                       covars.outcome = covars.outcome, 
                                       outcome = outcome,
                                       outcome.type = outcome.type, 
                                       cate.names = cate.names, 
                                       uncal.flag = TRUE)
  quant_uncal <- compute_quantiles(uncal.tau.test.preds, tau.0.test, cate.names)

  ## Get calibration measure and quantiles for calibrated learners
  cal.measure.isocal <- xgboost_cal_measure(tau.hat = tau.fit, 
                                     tau.0.test = tau.0.test, 
                                     sim_data_cal = sim_data_cal, 
                                     n.gb = n.gb, 
                                     nk = nk,
                                     covars.outcome = covars.outcome, 
                                     outcome = outcome, 
                                     outcome.type = outcome.type, 
                                     iso.flag = TRUE, 
                                     iso.cal.fit.unpooled = iso.cal.fit.unpooled, 
                                     iso.cal.fit.pooled = iso.cal.fit.pooled, 
                                     iso.test.mat.unpooled = iso.test.unpooled.preds, 
                                     iso.test.mat.pooled = iso.test.pooled.preds, 
                                     cate.names = cate.names, 
                                     uncal.flag = FALSE)
  
  # Compute quantiles
  quant_unpooled <- compute_quantiles(iso.test.unpooled.preds, tau.0.test, cate.names)
  quant_pooled <- compute_quantiles(iso.test.pooled.preds, tau.0.test, cate.names)
 
  # Store calibration measure
  iso.cal <- do.call(cbind, cal.measure.isocal$iso.res.pooled)
  iso.cal.2 <- do.call(cbind, cal.measure.isocal$iso.res.2.pooled)
  iso.cal.unpooled <- do.call(cbind, cal.measure.isocal$iso.res.unpooled)
  iso.cal.2.unpooled <- do.call(cbind, cal.measure.isocal$iso.res.2.unpooled) 
  uncal.tau.cal <-  do.call(cbind, cal.measure.uncal$tau.res)
  uncal.tau.cal.2 <-  do.call(cbind, cal.measure.uncal$tau.res.2)
  
  # Store quantiles
  quantile.upper.cal <- do.call(cbind, quant_pooled$quantile.res.upper)
  quantile.lower.cal <- do.call(cbind, quant_pooled$quantile.res.lower)
  quantile.upper.tau0.cal <- do.call(cbind, quant_pooled$true.cate.upper)
  quantile.lower.tau0.cal <- do.call(cbind, quant_pooled$true.cate.lower)
  
  quantile.upper.cal.unpooled <- do.call(cbind, quant_unpooled$quantile.res.upper)
  quantile.lower.cal.unpooled <- do.call(cbind, quant_unpooled$quantile.res.lower)
  quantile.upper.tau0.cal.unpooled <- do.call(cbind, quant_unpooled$true.cate.upper)
  quantile.lower.tau0.cal.unpooled <- do.call(cbind, quant_unpooled$true.cate.lower)
  
  quantile.upper.uncal <- do.call(cbind, quant_uncal$quantile.res.upper)
  quantile.lower.uncal <- do.call(cbind, quant_uncal$quantile.res.lower)
  quantile.upper.tau0.uncal <- do.call(cbind, quant_uncal$true.cate.upper)
  quantile.lower.tau0.uncal <- do.call(cbind, quant_uncal$true.cate.lower)
  
  return(list(iso.cal = iso.cal, 
              iso.cal.2  = iso.cal.2, 
              iso.cal.unpooled = iso.cal.unpooled, 
              iso.cal.2.unpooled = iso.cal.2.unpooled, 
              uncal.tau.cal = uncal.tau.cal,
              uncal.tau.cal.2 = uncal.tau.cal.2,
              quantile.upper.cal = quantile.upper.cal,
              quantile.lower.cal = quantile.lower.cal,
              quantile.upper.tau0.cal = quantile.upper.tau0.cal,
              quantile.lower.tau0.cal = quantile.lower.tau0.cal,
              quantile.upper.cal.unpooled = quantile.upper.cal.unpooled, 
              quantile.lower.cal.unpooled = quantile.lower.cal.unpooled,
              quantile.upper.tau0.cal.unpooled = quantile.upper.tau0.cal.unpooled,
              quantile.lower.tau0.cal.unpooled = quantile.lower.tau0.cal.unpooled,
              quantile.upper.uncal =  quantile.upper.uncal, 
              quantile.lower.uncal= quantile.lower.uncal, 
              quantile.upper.tau0.uncal = quantile.upper.tau0.uncal,
              quantile.lower.tau0.uncal = quantile.lower.tau0.uncal,
              uncal.tau.mse = uncal.tau.mse,
              unpooled.mse = unpooled.mse,
              pooled.mse = pooled.mse)) }

do.many.sim <- function(nsim, n.all = n.all, nk = 10, n.test = n.test, n.gb = n.gb,
                         sim_data = sim_data, sim_data_cal = sim_data_cal,
                         Qbar = Qbar, gbar = gbar, outcome = 'Y', 
                         covars.outcome = covars.outcome, covars.trt = covars.trt,
                         cate.estimators = cate.estimators, cate.names = cate.names,
                         learners.outcome = learners.outcome, learners.trt = learners.trt ,
                         outcome.type =  outcome.type, thresh = 0.01){
   
   all.iso.cal <- matrix(-1, nrow = nsim, ncol = length(cate.names))
   all.iso.cal.2 <- matrix(-1, nrow = nsim, ncol = length(cate.names))  
   all.iso.cal.unpooled <- matrix(-1, nrow = nsim, ncol = length(cate.names)) 
   all.iso.cal.2.unpooled <- matrix(-1, nrow = nsim, ncol = length(cate.names))
   all.uncal.tau.cal <- matrix(-1, nrow = nsim, ncol = length(cate.names))
   all.uncal.tau.cal.2 <- matrix(-1, nrow = nsim, ncol = length(cate.names))
   all.quantile.upper.cal <- matrix(-1, nrow = nsim, ncol = length(cate.names))
   all.quantile.lower.cal <- matrix(-1, nrow = nsim, ncol = length(cate.names))
   all.quantile.upper.tau0.cal <- matrix(-1, nrow = nsim, ncol = length(cate.names))
   all.quantile.lower.tau0.cal <- matrix(-1, nrow = nsim, ncol = length(cate.names))
   all.quantile.upper.cal.unpooled <- matrix(-1, nrow = nsim, ncol = length(cate.names))
   all.quantile.lower.cal.unpooled <- matrix(-1, nrow = nsim, ncol = length(cate.names))
   all.quantile.upper.tau0.cal.unpooled <- matrix(-1, nrow = nsim, ncol = length(cate.names))
   all.quantile.lower.tau0.cal.unpooled <- matrix(-1, nrow = nsim, ncol = length(cate.names))
   all.quantile.upper.uncal <- matrix(-1, nrow = nsim, ncol = length(cate.names))
   all.quantile.lower.uncal <- matrix(-1, nrow = nsim, ncol = length(cate.names)) 
   all.quantile.upper.tau0.uncal <- matrix(-1, nrow = nsim, ncol = length(cate.names))
   all.quantile.lower.tau0.uncal <- matrix(-1, nrow = nsim, ncol = length(cate.names))
   all.uncal.tau.mse <- matrix(-1, nrow = nsim, ncol = length(cate.names))
   all.unpooled.mse <- matrix(-1, nrow = nsim, ncol = length(cate.names))
   all.pooled.mse <- matrix(-1, nrow = nsim, ncol = length(cate.names))
   
   
   for(i in 1:nsim){
     res.one <- do.one.sim(n.all = n.all, nk = 10, n.test = n.test, n.gb = n.gb,
                           sim_data = sim_data, sim_data_cal = sim_data_cal,
                           Qbar = Qbar, gbar = gbar, outcome = 'Y', 
                           covars.outcome = covars.outcome, covars.trt = covars.trt,
                           cate.estimators = cate.estimators, cate.names = cate.names,
                           learners.outcome = learners.outcome, learners.trt = learners.trt ,
                           outcome.type =  outcome.type, thresh = 0.01)
     
     all.iso.cal[i, ] <- res.one$iso.cal
     all.iso.cal.2[i, ] <- res.one$iso.cal.2 
     all.iso.cal.unpooled[i, ] <- res.one$iso.cal.unpooled
     all.iso.cal.2.unpooled[i, ] <- res.one$iso.cal.2.unpooled 
     all.uncal.tau.cal[i, ] <- res.one$uncal.tau.cal
     all.uncal.tau.cal.2[i, ] <- res.one$uncal.tau.cal.2
     
     all.quantile.upper.cal[i, ] <- res.one$quantile.upper.cal
     all.quantile.lower.cal[i, ] <- res.one$quantile.lower.cal
     all.quantile.upper.tau0.cal[i, ] <- res.one$quantile.upper.tau0.cal
     all.quantile.lower.tau0.cal[i, ] <- res.one$quantile.lower.tau0.cal
     
     all.quantile.upper.cal.unpooled[i, ] <- res.one$quantile.upper.cal.unpooled
     all.quantile.lower.cal.unpooled[i, ] <- res.one$quantile.lower.cal.unpooled
     all.quantile.upper.tau0.cal.unpooled[i, ] <- res.one$quantile.upper.tau0.cal.unpooled
     all.quantile.lower.tau0.cal.unpooled[i, ] <- res.one$quantile.lower.tau0.cal.unpooled
     
     all.quantile.upper.uncal[i, ] <- res.one$quantile.upper.uncal
     all.quantile.lower.uncal[i, ] <- res.one$quantile.lower.uncal
     all.quantile.upper.tau0.uncal[i, ] <- res.one$quantile.upper.tau0.uncal
     all.quantile.lower.tau0.uncal[i, ] <- res.one$quantile.lower.tau0.uncal
     
     all.uncal.tau.mse[i, ] <- res.one$uncal.tau.mse
     all.unpooled.mse[i, ] <- res.one$unpooled.mse
     all.pooled.mse[i, ] <- res.one$pooled.mse
   }
   
   names(all.iso.cal) <- cate.names
   names(all.iso.cal.2) <- cate.names
   names(all.iso.cal.unpooled) <- cate.names
   names(all.iso.cal.2.unpooled) <- cate.names
   names(all.uncal.tau.cal) <- cate.names
   names(all.uncal.tau.cal.2) <- cate.names
   names(all.quantile.upper.cal) <- cate.names
   names(all.quantile.lower.cal) <- cate.names
   names(all.quantile.upper.tau0.cal) <- cate.names
   names(all.quantile.lower.tau0.cal) <- cate.names
   names(all.quantile.upper.cal.unpooled) <- cate.names
   names(all.quantile.lower.cal.unpooled) <- cate.names
   names(all.quantile.upper.tau0.cal.unpooled) <- cate.names
   names(all.quantile.lower.tau0.cal.unpooled) <- cate.names
   names(all.quantile.upper.uncal) <- cate.names
   names(all.quantile.lower.uncal) <- cate.names
   names(all.quantile.upper.tau0.uncal) <- cate.names
   names(all.quantile.lower.tau0.uncal) <- cate.names
   names(all.uncal.tau.mse) <- cate.names
   names(all.unpooled.mse) <- cate.names
   names(all.pooled.mse) <- cate.names
   
  return(list(all.iso.cal = all.iso.cal, 
              all.iso.cal.2  = all.iso.cal.2, 
              all.iso.cal.unpooled = all.iso.cal.unpooled, 
              all.iso.cal.2.unpooled = all.iso.cal.2.unpooled, 
              all.uncal.tau.cal = all.uncal.tau.cal,
              all.uncal.tau.cal.2 = all.uncal.tau.cal.2,
              all.quantile.upper.cal = all.quantile.upper.cal,
              all.quantile.lower.cal = all.quantile.lower.cal,
              all.quantile.upper.tau0.cal = all.quantile.upper.tau0.cal,
              all.quantile.lower.tau0.cal = all.quantile.lower.tau0.cal,
              all.quantile.upper.cal.unpooled = all.quantile.upper.cal.unpooled, 
              all.quantile.lower.cal.unpooled = all.quantile.lower.cal.unpooled,
              all.quantile.upper.tau0.cal.unpooled = all.quantile.upper.tau0.cal.unpooled,
              all.quantile.lower.tau0.cal.unpooled = all.quantile.lower.tau0.cal.unpooled,
              all.quantile.upper.uncal = all.quantile.upper.uncal, 
              all.quantile.lower.uncal = all.quantile.lower.uncal, 
              all.quantile.upper.tau0.uncal = all.quantile.upper.tau0.uncal,
              all.quantile.lower.tau0.uncal = all.quantile.lower.tau0.uncal,
              all.uncal.tau.mse = all.uncal.tau.mse,
              all.unpooled.mse = all.unpooled.mse,
              all.pooled.mse = all.pooled.mse)) 
   
}
 
wrapper_res <- function(all.pooled.res, all.unpooled.res, all.uncal.res, 
                        param, cate.names, n.all, nk){
  
  temp_fun <- function(res, cate.names, adj, param){
    temp.out <- sapply(1:ncol(res), function(x) mean(res[, x]))
    temp.out <- data.frame(temp.out)
    temp.out$adj <- adj
    temp.out$method <- cate.names
    names(temp.out)[1] <- param 
    return(temp.out)
  }
  
  res.final.pooled <- temp_fun(all.pooled.res, cate.names, 'pooled', param)
  res.final.unpooled <- temp_fun(all.unpooled.res, cate.names, 'unpooled', param)
  res.final.uncal <- temp_fun(all.uncal.res, cate.names, 'uncal', param)
  
  res.final <- rbind(res.final.pooled, res.final.unpooled, res.final.uncal)
  res.final$nall <- n.all
  res.final$nk <- nk
  return(res.final)
}

texttolatex <- function(form, outcome){
  form2 <- gsub('w\\$W', 'w_', form)
  form3 <- gsub('\\*', '', form2)
  form4 <- paste(outcome, form3)
}



