#!/usr/local/bin/Rscript
## Functions for data generating mechanism and defining learners
if(scen == 0){
  
  Qbar = function(a,w){plogis(.5 + 1.5*a + 2.2*w$W1 - 1.5*w$W2)}
  
  gbar = function(w){plogis(-0.25 -w$W1 + .5*w$W2)}
  
  sim_data = function(nn){
    W = data.frame(W1=runif(nn, min = -1, max = 1), 
                   W2=runif(nn, min = -1, max = 1))
    PS <- gbar(W)
    A = rbinom(nn, 1, PS)
    Q1 <- Qbar(rep(1, nn), W)
    Q0 <- Qbar(rep(0, nn), W)
    Y1 = rbinom(nn, 1, Q1)
    Y0 = rbinom(nn, 1, Q0)
    Y <- A * Y1 + (1-A) * Y0
    return(list(W=W, A=A, Y=Y, Y1=Y1, Y0=Y0, Q1=Q1, Q0=Q0, PS=PS))
  }
  
  sim_data_cal = sim_data  
  outcome = 'Y' 
  covars.outcome =  c('A', 'W1', 'W2')
  covars.trt = c('W1', 'W2')
  outcome.type = 'binomial'   
  
  # CATE Estimator
  lrnr.glm <- make_learner(Lrnr_glm)
  cate.estimators = list(glm = lrnr.glm)
  cate.names = names(cate.estimators)
  
  # For pseudo-outcomes:
    ## Treatment regression
  learners.trt = list(Lrnr_glm$new())
  
    ## Outcome regression
  learners.outcome = list(Lrnr_glm$new())
  
}

if(scen == 1){
  
  Qbar = function(a,w){plogis(.5 + 1.5*a + 2 * w$W1 - 1.5 * w$W2 + 
                                2.5 * w$W3 + 1.5 * w$W4)}
  
  gbar = function(w){plogis(-0.25 -w$W1 + .5*w$W2 - w$W3 + 0.5 * w$W4)}
  
}

if(scen == 2){
  Qbar = function(a,w){plogis(-1.5 + 1.5 * a + 2 * a * w$W1 - 1.5 * a * w$W2  
                              + 2.5 * w$W3 + 1.5 * w$W4)}
  
  gbar = function(w){plogis(-0.25 -w$W1 + .5*w$W2 - w$W3 + 0.5 * w$W4)}
}

if(scen == 3){
  Qbar = function(a,w){plogis(-1.5 + 1.5 * a + 3 * a * w$W1 - 2.5 * (1-a) * w$W2  
                              + 2.5 * a * w$W3 + 1.5 * (1-a) * w$W4)}
  
  gbar = function(w){plogis(-0.25 -w$W1 + .5*w$W2 - w$W3 + 0.5 * w$W4)}
}

if(scen == 4){
  #Qbar = function(a,w){plogis(-1.5 + 1.5 * a + 2 * a * tan(w$W1) - 1.5 * a * w$W2 
  #                            + 2.5 * w$W3 + 1.5 * tan(w$W4))}
  
  Qbar = function(a,w){plogis(-1.5 + 1.5 * a + 2 * a * abs(w$W1) * abs(w$W2) - 2.5 * (1-a) * abs(w$W2) * w$W3
                              + 2.5 * w$W3 + 2.5 * (1-a) * sqrt(abs(w$W4)) - 1.5 * a * I(w$W2 < .5) +
                                1.5 * (1-a) * I(w$W4 < 0))}
  
  gbar = function(w){plogis(-0.25 -w$W1 + .5*w$W2 - w$W3 + 0.5 * w$W4)}
  
}

if(scen == 5){
  Qbar = function(a,w){plogis(-0.5 + 2.5 * a + 2 * w$W1 + 1.5 * w$W2 + .5 * w$W3 - 
                               2 * w$W4 + 1.5 * w$W5 - 3 * w$W6 + 1 * w$W7 - .5 * w$W8 +
                               2 * w$W9  - 2 * w$W10 + .5 * w$W11 - 1.5  * w$W12 +
                               w$W13 - .5 * w$W14 + 2 * w$W15 - 1 * w$W16 + 
                               2 * w$W17 + w$W18 + .5 *w$W19 - w$W20)}
  
  gbar = function(w){plogis(.2 -.4 * w$W1 + 0.4 * w$W2 - .4 * w$W3 + .5 * w$W4 -
                            - .5 * w$W5 + 0.4 * w$W6 - .4 * w$W7 + .5 * w$W8 - 
                            - .5 * w$W9 + .4 * w$W10 - .4 * w$W11 + .5 * w$W12 +
                            - .5 * w$W13 + .4 * w$W14 - .4 * w$W15 + 0.5 * w$W16 + 
                            - .5 * w$W17 + 0.4 * w$W18 - 0.4 * w$W19 + 0.5 * w$W20)}
}

if(scen == 6){
  
  Qbar = function(a,w){plogis(-0.5 + 3.5 * a + 2 * a * w$W1 + 5.5 * (1-a) * w$W2 + 
                               1.5 * a * w$W3 + 3 * (1-a) * w$W4 + 2.5 * a * w$W5 - 
                               5 * (1-a) * w$W6 + 1 * a * w$W7 + 2.5 * (1-a) * w$W8 +
                               1 * a * w$W9  + 1.5 * (1-a) * w$W10 + 
                               1.5 * w$W11 - 2.5  * w$W12 + w$W13 -
                               1.5 * w$W14 + 3 * w$W15 - 2 * w$W16 + 
                               3 * w$W17 - w$W18 + 1.5 *w$W19 - 2 * w$W20)}
  
  gbar = function(w){plogis(.2 -.4 * w$W1 + 0.4 * w$W2 - .4 * w$W3 + .5 * w$W4 -
                              - .5 * w$W5 + 0.4 * w$W6 - .4 * w$W7 + .5 * w$W8 - 
                              - .5 * w$W9 + .4 * w$W10 - .4 * w$W11 + .5 * w$W12 +
                              - .5 * w$W13 + .4 * w$W14 - .4 * w$W15 + 0.5 * w$W16 + 
                              - .5 * w$W17 + 0.4 * w$W18 - 0.4 * w$W19 + 0.5 * w$W20)}
}

if(scen == 7){
  Qbar = function(a,w){plogis(-0.5 + 3.5 * a + 3 * a * w$W1 + 6.5 * (1-a) * w$W2 + 
                                1.5 * a * w$W3 + 4 * (1-a) * w$W4 + 2.5 * a * w$W5 - 
                                6 * (1-a) * w$W6 + 1 * a * w$W7 + 4.5 * (1-a) * w$W8 +
                                1 * a * w$W9  + 2.5 * (1-a) * w$W10 + 
                                1.5 * w$W11 - 2.5  * w$W12 + w$W13 -
                                1.5 * w$W14 + 3 * w$W15 - 2 * w$W16 + 
                                3 * w$W17 - w$W18 + 1.5 *w$W19 - 2 * w$W20)}
  
  gbar = function(w){plogis(.2 -.4 * w$W1 + 0.4 * w$W2 - .4 * w$W3 + .5 * w$W4 -
                              - .5 * w$W5 + 0.4 * w$W6 - .4 * w$W7 + .5 * w$W8 - 
                              - .5 * w$W9 + .4 * w$W10 - .4 * w$W11 + .5 * w$W12 +
                              - .5 * w$W13 + .4 * w$W14 - .4 * w$W15 + 0.5 * w$W16 + 
                              - .5 * w$W17 + 0.4 * w$W18 - 0.4 * w$W19 + 0.5 * w$W20)}
}

if(scen == 8){
  Qbar = function(a,w){plogis(+ 1.5 + 1.5 * a - 2.5 * (1-a) * abs(w$W1) - 3.5 * (1-a) * abs(w$W2) * abs(w$W3) - 
                                2.5 * (1-a) * abs(w$W3) + 2 * (1-a) - 2 * abs(w$W4) + 3.5 * (1-a) * I(w$W5 < 0.5) + 
                                3 * (1-a) * I(w$W6 >0.5) * abs(w$W7) + 1.5 * a * w$W7 + 3.5 * a * I(w$W8 > 0.5) +
                                2.5 * a * w$W9 + 2 * a * I(w$W10 < .2) + 
                                1.5 * a * abs(w$W11) - 2.5 * a * w$W12 + w$W13 -
                                1.5 * abs(w$W14) + 2.5 * a * I(w$W15 > .2) * abs(w$W14) - 
                                2 * a * abs(w$W16) * w$W17 + 
                                3 * w$W17 - w$W18 + 1.5 * w$W19 - w$W20)}
  
  gbar = function(w){plogis(.2 -.4 * w$W1 + 0.4 * w$W2 - .4 * w$W3 + .5 * w$W4 -
                              - .5 * w$W5 + 0.4 * w$W6 - .4 * w$W7 + .5 * w$W8 - 
                              - .5 * w$W9 + .4 * w$W10 - .4 * w$W11 + .5 * w$W12 +
                              - .5 * w$W13 + .4 * w$W14 - .4 * w$W15 + 0.5 * w$W16 + 
                              - .5 * w$W17 + 0.4 * w$W18 - 0.4 * w$W19 + 0.5 * w$W20)}
}

if(scen == 1 || scen == 2 || scen == 3 || scen == 4){
  sim_data = function(nn){
    W = data.frame(W1=runif(nn, min = -1, max = 1), 
                   W2=runif(nn, min = -1, max = 1), 
                   W3=runif(nn, min = -1, max = 1), 
                   W4=runif(nn, min = -1, max = 1))
    PS <- gbar(W)
    A = rbinom(nn, 1, PS)
    Q1 <- Qbar(rep(1, nn), W)
    Q0 <- Qbar(rep(0, nn), W)
    Y1 = rbinom(nn, 1, Q1)
    Y0 = rbinom(nn, 1, Q0)
    Y <- A * Y1 + (1-A) * Y0
    return(list(W=W,A=A,Y=Y,Y1=Y1,Y0=Y0,Q1=Q1,Q0=Q0,PS=PS))
}
  
  sim_data_cal = sim_data  
  outcome = 'Y' 
  covars.outcome =  c('A', 'W1', 'W2', 'W3', 'W4')
  covars.trt = c('W1', 'W2', 'W3', 'W4')
  outcome.type = 'binomial' 
}

if(scen == 5 || scen == 6 || scen == 7 || scen == 8){
  sim_data = function(nn){
    W = data.frame(W1=runif(nn, min = -1, max = 1), W2=runif(nn, min = -1, max = 1), W3=runif(nn, min = -1, max = 1), W4=runif(nn, min = -1, max = 1), W5=runif(nn, min = -1, max = 1), 
                   W6=runif(nn, min = -1, max = 1), W7=runif(nn, min = -1, max = 1), W8=runif(nn, min = -1, max = 1), W9=runif(nn, min = -1, max = 1), W10=runif(nn, min = -1, max = 1), 
                   W11=runif(nn, min = -1, max = 1), W12=runif(nn, min = -1, max = 1), W13=runif(nn, min = -1, max = 1), W14=runif(nn, min = -1, max = 1), W15=runif(nn, min = -1, max = 1), 
                   W16=runif(nn, min = -1, max = 1), W17=runif(nn, min = -1, max = 1), W18=runif(nn, min = -1, max = 1), W19=runif(nn, min = -1, max = 1), W20=runif(nn, min = -1, max = 1), 
                   W21=runif(nn, min = -1, max = 1), W22=runif(nn, min = -1, max = 1), W23=runif(nn, min = -1, max = 1), W24=runif(nn, min = -1, max = 1), W25=runif(nn, min = -1, max = 1), 
                   W26=runif(nn, min = -1, max = 1), W27=runif(nn, min = -1, max = 1), W28=runif(nn, min = -1, max = 1), W29=runif(nn, min = -1, max = 1), W30=runif(nn, min = -1, max = 1), 
                   W31=runif(nn, min = -1, max = 1), W32=runif(nn, min = -1, max = 1), W33=runif(nn, min = -1, max = 1), W34=runif(nn, min = -1, max = 1), W35=runif(nn, min = -1, max = 1), 
                   W36=runif(nn, min = -1, max = 1), W37=runif(nn, min = -1, max = 1), W38=runif(nn, min = -1, max = 1), W39=runif(nn, min = -1, max = 1), W40=runif(nn, min = -1, max = 1), 
                   W41=runif(nn, min = -1, max = 1), W42=runif(nn, min = -1, max = 1), W43=runif(nn, min = -1, max = 1), W44=runif(nn, min = -1, max = 1), W45=runif(nn, min = -1, max = 1), 
                   W46=runif(nn, min = -1, max = 1), W47=runif(nn, min = -1, max = 1), W48=runif(nn, min = -1, max = 1), W49=runif(nn, min = -1, max = 1), W50=runif(nn, min = -1, max = 1), 
                   W51=runif(nn, min = -1, max = 1), W52=runif(nn, min = -1, max = 1), W53=runif(nn, min = -1, max = 1), W54=runif(nn, min = -1, max = 1), W55=runif(nn, min = -1, max = 1), 
                   W56=runif(nn, min = -1, max = 1), W57=runif(nn, min = -1, max = 1), W58=runif(nn, min = -1, max = 1), W59=runif(nn, min = -1, max = 1), W60=runif(nn, min = -1, max = 1), 
                   W61=runif(nn, min = -1, max = 1), W62=runif(nn, min = -1, max = 1), W63=runif(nn, min = -1, max = 1), W64=runif(nn, min = -1, max = 1), W65=runif(nn, min = -1, max = 1), 
                   W66=runif(nn, min = -1, max = 1), W67=runif(nn, min = -1, max = 1), W68=runif(nn, min = -1, max = 1), W69=runif(nn, min = -1, max = 1), W70=runif(nn, min = -1, max = 1), 
                   W71=runif(nn, min = -1, max = 1), W72=runif(nn, min = -1, max = 1), W73=runif(nn, min = -1, max = 1), W74=runif(nn, min = -1, max = 1), W75=runif(nn, min = -1, max = 1), 
                   W76=runif(nn, min = -1, max = 1), W77=runif(nn, min = -1, max = 1), W78=runif(nn, min = -1, max = 1), W79=runif(nn, min = -1, max = 1), W80=runif(nn, min = -1, max = 1), 
                   W81=runif(nn, min = -1, max = 1), W82=runif(nn, min = -1, max = 1), W83=runif(nn, min = -1, max = 1), W84=runif(nn, min = -1, max = 1), W85=runif(nn, min = -1, max = 1), 
                   W86=runif(nn, min = -1, max = 1), W87=runif(nn, min = -1, max = 1), W88=runif(nn, min = -1, max = 1), W89=runif(nn, min = -1, max = 1), W90=runif(nn, min = -1, max = 1), 
                   W91=runif(nn, min = -1, max = 1), W92=runif(nn, min = -1, max = 1), W93=runif(nn, min = -1, max = 1), W94=runif(nn, min = -1, max = 1), W95=runif(nn, min = -1, max = 1), 
                   W96=runif(nn, min = -1, max = 1), W97=runif(nn, min = -1, max = 1), W98=runif(nn, min = -1, max = 1), W99=runif(nn, min = -1, max = 1), W100=runif(nn, min = -1, max = 1))
    PS <- gbar(W)
    A = rbinom(nn, 1, PS)
    Q1 <- Qbar(rep(1, nn), W)
    Q0 <- Qbar(rep(0, nn), W)
    Y1 = rbinom(nn, 1, Q1)
    Y0 = rbinom(nn, 1, Q0)
    Y <- A * Y1 + (1-A) * Y0
    return(list(W=W,A=A,Y=Y,Y1=Y1,Y0=Y0,Q1=Q1,Q0=Q0, PS=PS))
  }
  
  sim_data_cal = sim_data  
  outcome = 'Y' 
  covars.outcome = c('A','W1','W2','W3','W4','W5','W6','W7',
                     'W8','W9','W10','W11','W12','W13',
                     'W14','W15','W16','W17','W18','W19',
                     'W20','W21','W22','W23','W24','W25',
                     'W26','W27','W28','W29','W30','W31',
                     'W32','W33','W34','W35','W36','W37',
                     'W38','W39','W40','W41','W42','W43',
                     'W44','W45','W46','W47','W48','W49',
                     'W50','W51','W52','W53','W54','W55',
                     'W56','W57','W58','W59','W60','W61',
                     'W62','W63','W64','W65','W66','W67',
                     'W68','W69','W70','W71','W72','W73',
                     'W74','W75','W76','W77','W78','W79',
                     'W80','W81','W82','W83','W84','W85',
                     'W86','W87','W88','W89','W90','W91',
                     'W92','W93','W94','W95','W96','W97',
                     'W98','W99','W100')
  
  covars.trt = c('W1','W2','W3','W4','W5','W6','W7',
                 'W8','W9','W10','W11','W12','W13',
                 'W14','W15','W16','W17','W18','W19',
                 'W20','W21','W22','W23','W24','W25',
                 'W26','W27','W28','W29','W30','W31',
                 'W32','W33','W34','W35','W36','W37',
                 'W38','W39','W40','W41','W42','W43',
                 'W44','W45','W46','W47','W48','W49',
                 'W50','W51','W52','W53','W54','W55',
                 'W56','W57','W58','W59','W60','W61',
                 'W62','W63','W64','W65','W66','W67',
                 'W68','W69','W70','W71','W72','W73',
                 'W74','W75','W76','W77','W78','W79',
                 'W80','W81','W82','W83','W84','W85',
                 'W86','W87','W88','W89','W90','W91',
                 'W92','W93','W94','W95','W96','W97',
                 'W98','W99','W100')
  outcome.type = 'binomial'   
}

########################################
## Tau Learners and pseudo-outcomes #### 
#########################################

if(scen == 1 || scen == 2 || scen == 3 || scen == 4){
  
  # For pseudo-outcomes:
  ## Treatment regression
   learners.trt = list(Lrnr_glm$new(), Lrnr_gam$new(), 
                       Lrnr_xgboost$new(max_depth = 2), 
                       Lrnr_xgboost$new(max_depth = 4), 
                       Lrnr_xgboost$new(max_depth = 6), 
                       Lrnr_glmnet$new())
  
  #learners.trt = list(Lrnr_glm$new(), 
  #                    Lrnr_glmnet$new())
  
  ## Outcome regression
  learners.outcome = list(Lrnr_glm$new(), Lrnr_gam$new(), 
                          Lrnr_xgboost$new(max_depth = 2), 
                          Lrnr_xgboost$new(max_depth = 4), 
                          Lrnr_xgboost$new(max_depth = 6), 
                          Lrnr_glmnet$new())
  
  ## CATE learners
  lrnr.gam <- make_learner(Lrnr_gam)
  lrnr.ranger <- make_learner(Lrnr_ranger)
  lrnr.earth <- make_learner(Lrnr_earth)
  lrnr.glm <- make_learner(Lrnr_glm)
  lrnr.glmnet <- make_learner(Lrnr_glmnet)
  
  lrnr.xgboost2 <- make_learner(Lrnr_xgboost, eval_metric = 'logloss', 
                                max.depth = 2)
  lrnr.xgboost3 <- make_learner(Lrnr_xgboost, eval_metric = 'logloss', 
                                max.depth = 3)
  lrnr.xgboost5 <- make_learner(Lrnr_xgboost, eval_metric = 'logloss', 
                                max.depth = 5)
  lrnr.xgboost6 <- make_learner(Lrnr_xgboost, eval_metric = 'logloss', 
                                max.depth = 6)
  lrnr.xgboost8 <- make_learner(Lrnr_xgboost, eval_metric = 'logloss', 
                                max.depth = 8)
  
  cate.estimators = list(gam = lrnr.gam, 
                         ranger = lrnr.ranger,
                         earth = lrnr.earth, 
                         glm = lrnr.glm, 
                         glmnet = lrnr.glmnet, 
                         xgboost2 = lrnr.xgboost2, 
                         xgboost3 = lrnr.xgboost3, 
                         xgboost5 = lrnr.xgboost5, 
                         xgboost6 = lrnr.xgboost6,
                         xgboost8 = lrnr.xgboost8)
  
  cate.names = names(cate.estimators)
  #cate.names = c('gam', 'ranger', 'earth', 'glm', 'glmnet', 
                 #'xboost2', 'xboost3', 'xboost5', 'xboost6', 'xboost8')

}

if(scen == 5 || scen == 6 || scen == 7 || scen == 8){
  
  # For pseudo-outcomes:
  ## Treatment regression
  learners.trt = list(Lrnr_glmnet$new())
  
  ## Outcome regression
  learners.outcome = list(Lrnr_glmnet$new())
  
  ## CATE learners
  lrnr.glmnet <- make_learner(Lrnr_glmnet)
  lrnr.ranger <- make_learner(Lrnr_ranger)
  
  screen_cor <- Lrnr_pkg_SuperLearner_screener$new("screen.corP")
  lrnr.xgboost <-  Lrnr_xgboost$new()

  lrnr.glmnetxgb <- make_learner(Pipeline, screen_cor,lrnr.xgboost)
  
  cate.estimators = list(glmnet = lrnr.glmnet, 
                         glmnetxgb = lrnr.glmnetxgb,
                         ranger = lrnr.ranger)
  
  cate.names = names(cate.estimators)
  
}

# scenario 3 interactions
# scenario 4 nonlinearities 
# scenario 7 interactions
# scenario 8 nonlinearities

if(scen == 0 || scen == 3 || scen == 8 || scen == 9){
  flag.gld.std <- TRUE
}else{
  flag.gld.std <- FALSE
}

#####################
## Y is continuous ##
#####################

# scenario that is the same as 3 but with y continuous 
if(scen == 9){   
  Qbar = function(a,w){-1.5 + 1.5 * a + 3 * a * w$W1 - 2.5 * (1-a) * w$W2 + 2.5 * a * w$W3 + 1.5 * (1-a) * w$W4}
  gbar = function(w){plogis(-0.25 -w$W1 + .5*w$W2 - w$W3 + 0.5 * w$W4)}
}

# scenario that is the same as 4 but with y continuous  
if(scen == 10){   
  Qbar = function(a,w){-1.5 + 1.5 * a + 2 * a * abs(w$W1) * abs(w$W2) - 2.5 * (1-a) * abs(w$W2) * w$W3
                       + 2.5 * w$W3 + 2.5 * (1-a) * sqrt(abs(w$W4)) - 1.5 * a * I(w$W2 < .5) + 1.5 * (1-a) * I(w$W4 < 0)}
  
  gbar = function(w){plogis(-0.25 -w$W1 + .5*w$W2 - w$W3 + 0.5 * w$W4)}
  
}

if(scen == 9 || scen == 10){
  sim_data = function(nn){
    W = data.frame(W1=runif(nn, min = -1, max = 1), 
                   W2=runif(nn, min = -1, max = 1), 
                   W3=runif(nn, min = -1, max = 1), 
                   W4=runif(nn, min = -1, max = 1))
    PS <- gbar(W)
    A = rbinom(nn, 1, PS)
    Q1 <- Qbar(rep(1, nn), W)
    Q0 <- Qbar(rep(0, nn), W)
    Y1 = rnorm(nn, Q1, 1)
    Y0 = rnorm(nn, Q0, 1)
    Y <- A * Y1 + (1-A) * Y0
    return(list(W=W,A=A,Y=Y,Y1=Y1,Y0=Y0,Q1=Q1,Q0=Q0,PS=PS))
  }
  
  sim_data_cal = sim_data  
  outcome = 'Y' 
  covars.outcome =  c('A', 'W1', 'W2', 'W3', 'W4')
  covars.trt = c('W1', 'W2', 'W3', 'W4')
  outcome.type = 'continuous' 
}

if(scen == 9 || scen == 10){
  # For pseudo-outcomes:
  ## Treatment regression
  learners.trt = list(Lrnr_glm$new(), Lrnr_gam$new(), 
                      Lrnr_xgboost$new(max_depth = 2), 
                      Lrnr_xgboost$new(max_depth = 4), 
                      Lrnr_xgboost$new(max_depth = 6), 
                      Lrnr_glmnet$new())
  
  ## Outcome regression
  learners.outcome = list(Lrnr_glm$new(), Lrnr_gam$new(), 
                          Lrnr_xgboost$new(max_depth = 2), 
                          Lrnr_xgboost$new(max_depth = 4), 
                          Lrnr_xgboost$new(max_depth = 6), 
                          Lrnr_glmnet$new())
  
  ## CATE learners
  lrnr.gam <- make_learner(Lrnr_gam)
  lrnr.ranger <- make_learner(Lrnr_ranger)
  lrnr.earth <- make_learner(Lrnr_earth)
  lrnr.glm <- make_learner(Lrnr_glm)
  lrnr.glmnet <- make_learner(Lrnr_glmnet)
  
  lrnr.xgboost2 <- make_learner(Lrnr_xgboost, eval_metric = 'logloss', 
                                max.depth = 2)
  lrnr.xgboost3 <- make_learner(Lrnr_xgboost, eval_metric = 'logloss', 
                                max.depth = 3)
  lrnr.xgboost5 <- make_learner(Lrnr_xgboost, eval_metric = 'logloss', 
                                max.depth = 5)
  lrnr.xgboost6 <- make_learner(Lrnr_xgboost, eval_metric = 'logloss', 
                                max.depth = 6)
  lrnr.xgboost8 <- make_learner(Lrnr_xgboost, eval_metric = 'logloss', 
                                max.depth = 8)
  
  cate.estimators = list(gam = lrnr.gam, 
                         ranger = lrnr.ranger,
                         earth = lrnr.earth, 
                         glm = lrnr.glm, 
                         glmnet = lrnr.glmnet, 
                         xgboost2 = lrnr.xgboost2, 
                         xgboost3 = lrnr.xgboost3, 
                         xgboost5 = lrnr.xgboost5, 
                         xgboost6 = lrnr.xgboost6,
                         xgboost8 = lrnr.xgboost8)
  
  cate.names = names(cate.estimators)
  #cate.names = c('gam', 'ranger', 'earth', 'glm', 'glmnet', 
  #'xboost2', 'xboost3', 'xboost5', 'xboost6', 'xboost8')
}

## scenario that is the same as 7 but with y continuous 
if(scen == 11){
  Qbar = function(a,w){-0.5 + 3.5 * a + 3 * a * w$W1 + 6.5 * (1-a) * w$W2 + 
                       1.5 * a * w$W3 + 4 * (1-a) * w$W4 + 2.5 * a * w$W5 - 
                       6 * (1-a) * w$W6 + 1 * a * w$W7 + 4.5 * (1-a) * w$W8 +
                       1 * a * w$W9  + 2.5 * (1-a) * w$W10 + 
                       1.5 * w$W11 - 2.5  * w$W12 + w$W13 -
                       1.5 * w$W14 + 3 * w$W15 - 2 * w$W16 + 
                       3 * w$W17 - w$W18 + 1.5 *w$W19 - 2 * w$W20}
  
  gbar = function(w){plogis(.2 -.4 * w$W1 + 0.4 * w$W2 - .4 * w$W3 + .5 * w$W4 -
                              - .5 * w$W5 + 0.4 * w$W6 - .4 * w$W7 + .5 * w$W8 - 
                              - .5 * w$W9 + .4 * w$W10 - .4 * w$W11 + .5 * w$W12 +
                              - .5 * w$W13 + .4 * w$W14 - .4 * w$W15 + 0.5 * w$W16 + 
                              - .5 * w$W17 + 0.4 * w$W18 - 0.4 * w$W19 + 0.5 * w$W20)}
}

## scenario that is the same as 8 but with y continuous 
if(scen == 12){
  Qbar = function(a,w){+ 1.5 + 1.5 * a - 2.5 * (1-a) * abs(w$W1) - 3.5 * (1-a) * abs(w$W2) * abs(w$W3) - 
                        2.5 * (1-a) * abs(w$W3) + 2 * (1-a) - 2 * abs(w$W4) + 3.5 * (1-a) * I(w$W5 < 0.5) + 
                        3 * (1-a) * I(w$W6 >0.5) * abs(w$W7) + 1.5 * a * w$W7 + 3.5 * a * I(w$W8 > 0.5) +
                        2.5 * a * w$W9 + 2 * a * I(w$W10 < .2) + 
                        1.5 * a * abs(w$W11) - 2.5 * a * w$W12 + w$W13 -
                        1.5 * abs(w$W14) + 2.5 * a * I(w$W15 > .2) * abs(w$W14) - 
                        2 * a * abs(w$W16) * w$W17 + 
                        3 * w$W17 - w$W18 + 1.5 * w$W19 - w$W20}
  
  gbar = function(w){plogis(.2 -.4 * w$W1 + 0.4 * w$W2 - .4 * w$W3 + .5 * w$W4 -
                              - .5 * w$W5 + 0.4 * w$W6 - .4 * w$W7 + .5 * w$W8 - 
                              - .5 * w$W9 + .4 * w$W10 - .4 * w$W11 + .5 * w$W12 +
                              - .5 * w$W13 + .4 * w$W14 - .4 * w$W15 + 0.5 * w$W16 + 
                              - .5 * w$W17 + 0.4 * w$W18 - 0.4 * w$W19 + 0.5 * w$W20)}
}

if(scen == 11 || scen == 12){
  # For pseudo-outcomes:
  ## Treatment regression
  learners.trt = list(Lrnr_glmnet$new(), Lrnr_glmnet$new())
  
  ## Outcome regression
  learners.outcome = list(Lrnr_glmnet$new(), Lrnr_glmnet$new())
  
  ## CATE learners
  lrnr.glmnet <- make_learner(Lrnr_glmnet)
  lrnr.ranger <- make_learner(Lrnr_ranger)
  
  screen_cor <- Lrnr_pkg_SuperLearner_screener$new("screen.corP")
  lrnr.xgboost <-  Lrnr_xgboost$new()
  
  lrnr.glmnetxgb <- make_learner(Pipeline, screen_cor,lrnr.xgboost)
  
  cate.estimators = list(glmnet = lrnr.glmnet, 
                         glmnetxgb = lrnr.glmnetxgb, 
                         ranger = lrnr.ranger)
  
  cate.names = names(cate.estimators)
}

if(scen == 11 || scen == 12){
  sim_data = function(nn){
    W = data.frame(W1=runif(nn, min = -1, max = 1), W2=runif(nn, min = -1, max = 1), W3=runif(nn, min = -1, max = 1), W4=runif(nn, min = -1, max = 1), W5=runif(nn, min = -1, max = 1), 
                   W6=runif(nn, min = -1, max = 1), W7=runif(nn, min = -1, max = 1), W8=runif(nn, min = -1, max = 1), W9=runif(nn, min = -1, max = 1), W10=runif(nn, min = -1, max = 1), 
                   W11=runif(nn, min = -1, max = 1), W12=runif(nn, min = -1, max = 1), W13=runif(nn, min = -1, max = 1), W14=runif(nn, min = -1, max = 1), W15=runif(nn, min = -1, max = 1), 
                   W16=runif(nn, min = -1, max = 1), W17=runif(nn, min = -1, max = 1), W18=runif(nn, min = -1, max = 1), W19=runif(nn, min = -1, max = 1), W20=runif(nn, min = -1, max = 1), 
                   W21=runif(nn, min = -1, max = 1), W22=runif(nn, min = -1, max = 1), W23=runif(nn, min = -1, max = 1), W24=runif(nn, min = -1, max = 1), W25=runif(nn, min = -1, max = 1), 
                   W26=runif(nn, min = -1, max = 1), W27=runif(nn, min = -1, max = 1), W28=runif(nn, min = -1, max = 1), W29=runif(nn, min = -1, max = 1), W30=runif(nn, min = -1, max = 1), 
                   W31=runif(nn, min = -1, max = 1), W32=runif(nn, min = -1, max = 1), W33=runif(nn, min = -1, max = 1), W34=runif(nn, min = -1, max = 1), W35=runif(nn, min = -1, max = 1), 
                   W36=runif(nn, min = -1, max = 1), W37=runif(nn, min = -1, max = 1), W38=runif(nn, min = -1, max = 1), W39=runif(nn, min = -1, max = 1), W40=runif(nn, min = -1, max = 1), 
                   W41=runif(nn, min = -1, max = 1), W42=runif(nn, min = -1, max = 1), W43=runif(nn, min = -1, max = 1), W44=runif(nn, min = -1, max = 1), W45=runif(nn, min = -1, max = 1), 
                   W46=runif(nn, min = -1, max = 1), W47=runif(nn, min = -1, max = 1), W48=runif(nn, min = -1, max = 1), W49=runif(nn, min = -1, max = 1), W50=runif(nn, min = -1, max = 1), 
                   W51=runif(nn, min = -1, max = 1), W52=runif(nn, min = -1, max = 1), W53=runif(nn, min = -1, max = 1), W54=runif(nn, min = -1, max = 1), W55=runif(nn, min = -1, max = 1), 
                   W56=runif(nn, min = -1, max = 1), W57=runif(nn, min = -1, max = 1), W58=runif(nn, min = -1, max = 1), W59=runif(nn, min = -1, max = 1), W60=runif(nn, min = -1, max = 1), 
                   W61=runif(nn, min = -1, max = 1), W62=runif(nn, min = -1, max = 1), W63=runif(nn, min = -1, max = 1), W64=runif(nn, min = -1, max = 1), W65=runif(nn, min = -1, max = 1), 
                   W66=runif(nn, min = -1, max = 1), W67=runif(nn, min = -1, max = 1), W68=runif(nn, min = -1, max = 1), W69=runif(nn, min = -1, max = 1), W70=runif(nn, min = -1, max = 1), 
                   W71=runif(nn, min = -1, max = 1), W72=runif(nn, min = -1, max = 1), W73=runif(nn, min = -1, max = 1), W74=runif(nn, min = -1, max = 1), W75=runif(nn, min = -1, max = 1), 
                   W76=runif(nn, min = -1, max = 1), W77=runif(nn, min = -1, max = 1), W78=runif(nn, min = -1, max = 1), W79=runif(nn, min = -1, max = 1), W80=runif(nn, min = -1, max = 1), 
                   W81=runif(nn, min = -1, max = 1), W82=runif(nn, min = -1, max = 1), W83=runif(nn, min = -1, max = 1), W84=runif(nn, min = -1, max = 1), W85=runif(nn, min = -1, max = 1), 
                   W86=runif(nn, min = -1, max = 1), W87=runif(nn, min = -1, max = 1), W88=runif(nn, min = -1, max = 1), W89=runif(nn, min = -1, max = 1), W90=runif(nn, min = -1, max = 1), 
                   W91=runif(nn, min = -1, max = 1), W92=runif(nn, min = -1, max = 1), W93=runif(nn, min = -1, max = 1), W94=runif(nn, min = -1, max = 1), W95=runif(nn, min = -1, max = 1), 
                   W96=runif(nn, min = -1, max = 1), W97=runif(nn, min = -1, max = 1), W98=runif(nn, min = -1, max = 1), W99=runif(nn, min = -1, max = 1), W100=runif(nn, min = -1, max = 1))
    PS <- gbar(W)
    A = rbinom(nn, 1, PS)
    Q1 <- Qbar(rep(1, nn), W)
    Q0 <- Qbar(rep(0, nn), W)
    Y1 = rnorm(nn, Q1, 1)
    Y0 = rnorm(nn, Q0, 1)
    Y <- A * Y1 + (1-A) * Y0
    return(list(W=W,A=A,Y=Y,Y1=Y1,Y0=Y0,Q1=Q1,Q0=Q0, PS=PS))
  }
  
  sim_data_cal = sim_data  
  outcome = 'Y' 
  covars.outcome = c('A','W1','W2','W3','W4','W5','W6','W7',
                     'W8','W9','W10','W11','W12','W13',
                     'W14','W15','W16','W17','W18','W19',
                     'W20','W21','W22','W23','W24','W25',
                     'W26','W27','W28','W29','W30','W31',
                     'W32','W33','W34','W35','W36','W37',
                     'W38','W39','W40','W41','W42','W43',
                     'W44','W45','W46','W47','W48','W49',
                     'W50','W51','W52','W53','W54','W55',
                     'W56','W57','W58','W59','W60','W61',
                     'W62','W63','W64','W65','W66','W67',
                     'W68','W69','W70','W71','W72','W73',
                     'W74','W75','W76','W77','W78','W79',
                     'W80','W81','W82','W83','W84','W85',
                     'W86','W87','W88','W89','W90','W91',
                     'W92','W93','W94','W95','W96','W97',
                     'W98','W99','W100')
  
  covars.trt = c('W1','W2','W3','W4','W5','W6','W7',
                 'W8','W9','W10','W11','W12','W13',
                 'W14','W15','W16','W17','W18','W19',
                 'W20','W21','W22','W23','W24','W25',
                 'W26','W27','W28','W29','W30','W31',
                 'W32','W33','W34','W35','W36','W37',
                 'W38','W39','W40','W41','W42','W43',
                 'W44','W45','W46','W47','W48','W49',
                 'W50','W51','W52','W53','W54','W55',
                 'W56','W57','W58','W59','W60','W61',
                 'W62','W63','W64','W65','W66','W67',
                 'W68','W69','W70','W71','W72','W73',
                 'W74','W75','W76','W77','W78','W79',
                 'W80','W81','W82','W83','W84','W85',
                 'W86','W87','W88','W89','W90','W91',
                 'W92','W93','W94','W95','W96','W97',
                 'W98','W99','W100')
  outcome.type = 'continuous'   
}



