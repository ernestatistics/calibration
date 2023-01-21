
library(tidyr)
library(dplyr)
library(ggplot2)
library(patchwork)

plot_res <- function(methods, yvar, xvar = 'ntotal', 
                     propt, dat, scenario, xl = '', yl = '', 
                     all.methods, idx.scen){
  # one proportion 
  # one scenario
  # subset of methods 
  # data we plot 
  dat.p <- dat %>% dplyr::filter(method %in% methods)
  dat.p <- dat.p %>% dplyr::filter(scen == scenario)
  dat.p <- dat.p %>% dplyr::filter(proptau == propt)
  
  dat.p$Y <- dat.p[[yvar]]
  dat.p$X <- dat.p[[xvar]]
  
  labs <- as.vector(all.methods[all.methods$type %in% methods, ]$label)
  orig_labs <- all.methods[all.methods$type %in% methods, ]$type
  
  if('Splines' %in% labs){
    labs <- labs[order(match(labs,c('Splines', 'GAM', 'Glm', 'Glmnet', 'Random Forest')))]
    orig_labs <- all.methods[all.methods$type %in% methods, ]$type
  }
  
  scenpaper <- idx.scen$scenpaper[idx.scen$scen == scenario]
  

  
  
  #p <- ggplot(dat.p, aes(x = X, y = Y, group = method, color = method)) + 
    #geom_point() + geom_line() + 
    #ggtitle(paste0('Scenario', scenpaper, ':','Proportion of sample to fit CATE estimator:', propt)) + xlab(xl) + ylab(yl) + theme_bw() + 
    #theme(axis.text.x=element_text(angle=30,hjust=1)) + scale_colour_discrete("CATE estimator", 
    #                                                                          breaks = orig_labs,
    #                                                                          labels = labs) 
  
  
  p <- ggplot(dat.p, aes(x = X, y = Y, group = method, color = method)) + 
  geom_point() + geom_line() + 
  ggtitle(expression(paste('Proportion of sample to fit ', hat(tau), ': 0.8'))) + xlab(xl) + ylab(yl) + theme_bw() + 
  theme(axis.text.x=element_text(angle=30,hjust=1)) + 
    scale_colour_discrete("CATE estimator", breaks = orig_labs,labels = labs) 
  
  
  if(yvar == 'mseratio'){
    p <- p + geom_hline(yintercept = 1, color = 'red', linetype = 'dashed') 
  }
  if(yvar == 'cal'|| yvar == 'mse'){
    # New facet label names for adj variable
    adj.labs <- c("Causal Isotonic calibrator", "Uncalibrated method")
    names(adj.labs) <- c("cal", "uncal")
    
    p <- p + geom_hline(yintercept = 0, color = 'red', linetype = 'dashed') + 
      facet_grid(~ adj, labeller = labeller(adj = adj.labs)) 
    
  }
  if(yvar == 'scaledcal'){
    dat.p <- dat.p %>% filter(adj == 'cal')
    p <- ggplot(dat.p, aes(x = X, y = Y, group = method, color = method)) + 
      geom_point() + geom_line() + 
      ggtitle(ggtitle(expression(paste('Proportion of sample to fit ', hat(tau), ': 0.8')))) + 
      xlab(xl) + ylab(yl) + theme_bw() + 
      geom_hline(yintercept = 0, color = 'red', linetype = 'dashed')  + c
  }
  if(yvar == 'msedif'){
    p <- p + geom_hline(yintercept = 0, color = 'red', linetype = 'dashed')
  }
  if(xvar == 'ntotal'){
    p <- p + scale_x_continuous(name = xl, breaks = c(1000, 2000, 5000, 10000), 
                                labels = c('1000', '2000', '5000', '10000'))
  }
  return(p)
}


## Update Accordingly
## Create list of results 
all.ss <- c(1000, 2000, 5000, 10000)
#all.proptau <- c(20, 50, 80)
all.proptau <- c(80)
#all.scen <- c(3, 4, 7, 8, 9, 10, 11, 12)
all.scen <- c(4,11)
idx.scen <- data.frame(scen = all.scen, idxscen = 1:length(all.scen), 
                       scenpaper = 1:length(all.scen))
#idx.scen$idxxgb <- c(1, 2, -1, -1, 3, 4, -1,-1)
idx.scen$idxxgb <- c(2, -1)
dir.results <- '/Users/Ernesto/Dropbox/UW/Research/Dissertation/Project 3/Simulation Results/Results 05 09/'

## Load results
#list.files(dir.results)

res.mse <- list()
res.cal <- list()
r <- 1

for(i in all.ss){
  for(j in all.proptau){
    for(k in all.scen){
      load(paste0(dir.results, paste0(paste(i, j, k, sep = '-'),'.RData')))
      res.mse[[r]] <- res_ss$all.res.mse
      res_ss$all.res.cal <- res_ss$all.res.cal.2 
      res.cal[[r]] <- res_ss$all.res.cal
      res.mse[[r]]$scen <- k
      res.cal[[r]]$scen <- k
      r <- r+1
    }
  }
}
rm(res_ss)

res.mse <- do.call(rbind, res.mse)
res.cal <- do.call(rbind, res.cal)
names(res.cal)[1] <- 'cal'

## Reconfigure MSE results
res.mse  <- res.mse %>% mutate(ntotal = ntau + niso)
res.mse  <- res.mse %>% dplyr::filter(adj != 'init')

## Subset data to methods we want to compare
res.mse.cal <- res.mse %>% dplyr::filter(adj == 'cal')
res.mse.uncal <- res.mse %>% dplyr::filter(adj == 'uncal')
temp <- merge(res.mse.cal, res.mse.uncal, by = c('method', 'ntau', 
                                                 'niso', 'scrossf', 
                                                 'proptau', 'scen'))
temp$msedif <- temp$mse.x - temp$mse.y
temp <- temp %>% dplyr::select(method, ntau, niso, scrossf, proptau, msedif, scen)
res.mse.2 <- temp 
res.mse.2  <- res.mse.2 %>% mutate(ntotal = ntau + niso)

######################
## Calibration measure 
######################

## remove init adjustment 
res.cal <- res.cal %>% filter(adj != 'init')
## compute total sample size 
res.cal <- res.cal %>% mutate(ntotal = ntau + niso)
## multiply calibration measure by n^1/3
res.cal$scaledcal <- res.cal$cal * (res.cal$niso ^ 1/3)

## separate into calibrated and uncalibrated
res.cal.cal <- res.cal %>% dplyr::filter(adj == 'cal')
res.cal.uncal <- res.cal %>% dplyr::filter(adj == 'uncal')
# merge uncalibrated with calibrated results to compute difference between MSE 
temp <- merge(res.cal.cal, res.cal.uncal, by = c('method', 'ntau', 'niso', 'scrossf', 'proptau', 'scen', 'ntotal'))
# Calibrated / uncalibrated
temp$calratio <- temp$cal.x/temp$cal.y
# Calibrated - uncalibrated
temp$caldif <- temp$cal.x - temp$cal.y
temp <- temp %>% dplyr::select(method, ntau, niso, scrossf, proptau, caldif, calratio, scen)
res.cal.2 <- temp
res.cal.2 <- gather(res.cal.2, adj, cal, caldif:calratio, factor_key=TRUE)
res.cal.plot <- bind_rows(res.cal,res.cal.2)

all.methods <- data.frame(type = c('gam', 'ranger', 'earth', 'glm', 'glmnet',
                                   'glmnetxgb','xgboost2', 
                                   'xgboost3', 'xgboost5', 
                                   'xgboost6', 'xgboost8'), 
                          label = c('GAM', 'Random Forest', 'Splines', 'Glm', 'Glmnet', 
                                    'Glmnet + \n Boosted Trees',
                                    'Boosted Trees 2', 
                                    'Boosted Trees 3', 
                                    'Boosted Trees 5', 
                                    'Boosted Trees 6', 
                                    'Boosted Trees 8'))

# New facet label names for scenario variable
scenario.labs <- c("Scenario 1", "Scenario 2")
names(scenario.labs) <- c("4", "11")

# New facet label names for supp variable
prop.labs <- c('0.8')
names(prop.labs) <- c('0.8')

# New facet label names for method variable
adj.labs <- c("Causal Isotonic calibrator Sample Splitting", "Uncalibrated method")
names(adj.labs) <- c("cal", "uncal")
mthds <- all.methods$label

plots.mse <- list()
plots.cal <- list()
for(j in 1:length(all.scen)){
  
  scen <- all.scen[j]
  plots.cal[[j]] <- list()
  plots.mse[[j]] <- list()
  
  if(scen == 0){
    methods.cal <- 'glm'
    methods.mse <- 'glm'}
  if(scen == 1 || scen == 2 || scen == 3 || scen == 4 ||  scen == 9 || scen == 10){
    methods.cal <- c('gam', 'ranger', 'earth', 'glm', 'glmnet',
                     'xgboost2', 'xgboost3', 'xgboost5', 'xgboost6', 'xgboost8')
    methods.mse<- c('gam', 'ranger', 'earth', 'glm', 'glmnet',
                    'xgboost2', 'xgboost3', 'xgboost5', 'xgboost6', 'xgboost8')}
  
  if(scen == 5 || scen == 6 || scen == 7 || scen == 8 || scen == 11 || scen == 12){
    methods.cal <- c('glmnet', 'glmnetxgb', 'ranger')
    methods.mse <- c('glmnet', 'glmnetxgb', 'ranger')
  }
  
  for(i in 1:length(all.proptau)){
    plots.cal[[j]][[i]] <- plot_res(methods.cal, 'cal', propt = all.proptau[i]/100, 
                                    dat = res.cal, scenario = scen, 
                                    xl = 'Sample Size', yl = 'Calibration measure', 
                                    all.methods = all.methods, 
                                    idx.scen = idx.scen)
    
    plots.mse[[j]][[i]] <- plot_res(methods.mse, 'mse', propt = all.proptau[i]/100,
                                    dat = res.mse, scenario = scen, 
                                    xl = 'Sample Size', yl = 'Mean Squared Error', 
                                    all.methods = all.methods, 
                                    idx.scen = idx.scen)
    #plots.ate[[i]] <- plot_by_ss(methods, all.proptau[i], res.ate)
  }
}


p.all <- plots.cal[[1]][[1]] + plots.mse[[1]][[1]] + plots.cal[[2]][[1]] + plots.mse[[2]][[1]] + 
  plot_layout(ncol = 2,  guides = "collect") 
     
ggsave(paste0('results-ss-80.pdf'), p.all, 
       device = 'pdf', path = paste0('Figures/Manuscript Figures ICML/'), 
       width = 10, height = 7)
  
  


