
library(tidyr)
library(dplyr)
library(ggplot2)
library(patchwork)

plot_res <- function(methods, yvar, xvar = 'ntotal', dat, 
                     scenario, xl = '', yl = '', 
                     all.methods, idx.scen){
  # one proportion 
  # one scenario
  # subset of methods 
  # data we plot 
  dat.p <- dat %>% dplyr::filter(method %in% methods)
  dat.p <- dat.p %>% dplyr::filter(scen == scenario)

  dat.p$Y <- dat.p[[yvar]]
  dat.p$X <- dat.p[[xvar]]
  
  labs <- as.vector(all.methods[all.methods$type %in% methods, ]$label)
  
  if('Splines' %in% labs){
    labs <- labs[order(match(labs,c('Splines', 'GAM', 'Glm', 'Glmnet', 'Random Forest')))]
  }
  
  scenpaper <- idx.scen$scenpaper[idx.scen$scen == scenario]
  
  p <- ggplot(dat.p, aes(x = X, y = Y, group = method, color = method)) + 
    geom_point() + geom_line() + ggtitle(paste('Scenario', scenpaper)) + 
    xlab(xl) + ylab(yl) + theme_bw() + scale_colour_discrete("CATE estimator", labels = labs)
  
  if(yvar == 'mseratio'){
    p <- p + geom_hline(yintercept = 1, color = 'red', linetype = 'dashed') 
  }
  if(yvar == 'cal' || yvar == 'mse'){
    # New facet label names for adj variable
    adj.labs <- c("Causal Isotonic calibrator", "Uncalibrated method")
    names(adj.labs) <- c("pooled", "uncal")
    
    p <- p + geom_hline(yintercept = 0, color = 'red', linetype = 'dashed') + 
      facet_grid(~ adj, labeller = labeller(adj = adj.labs)) 
    
  }
  if(yvar == 'scaledcal'){
    dat.p <- dat.p %>% filter(adj == 'pooled')
    p <- ggplot(dat.p, aes(x = X, y = Y, group = method, color = method)) + 
      geom_point() + geom_line() + xlab(xl) + ylab(yl) + theme_bw() + 
      ggtitle(paste('Scenario', scenpaper)) + 
      geom_hline(yintercept = 0, color = 'red', linetype = 'dashed') + c
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
all.scen <- c(4, 11)
#all.scen <- c(4)

#scen.xgb <- c(4)
idx.scen <- data.frame(scen = all.scen, idxscen = 1:length(all.scen), 
                       scenpaper = 1:length(all.scen))

#idx.scen$idxxgb <- c(1, 2, -1, -1, 3, 4, -1,-1)
idx.scen$idxxgb <- c(2, -1)
#idx.scen$idxxgb <- c(2)

dir.results <- '/Users/Ernesto/Dropbox/UW/Research/Dissertation/Project 3/Alg2sims/'

## Load results
#list.files(dir.results)

res.mse <- list()
res.cal <- list()
r <- 1

for(i in all.ss){
    for(k in all.scen){
      load(paste0(dir.results, paste0(paste(i, k, sep = '-'),'-alg2.RData')))
      res.mse[[r]] <- all.res.mse
      res.cal[[r]] <- all.res.cal
      res.mse[[r]]$scen <- k
      res.cal[[r]]$scen <- k
      r <- r+1
    }
}

res.mse <- do.call(rbind, res.mse)
res.cal <- do.call(rbind, res.cal)
names(res.cal)[1] <- 'cal'
res.mse  <- res.mse %>% mutate(ntotal = nall)
res.mse  <- res.mse %>% dplyr::filter(adj != 'unpooled')

## Reconfigure MSE results
## Subset data to methods we want to compare
res.mse.cal <- res.mse %>% dplyr::filter(adj == 'pooled')
res.mse.uncal <- res.mse %>% dplyr::filter(adj == 'uncal')
temp <- merge(res.mse.cal, res.mse.uncal, by = c('method', 'ntotal', 'nk', 'scen'))

# pooled - uncal (if difference negative then our method is doing better)
temp$msedif <- temp$mse.x - temp$mse.y
temp <- temp %>% dplyr::select(method, ntotal, nk, msedif, scen)
res.mse.2 <- temp 

######################
## Calibration measure 
######################

## remove init adjustment 
res.cal <- res.cal %>% filter(adj != 'unpooled')
## compute total sample size 
res.cal <- res.cal %>% mutate(ntotal = nall)
## multiply calibration measure by n^1/3
res.cal$scaledcal <- res.cal$cal * (res.cal$nall ^ 1/3)

## separate into calibrated and uncalibrated
res.cal.cal <- res.cal %>% dplyr::filter(adj == 'pooled')
res.cal.uncal <- res.cal %>% dplyr::filter(adj == 'uncal')
# merge uncalibrated with calibrated results to compute difference between MSE 
temp <- merge(res.cal.cal, res.cal.uncal, by = c('method', 'ntotal', 'nk', 'scen'))
# Calibrated / uncalibrated
temp$calratio <- temp$cal.x/temp$cal.y
# Calibrated - uncalibrated
temp$caldif <- temp$cal.x - temp$cal.y
temp <- temp %>% dplyr::select(method, ntotal, nk, caldif, calratio, scen)

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

plots.mse <- list()
plots.cal <- list()
plots.scaled.cal <- list()

for(j in 1:length(all.scen)){
  scen <- all.scen[j]
  plots.cal[[j]] <- list()
  plots.mse[[j]] <- list()
  plots.scaled.cal[[j]] <- list()
  
  if(scen == 0){
    methods.cal <- 'glm'
    methods.mse <- 'glm'}
  if(scen == 1 || scen == 2 || scen == 3 || scen == 4 ||  scen == 9 || scen == 10){
    methods.cal <- c('gam', 'ranger', 'earth', 'glm', 'glmnet',
                     'xgboost2', 'xgboost3', 'xgboost5', 'xgboost6', 'xgboost8')
    methods.mse<- c('gam', 'ranger', 'earth', 'glm', 'glmnet',
                    'xgboost2', 'xgboost3', 'xgboost5', 'xgboost6', 'xgboost8')}
  if(scen == 5 || scen == 6 || scen == 7 || scen == 8 || scen == 11 || scen == 12){
    methods.cal <- c('glmnet', 'glmnetxgb')
    methods.mse <- c('glmnet', 'glmnetxgb')
  }
  
    plots.cal[[j]] <- plot_res(methods.cal, 'cal',
                                    dat = res.cal, scenario = scen, 
                                    xl = 'Sample Size', yl = 'Calibration measure', 
                                    all.methods = all.methods, 
                                    idx.scen = idx.scen)
    
    #plots.scaled.cal[[j]] <- plot_res(methods.cal, yvar = 'scaledcal', xvar = 'niso', 
    #                                       propt = all.proptau[i]/100, dat = res.cal, 
    #                                       scenario = scen, 
    #                                       xl = 'Sample size used to calibrate', 
    #                                       yl = 'Calibration measure x n^1/3', 
    #                                       all.methods = all.methods, 
    #                                       idx.scen = idx.scen)
    
    plots.mse[[j]] <- plot_res(methods.mse, 'mse',
                                    dat = res.mse, scenario = scen, 
                                    xl = 'Sample Size', yl = 'Mean Squared Error', 
                                    all.methods = all.methods, 
                                    idx.scen = idx.scen)
}

plots.cal[[1]]
plots.cal[[2]]
plots.mse[[1]]
plots.mse[[2]]


p.all <- plots.cal[[1]] + plots.mse[[1]] + plots.cal[[2]] + plots.mse[[2]] + plot_layout(ncol = 2,  guides = "collect")
ggsave(paste0('resultsalg2.pdf'), p.scen.1, 
       device = 'pdf', path = paste0('Figures/Manuscript Figures ICML/'), 
       width = 10, height = 10)


p.scen.1 <-  plots.cal[[1]] + plots.mse[[1]] + plot_layout(ncol = 2,  guides = "collect")

ggsave(paste0('scenario', 1,'.pdf'), p.scen.1, 
       device = 'pdf', path = paste0('Figures/Manuscript Figures ICML/'), 
       width = 10, height = 10)

p.scen.2 <- plots.cal[[2]] + plots.mse[[2]] + plot_layout(ncol = 2,  guides = "collect")
ggsave(paste0('scenario', 2,'.pdf'), p.scen.2, 
       device = 'pdf', path = paste0('Figures/Manuscript Figures ICML/'), 
       width = 10, height = 10)

