

library(tidyr)
library(dplyr)
library(ggplot2)
library(patchwork)

## Update Accordingly
## Create list of results 
all.ss <- c(1000, 2000, 5000, 10000)

all.scen <- c(4, 11)
scen.xgb <- c(4)
idx.scen <- data.frame(scen = all.scen, idxscen = 1:length(all.scen), 
                       scenpaper = 1:length(all.scen))
idx.scen$idxxgb <- c(2,-1)

dir.results <- '/Users/Ernesto/Dropbox/UW/Research/Dissertation/Project 3/Alg2sims/'

## Load results
#list.files(dir.results)

res.quant.upper <- list()      
res.quant.lower <- list()     

res.quant.upper.tau0 <- list()
res.quant.lower.tau0 <- list()

all.quant.upper.cal <- list()
all.quant.lower.cal <- list()

all.quant.upper.uncal <- list()
all.quant.lower.uncal <- list()

all.quant.lower.tau0 <- list()
all.quant.upper.tau0 <- list()

r <- 1

for(i in all.ss){
    for(k in all.scen){
      load(paste0(dir.results, paste0(paste(i, k, sep = '-'),'-alg2.RData')))
      
      res.quant.upper[[r]] <- all.res.quant.upper         
      res.quant.lower[[r]] <-  all.res.quant.lower        
      res.quant.upper.tau0[[r]] <- all.res.quant.upper.tau0
      res.quant.lower.tau0[[r]] <- all.res.quant.lower.tau0
      
      res.quant.upper[[r]]$scen <- k
      res.quant.lower[[r]]$scen <- k
      res.quant.upper.tau0[[r]]$scen <- k
      res.quant.lower.tau0[[r]]$scen <- k 
      
      # for plots
      all.quant.upper.cal[[r]] <- data.frame(all_res$all.quantile.upper.cal)     
      all.quant.lower.cal[[r]] <- data.frame(all_res$all.quantile.lower.cal)      
      all.quant.upper.uncal[[r]] <- data.frame(all_res$all.quantile.upper.uncal)    
      all.quant.lower.uncal[[r]] <- data.frame(all_res$all.quantile.lower.uncal)   
      all.quant.lower.tau0[[r]] <- data.frame(all_res$all.quantile.lower.tau0.cal)
      all.quant.upper.tau0[[r]] <- data.frame(all_res$all.quantile.upper.tau0.cal)
      
      all.quant.upper.cal[[r]]$scen <- k 
      all.quant.lower.cal[[r]]$scen <- k  
      all.quant.upper.uncal[[r]]$scen <- k 
      all.quant.lower.uncal[[r]]$scen <- k 
      all.quant.lower.tau0[[r]]$scen <- k 
      all.quant.upper.tau0[[r]]$scen <- k 
      
      all.quant.upper.cal[[r]]$all.ss <- i 
      all.quant.lower.cal[[r]]$all.ss <- i  
      all.quant.upper.uncal[[r]]$all.ss <- i 
      all.quant.lower.uncal[[r]]$all.ss <- i 
      all.quant.lower.tau0[[r]]$all.ss <- i
      all.quant.upper.tau0[[r]]$all.ss <- i 
      
      r <- r+1
    }
  }

rm(all_res)


res.quant.upper <- do.call(rbind, res.quant.upper)     
res.quant.lower <- do.call(rbind, res.quant.lower)       
res.quant.upper.tau0 <- do.call(rbind, res.quant.upper.tau0)   
res.quant.lower.tau0 <- do.call(rbind, res.quant.lower.tau0)   
names(res.quant.lower.tau0)[1] <- 'lowerqtau0'

# merge res.quant.upper and res.quant.upper.tau0 by proptau, scen, niso, adj, method
all.upper <- merge(res.quant.upper, res.quant.upper.tau0, 
                   by = c('scen', 'nall', 'adj', 'method', 'nk'))

# merge res.quant.lower and res.quant.lower.tau0 by proptau, scen, niso, adj, method 
all.lower <- merge(res.quant.lower, res.quant.lower.tau0, 
                   by = c('scen', 'nall', 'adj', 'method', 'nk'))

##### TO DO: Remove filter by proptau

# compute difference 
all.upper <- all.upper %>% mutate(difupperq = upperq - upperqtau0)
all.lower <- all.lower %>% mutate(diflowerq = lowerq - lowerqtau0)
all.upper <- all.upper %>% mutate(reldifupperq = (upperq - upperqtau0)/upperqtau0)
all.lower <- all.lower %>% mutate(reldiflowerq = (lowerq - lowerqtau0)/lowerqtau0)

# round to two decimals
all.upper <- all.upper %>% mutate(difupperq = round(difupperq, 2))
all.lower <- all.lower %>% mutate(diflowerq = round(diflowerq, 2 ))
all.upper <- all.upper %>% mutate(reldifupperq = round(reldifupperq, 2))
all.lower <- all.lower %>% mutate(reldiflowerq = round(reldiflowerq, 2 ))

# merge tables 
all.quant <- merge(all.upper, all.lower, by = c('scen', 'nall', 'adj', 'method', 'nk'))
# remove 'init'
all.quant <- all.quant %>% filter(adj != 'unpooled')

# add total sample size 

# remove unecessary columns 
# subset to keep only difference 
# otherwise table is too big [(upper, upper, dif) + (lowe, lower, dif)] x 4 = 24
all.quant <- all.quant %>% dplyr::select(-c(nk, upperq, upperqtau0, lowerq, lowerqtau0))

# make table 
small.dim.bin <- all.quant %>% filter(scen %in% c(4)) %>% dplyr::select(-c(reldiflowerq, reldifupperq))
large.dim.bin <- all.quant %>% filter(scen %in% c(11)) %>% dplyr::select(-c(diflowerq, difupperq))

small.dim.bin <- small.dim.bin %>%  pivot_wider(names_from = nall, 
                                                values_from = c(difupperq, diflowerq))

large.dim.bin <- large.dim.bin %>%  pivot_wider(names_from = nall, 
                                                values_from = c(reldifupperq, reldiflowerq))

# select methods we are interested in showing [glmnet, splines, ranger]

# small.dim.bin <- small.dim.bin %>% filter(method %in% c('ranger', 'earth', 'glmnet')) 
# large.dim.bin <- large.dim.bin %>% filter(method %in% c('glmnet', 'glmnetxgb')) 

# rearrange columns and rows 
out.small.bin <- data.frame(small.dim.bin)
out.small.bin <- out.small.bin[order(out.small.bin$scen, out.small.bin$method),]
out.small.bin <- out.small.bin[,c('scen', 'adj', 'method', 
                                  'diflowerq_1000', 'difupperq_1000', 
                                  'diflowerq_2000', 'difupperq_2000',
                                  'diflowerq_5000','difupperq_5000',
                                  'diflowerq_10000',  'difupperq_10000')]


out.large.bin <- data.frame(large.dim.bin)
out.large.bin <- out.large.bin[order(out.large.bin$scen, out.large.bin$method),]
out.large.bin <- out.large.bin[, c('scen', 'adj', 'method', 
                                  'reldiflowerq_1000', 'reldifupperq_1000',
                                  'reldiflowerq_2000', 'reldifupperq_2000',
                                  'reldiflowerq_5000', 'reldifupperq_5000',
                                  'reldiflowerq_10000', 'reldifupperq_10000')]

# make table as csv

names(out.small.bin)[4:11] <- c('lower10k', 'upper10k', 'lower5k', 'upper5k',
                                'lower1k', 'upper1k', 'lower2k', 'upper2k')
names(out.large.bin)[4:11] <- c('lower10k', 'upper10k', 'lower5k', 'upper5k',
                                'lower1k', 'upper1k', 'lower2k', 'upper2k')

all.bins <- rbind(out.small.bin, out.large.bin)
all.bins <- all.bins[,c(1,2,3,8,9,10,11,6,7,4,5)]
all.bins

write.csv(all.bins, 'Tables/Tables ICML/allbins.csv', row.names = FALSE)



