n.seq <- seq(5,11,1)
library(tidyr)
library(magrittr)
library(ggplot2)

mcmc.res <- c()
tmb.est <- tmbstan.est <- stan.est <- admb.est <- c()
tmb.se <- tmbstan.se <- stan.se <- admb.se <- c()
time.est <- c()
for(i in 1:length(n.seq)){
  n <- 2^n.seq[i]
  load( paste0('results/logistic/logistic', '_n', n, '.RData'))
  
  mcmc.res <- rbind(mcmc.res, sapply(logistic.results[2:3], function(x) x$meanESS)/
                      sapply(logistic.results[2:3], function(x) x$time))
  mcmc.res <- rbind(mcmc.res, sapply(logistic.results[2:3], function(x) x$minESS)/
                      sapply(logistic.results[2:3], function(x) x$time))
  mcmc.res <- rbind(mcmc.res, sapply(logistic.results[2:3], function(x) x$time))
  tmb.est <- rbind(tmb.est, logistic.results$tmb$par.est)
  tmbstan.est <- rbind(tmbstan.est, logistic.results$tmbstan$par.est)
  stan.est <- rbind(stan.est, logistic.results$stan$par.est)
  tmb.se <- rbind(tmb.se, logistic.results$tmb$se.est)
  tmbstan.se <- rbind(tmbstan.se, logistic.results$tmbstan$se.est)
  stan.se <- rbind(stan.se, logistic.results$stan$se.est)
  
  admb.est <- rbind(admb.est, c(logistic.results$admb$logisticGrowth$par.est, 'logisticGrowth'))
  admb.est <- rbind(admb.est, c(logistic.results$admb$logisticGrowth_nonsep$par.est[1:4], 'logisticGrowth_nonsep'))
  admb.est <- rbind(admb.est, c(logistic.results$admb$logisticGrowth_PQL$par.est, 'logisticGrowth_PQL'))
  time.vec <- c(sapply(logistic.results[1:3], function(x) x$time), 
                logistic.results$admb$logisticGrowth$time, 
                logistic.results$admb$logisticGrowth_nonsep$time, 
                logistic.results$admb$logisticGrowth_PQL$time)
  
  #admb.se <- rbind(admb.est, logistic.results$admb$se.est)
  time.est <- rbind(time.est, time.vec )
}
colnames(mcmc.res) <- names(logistic.results)[2:3]
mcmc.res %<>% as.data.frame()
mcmc.res$metric <- rep(c('mean MCMC eff', 'min MCMC eff', 'minutes'), length(n.seq))
mcmc.res$nsamp <- rep(2^n.seq, each = 3)
png(filename = 'results/plots/logistic_mcmc.png', width = 1200, height = 500, res=200)
mcmc.res %>% 
  pivot_longer(., 1:2, names_to = 'model', values_to = 'value') %>%
  ggplot(., aes(x=nsamp, y=value,col=model)) + geom_line() + 
  theme_classic() + facet_wrap(~metric, scales = 'free', nrow=1) 
dev.off()

tmb.est <- as.data.frame(tmb.est)
tmbstan.est <- as.data.frame(tmbstan.est)
stan.est <- as.data.frame(stan.est)
admb.est <- as.data.frame(admb.est)
colnames(admb.est) <-  c(colnames(tmb.est),'model')
tmb.est$model <- rep('tmb-re',8)
tmb.est$nsamp <- 2^n.seq
tmbstan.est$model <- rep('tmbstan',8)
tmbstan.est$nsamp <- 2^n.seq
stan.est$model <- rep('stan',8)
stan.est$nsamp <- 2^n.seq
admb.est$nsamp <- rep(2^n.seq, each = 3)
colnames(stan.est) <- colnames(tmb.est)
colnames(tmbstan.est) <- colnames(tmb.est)

est.res <- rbind(tmb.est, tmbstan.est, stan.est, admb.est[admb.est$model == 'logisticGrowth',])
est.res %<>% pivot_longer(.,1:4,names_to = 'parm', values_to = 'est') %>% as.data.frame()
levels(est.res$model)
png(filename = 'results/plots/logistic_est.png', width = 1500, height = 1000, res = 200)
ggplot(est.res, aes(x=nsamp, y = as.numeric(est))) + 
  facet_grid(parm~model, scales = "free_y") + 
  geom_point() + theme_bw() +
  geom_hline(data=dplyr::filter(est.res, parm=="K"), aes(yintercept=100), color = 'blue')+
  geom_hline(data=dplyr::filter(est.res, parm=="r"), aes(yintercept=0.2), color = 'blue')+
  geom_hline(data=dplyr::filter(est.res, parm=="sigma"), aes(yintercept=sqrt(0.01)), color = 'blue')+
  geom_hline(data=dplyr::filter(est.res, parm=="tau"), aes(yintercept=sqrt(0.001)), color = 'blue')
dev.off()

est.res <- admb.est[admb.est$model != 'logisticGrowth_PQL',]#rbind(tmb.est, tmbstan.est, stan.est, admb.est)
est.res %<>% pivot_longer(.,1:4,names_to = 'parm', values_to = 'est') %>% as.data.frame()
est.res$model <- factor(est.res$model, levels = c("logisticGrowth", "logisticGrowth_nonsep"), 
                        labels = c("Separable", "Non-Separable"))

png(filename = 'results/plots/logistic_admb_sep_est.png', width = 800, height = 1000, res = 200)
ggplot(est.res, aes(x=nsamp, y = as.numeric(est))) + 
  facet_grid(parm~model, scales = "free_y") + 
  geom_point() + theme_bw() + ylab("")+
  geom_hline(data=dplyr::filter(est.res, parm=="K"), aes(yintercept=100), color = 'blue')+
  geom_hline(data=dplyr::filter(est.res, parm=="r"), aes(yintercept=0.2), color = 'blue')+
  geom_hline(data=dplyr::filter(est.res, parm=="sigma"), aes(yintercept=sqrt(0.01)), color = 'blue')+
  geom_hline(data=dplyr::filter(est.res, parm=="tau"), aes(yintercept=sqrt(0.001)), color = 'blue')
dev.off()
time.df <- as.data.frame(time.est)
colnames(time.df) <- c(colnames(time.est)[1:3], 'Separable', 'Non-separable', 'PQL')
admb_sep_time <- time.df[,4:5]
admb_sep_time$nsamp <- 2^(5:12)

png(filename = 'results/plots/logistic_admb_sep_time.png', width = 1000, height = 700, res = 200)
admb_sep_time %>% pivot_longer(.,1:2, names_to='model', values_to='minutes') %>%
  ggplot(., aes(x=nsamp, y=minutes, color = model)) + geom_line() + theme_classic()
dev.off()

est.res <- rbind(tmb.est,admb.est[admb.est$model == 'logisticGrowth',])
est.res %<>% pivot_longer(.,1:4,names_to = 'parm', values_to = 'est') %>% as.data.frame()
est.res$model <- factor(est.res$model, levels = c("logisticGrowth", "tmb-re"), 
                        labels = c("ADMB", "TMB"))
png(filename = 'results/plots/logistic_admb_tmb_est.png', width = 800, height = 1000, res = 200)
ggplot(est.res, aes(x=nsamp, y = as.numeric(est))) + 
  facet_grid(parm~model, scales = "free_y") + 
  geom_point() + theme_bw() + ylab("")+
  geom_hline(data=dplyr::filter(est.res, parm=="K"), aes(yintercept=100), color = 'blue')+
  geom_hline(data=dplyr::filter(est.res, parm=="r"), aes(yintercept=0.2), color = 'blue')+
  geom_hline(data=dplyr::filter(est.res, parm=="sigma"), aes(yintercept=sqrt(0.01)), color = 'blue')+
  geom_hline(data=dplyr::filter(est.res, parm=="tau"), aes(yintercept=sqrt(0.001)), color = 'blue')
dev.off()
admb_sep_time
time.est
time.est <- data.frame(TMB = time.est$TMB, ADMB = admb_sep_time$Separable, nsamp = admb_sep_time$nsamp)
png(filename = 'results/plots/logistic_admb_tmb_time.png', width = 1000, height = 700, res = 200)
time.est %>% pivot_longer(.,1:2, names_to='model', values_to='minutes') %>%
  ggplot(., aes(x=nsamp, y=minutes, color = model)) + geom_line() + theme_classic()
dev.off()


est.res <- rbind(tmbstan.est, stan.est)
est.res %<>% pivot_longer(.,1:4,names_to = 'parm', values_to = 'est') %>% as.data.frame()
est.res$model <- factor(est.res$model, levels = c("tmbstan", "stan"), 
                        labels = c("TMB", "Stan"))
png(filename = 'results/plots/logistic_stan_tmb_est.png', width = 800, height = 1000, res = 200)
ggplot(est.res, aes(x=nsamp, y = as.numeric(est))) + 
  facet_grid(parm~model, scales = "free_y") + 
  geom_point() + theme_bw() + ylab("")+
  geom_hline(data=dplyr::filter(est.res, parm=="K"), aes(yintercept=100), color = 'blue')+
  geom_hline(data=dplyr::filter(est.res, parm=="r"), aes(yintercept=0.2), color = 'blue')+
  geom_hline(data=dplyr::filter(est.res, parm=="sigma"), aes(yintercept=sqrt(0.01)), color = 'blue')+
  geom_hline(data=dplyr::filter(est.res, parm=="tau"), aes(yintercept=sqrt(0.001)), color = 'blue')
dev.off()


