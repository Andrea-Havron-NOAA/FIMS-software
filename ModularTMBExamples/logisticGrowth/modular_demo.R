library(TMB)
source('data/simdata.R')
TMB::compile("ModularTMBExamples/logisticGrowth/logisticGrowth.cpp",
             flags="-DTMB_MODEL -w")

dyn.load(dynlib("ModularTMBExamples/logisticGrowth/logisticGrowth"))
y <- gendat(seed=123,
                  N=100,
                  theta = c(0.2,100),
                  u1 = 4,
                  var = list(proc=0.01,obs=0.001),
                  mod.name = 'logistic')
dat <- list(y=y)
par <- list(theta = c(log(0.5), log(80)), ln_sig=-1,ln_tau=-1,
            u = rep(1,length(y)))
a <- Sys.time()
obj <- MakeADFun(dat, par, random="u", DLL="logisticGrowth")
newtonOption(obj, smartsearch=FALSE)
obj$fn()
obj$gr()
system.time(opt <- nlminb(obj$par, obj$fn, obj$gr))
sdr <- sdreport(obj)
b <- Sys.time()
difftime(b,a)
summary(sdr, 'fixed')
summary(sdr, 'random')[,1]
summary(sdr, 'report')

report <- obj$report()
report$u
report$K
report$r

#simulate data
set.seed(123)

dyn.unload(dynlib("ModularTMBExamples/logisticGrowth/logisticGrowth"))

## not working yet
# detach("package:TMB", unload=TRUE)
# library(R2admb)
# setwd('ModularTMBExamples/logisticGrowth/admb/')
# write_dat('logisticGrowth', dat)
# par <- list(ln_r = log(0.5), ln_K = log(80), ln_sig=-1,ln_tau=-1,
#             u = rep(1,length(y)))
# write_pin('logisticGrowth', par)
# compile_admb('logisticGrowth')#, re = TRUE, verbose = TRUE)
# a <- Sys.time()
# admb.mod <- run_admb('logisticGrowth', verbose = TRUE, extra.args = '-noinit')
# b <- Sys.time()
