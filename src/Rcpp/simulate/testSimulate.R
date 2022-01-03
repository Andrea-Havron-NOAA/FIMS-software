library(TMB)
TMB::compile('src/Rcpp/simulate/testSimulate.cpp')

x <- 1:100
a <- 2
b <- -2
sig <- 1

mu <- a + x*b
ysim <- rnorm(100,mu,sig)

dyn.load(dynlib('src/Rcpp/simulate/testSimulate'))
obj <- MakeADFun(
  data = list(
    y=ysim, x=x
  ),
  parameters = list(
    a=0, b=0, sd = 1
  ) )
opt <- nlminb(obj$par, obj$fn, obj$gr)
opt$par
plot(ysim, obj$simulate()$y)

dyn.unload(dynlib('src/Rcpp/simulate/testSimulate'))

