library(TMB)
compile("src/Rcpp/logisticGrowth.cpp", flags= "-w")
dyn.load(dynlib("logisticGrowth"))