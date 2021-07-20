tuna.dat
15.9 61.89
25.7 78.98
28.5 55.59
23.7 44.61
25.0 56.89
33.3 38.27
28.2 33.84
19.7 36.13
17.5 41.95
19.3 36.63
21.6 36.33
23.1 38.82
22.5 34.32
22.5 37.64
23.6 34.01
29.1 32.16
14.4 26.88
13.2 36.61
28.4 30.07
34.6 30.75
37.5 23.36
25.9 22.36
25.3 21.91
tuna_S.dat
list(C=c(15.9,25,7,...,25.3), I=c(61.89,78.98,...,21.91))
surplus.in
list(P=c(0.99,0.98,0.96,0.94,0.92,0.90,0.88,0.86,0.84,0.82,
         0.80,0.78,0.76,0.74,0.72,0.70,0.68,0.66,0.64,0.62,0.60,0.5
         8,0.56),
     r=0.8, K=200, iq=5, isigma2=100, itau2=100)
surplus.bug
model surplusproduction;
const N=23;
var C[N], I[N], Imed[N], P[N], Pmed[N],
r, K, q, iq, sigma2, isigma2, tau2, itau2, MSP, EMSP,
P1990, B1990;
data C, I in “tuna.dat”;
inits in “surplus.in”;
{
  # prior distribution of K: lognormal with 10% and 90%
  quantile at 80 and 300
  K ~ dlnorm(5.042905,3.7603664)I(10,1000);
  # prior distribution of r: lognormal with 10% and 90%
  quantile at 0.13 and 0.48
  r ~ dlnorm(–1.38,3.845)I(0.01,1.2);
  # prior distribution of q: instead of improper (prop. to 1/q)
  use just proper IG
  iq ~ dgamma(0.001,0.001)I(0.5,100);
  q <- 1/iq;
  # prior distribution of sigma2: inv. gamma with 10% and
  90% qu. at 0.04 and 0.08
  isigma2 ~ dgamma(3.785518,0.010223);
  sigma2 <- 1/isigma2;
  # prior distribution of tau2: inv. gamma with 10% and
  90% qu. at 0.05 and 0.15
  itau2 ~ dgamma(1.708603,0.008613854);
  tau2 <- 1/itau2;
  # (conditional) prior distribution of Ps (from state equations):
  Pmed[1] <- 0;
  P[1] ~ dlnorm(Pmed[1],isigma2) I(0.001,2.0)
  for (t in 2:N) { Pmed[t] <- log(P[t–1] + r*P[t–1]*(1-P[t–1])
                                  - C[t–1]/K);
  P[t] ~ dlnorm(Pmed[t],isigma2)I(0.001,2.0) }
  # sampling distribution:
  for (t in 1:N) { Imed[t] <- log(q*K*P[t]);
  I[t] ~ dlnorm(Imed[t],itau2);
  }
  # further management parameters and predictions:
  MSP <- r*K/4;
  EMSP <- r/(2*q);
  P1990 <- P[N]+r*P[N]*(1-P[N]) -C[N]/K;
  B1990 <- P1990*K;
}
surplus.cmd
compile(“surplus.bug”)
update(25000)
monitor(K,25)
monitor(r,25)
monitor(q,25)
monitor(sigma2,25)
monitor(tau2,25)
monitor(P[],25)
monitor(P1990,25)
monitor(B1990,25)
monitor(MSP,25)
monitor(EMSP,25)
update(225000)
stats(K)
stats(r)
stats(q)
stats(sigma2)
stats(tau2)
stats(P[])
stats(P1990)
stats(B1990)
stats(MSP)
stats(EMSP)
q()