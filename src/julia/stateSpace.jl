using CSV
using DataFrames
using RData
using Turing
using Distributions
using Optim
using Memoization
using Zygote
using ReverseDiff
using FileIO

#Gompertz model
@model function gompertz(y, pr_type)
    if pr_type == 0
        #these priors work
        alpha ~ Uniform(-5,5)
        beta ~ Uniform(-1,1)
        sigma ~ Uniform(0.01,1)
        tau ~ Uniform(0.01,1)
    else
        alpha ~ Normal(0,5)
        beta ~ Normal(0,1)
        sigma ~ InverseGamma(0.01, 0.01)
        tau ~ InverseGamma(0.01, 0.01)
    end
    #Unifrom(-10,10) prior works with pr_type=0
    u_init ~ Uniform(-10,10) #not the ideal way to specify an improper prior?
    n = size(y)[1]
    umed = Vector(undef, n)
    u = Vector(undef, n)
    umed[1] = u_init
    u[1] ~ Normal(umed[1], sigma)
    for t in 2:n
        umed[t] = alpha + beta*u[t-1]
        u[t] ~ Normal(umed[t], sigma)
    end
    for t in 1:n
        y[t] ~ Normal(u[t], tau)
    end
end


#Logistic growth model
@model function logistic(y, pr_type)
    if pr_type == 0
        r ~ Uniform(0.01,5)
        K ~ Uniform(0.01,1000)
        sigma ~ Uniform(0.0001,1)
        tau ~ Uniform(0.0001,1)
    else
        r ~ LogNormal(-1,4)
        K ~ LogNormal(5,4)
        sigma ~ Exponential(0.1)
        tau ~ Exponential(0.1)
    end
    u_init ~ Uniform(0.01,100)
    n = size(y)[1]
    umed = Vector(undef, n)
    u = Vector(undef, n)
    umed[1] = u_init
    u[1] ~ LogNormal(umed[1], sigma)
    for t in 2:n
        umed[t] = u[t-1] + r*u[t-1]*(1-u[t-1]/K)
        u[t] ~ LogNormal(log(umed[t]), sigma)
    end
    for t in 1:n
        y[t] ~ LogNormal(log(u[t]), tau)
    end
end


#read in csv file: requires CSV and DataFrames
df = DataFrame(CSV.File("C:\\Users\\Andrea.Havron\\Documents\\github\\FIMSsoftware\\data\\gompertz\\gompertz_n100.csv"))
# . syntax extracts column from dataframe and casts as vector
gompertzDat = df.y
gompertzInits = load("C:\\Users\\Andrea.Havron\\Documents\\github\\FIMSsoftware\\data\\gompertz\\gompertzInits_n100.RData")
#read in csv file: requires CSV and DataFrames
df = DataFrame(CSV.File("C:\\Users\\Andrea.Havron\\Documents\\github\\FIMSsoftware\\data\\logistic\\logistic_n100.csv"))
# . syntax extracts column from dataframe and casts as vector
logisticDat = df.y
logisticInits = load("C:\\Users\\Andrea.Havron\\Documents\\github\\FIMSsoftware\\data\\logistic\\logisticInits_n100.RData")


model = gompertz(gompertzDat, 0)
mle_estimate = optimize(model,MLE())
#Default Forward mode, not timed but takes long time (~1hr)
#TrackerAD: 1791 seconds - now 754 sec using 1000 burn-in and 0.8 adaptive sampling
#Zygote: didn't work: mutating arrays not supported
#ReverseDiff: 505 sec using 1000 burn-in and 0.8 adaptive sampling
#NUTS() working but not same as default STAN/tmbstan settings (2000,0.8)
Turing.setadbackend(:reversediff)
@time sampNUTSgompertzP0 = sample(gompertz(gompertzDat,0), NUTS(2000,0.8), 4000, nchains = 4, init_theta = gompertzInits)
describe(sampNUTSgompertzP0)
save('MCMCgompertzPO.jl',sampNUTSgompertzP0)
#TrackerAD: 44.5 sec
#Zygote: didn't work: mutating arrays not supported
#ReverseDiff: 20 sec, 3140 sec after chaging settings and running first, 20 sec after running a third time: 3120 sec to compile?
@time sampNUTSgompertzP1 = sample(gompertz(gompertzDat,1), NUTS(1000,0.8), 2000, nchains = 4, init_theta = gompertzInits)
describe(sampNUTSgompertzP0)
alpha = mean(sampNUTSgompertzP0[:alpha])
beta = mean(sampNUTSgompertzP0[:beta])
sigma = mean(sampNUTSgompertzP0[:sigma])
tau = mean(sampNUTSgompertzP0[:tau])
CSV.write("MCMCgompertzP0.csv", describe(sampNUTSgompertzP0))
params = summarystats(sampNUTSgompertzP0)

ess(sampNUTSgompertzP0)
ess = ess(describe(sampNUTSgompertzP0)[1])
@time sample(gompertz(gompertzDat,0), NUTS(2000, 0.8), 2000, nchains = 4, init_theta = gompertzInits)


sampNUTSlogistic = sample(logistic(logisticDat,1), NUTS(0.65), 1000)

# MLE models

using NLsolve

@model function gompertz_mle(y, alpha, beta, u_init, u)
    n = size(y)[1]
    umed = Vector(undef, n)
    u = Vector(undef, n)
    umed[1] = u_init
    u[1] ~ Normal(umed[1], sigma)
    for t in 2:n
        umed[t] = alpha + beta*u[t-1]
        u[t] ~ Normal(umed[t], sigma)
    end
    for t in 1:n
        y[t] ~ Normal(u[t], tau)
    end
end

nlsolve(gompertz(gompertzDat, 0), gompertzInits, autodiff = :forward)
