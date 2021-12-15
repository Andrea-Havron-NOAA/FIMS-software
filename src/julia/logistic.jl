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
    u_init ~ Uniform(2,100)
    n = size(y)[1]
    umed = Vector(undef, n)
    ln_u = Vector(undef, n)
    umed[1] = u_init
    ln_u[1] ~ Normal(umed[1], sigma)
    for t in 2:n
        umed[t] = exp(ln_u[t-1]) + r*exp(ln_u[t-1])*(1-exp(ln_u[t-1])/K)
        ln_u[t] ~ Normal(umed[t], sigma)
    end
    for t in 1:n
        y[t] ~ LogNormal(ln_u[t], tau)
    end
end

#read in csv file: requires CSV and DataFrames
df = DataFrame(CSV.File("C:\\Users\\Andrea.Havron\\Documents\\github\\FIMSsoftware\\data\\logistic\\logistic_n32.csv"))
# . syntax extracts column from dataframe and casts as vector
logisticDat = df.y
logisticInits = load("C:\\Users\\Andrea.Havron\\Documents\\github\\FIMSsoftware\\data\\logistic\\logisticInits_n32.RData")

model = logistic(logisticDat, 1)
mle_estimate = optimize(model,MLE())

#Default Forward mode, not timed but takes long time (~1hr)
#TrackerAD: 1791 seconds - now 754 sec using 1000 burn-in and 0.8 adaptive sampling
#Zygote: didn't work: mutating arrays not supported
#ReverseDiff: 505 sec using 1000 burn-in and 0.8 adaptive sampling
#NUTS() working but not same as default STAN/tmbstan settings (2000,0.8)
Turing.setadbackend(:reversediff)

@time sampNUTSlogistic = sample(logistic(logisticDat,1), NUTS(2000,0.8), 4000, nchains = 4, init_theta = logisticInits)
describe(sampNUTSlogistic)
CSV.write("MCMClogistic_n256.csv", describe(sampNUTSlogistic)[1])
