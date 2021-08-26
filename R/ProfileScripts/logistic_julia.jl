using CSV
using DataFrames
using RData
using Turing
using Distributions
using ReverseDiff

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
df = DataFrame(CSV.File("C:\\Users\\Andrea.Havron\\Documents\\github\\FIMSsoftware\\data\\logistic\\logistic_n100.csv"))
# . syntax extracts column from dataframe and casts as vector
logisticDat = df.y
logisticInits = load("C:\\Users\\Andrea.Havron\\Documents\\github\\FIMSsoftware\\data\\logistic\\logisticInits_n100.RData")

#Set AD to ReverseDiff
Turing.setadbackend(:reversediff)
#run model using NUTS and ReverseDiff
@time sampNUTSlogistic = sample(logistic(logisticDat,1), NUTS(2000,0.8), 4000, nchains = 4, init_theta = logisticInits)
#results
r = mean(sampNUTSlogistic[:r])
K = mean(sampNUTSlogistic[:K])
sigma = mean(sampNUTSlogistic[:sigma])
tau = mean(sampNUTSlogistic[:tau])
#ess = ess(describe(sampNUTSlogistic)[1])
