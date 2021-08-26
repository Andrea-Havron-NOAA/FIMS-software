using CSV
using DataFrames
using RData
using Turing
using Distributions
using ReverseDiff

#Logistic growth model
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

#read in csv file: requires CSV and DataFrames
df = DataFrame(CSV.File("C:\\Users\\Andrea.Havron\\Documents\\github\\FIMSsoftware\\data\\gompertz\\gompertz_n100.csv"))
# . syntax extracts column from dataframe and casts as vector
logisticDat = df.y
logisticInits = load("C:\\Users\\Andrea.Havron\\Documents\\github\\FIMSsoftware\\data\\gompertz\\gompertzInits_n100.RData")

#Set AD to ReverseDiff
Turing.setadbackend(:reversediff)
#run model using NUTS and ReverseDiff
@time sampNUTSgompertz = sample(logistic(gompertzDat,1), NUTS(2000,0.8), 4000, nchains = 4, init_theta = gompertzInits)
#results
r = mean(sampNUTSgompertz[:r])
K = mean(sampNUTSgompertz[:K])
sigma = mean(sampNUTSgompertz[:sigma])
tau = mean(sampNUTSgompertz[:tau])
ess = ess(describe(sampNUTSgompertz)[1])
