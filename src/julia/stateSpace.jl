using CSV
using DataFrames
using Turing
using Distributions


#Gompertz model
@model function gompertz(y, pr_type)
    if pr_type == 0
        σ ~ Uniform(0,1)
        τ ~ Uniform(0,1)
        α ~ Uniform(0,1)
        β ~ Uniform(0,1)
    else
        σ ~ InverseGamma(0.001, 0.001)
        τ ~ InverseGamma(0.001, 0.001)
        α ~ Normal(0,100)
        β ~ Normal(0,100)
    end
    u_init ~ Uniform(0,1) #not the ideal way to specify an improper prior?
    n = size(y)[1]
    umed = Vector(undef, n)
    u = Vector(undef, n)
    umed[1] = u_init
    u[1] ~ Normal(umed[1], σ)
    for t in 2:n
        umed[t] = α + β*u[t-1]
        u[t] ~ Normal(umed[t], σ)
    end
    for t in 1:n
        y[t] ~ Normal(u[t], τ)
    end
end


#Logistic growth model
@model function logistic(y)
    if pr_type == 0
        σ ~ Uniform(0,1)
        τ ~ Uniform(0,1)
        α ~ Uniform(0,1)
        β ~ Uniform(0,1)
    else
        σ ~ InverseGamma(0.001, 0.001)
        τ ~ InverseGamma(0.001, 0.001)
        α ~ Normal(0,100)
        β ~ Normal(0,100)
    end
    u_init ~ Uniform(0,1) #not the ideal way to specify an improper prior?
    n = size(y)[1]
    umed = Vector(undef, n)
    u = Vector(undef, n)
    umed[1] = u_init
    u[1] ~ LogNormal(umed[1], σ)
    for t in 2:n
        umed[t] = log(u[t-1] + r*u[t-1]*(1-u[t-1]/K))
        u[t] ~ LogNormal(umed[t], σ)
    end
    y ~ LogNormal(log(u), τ)
end

#read in csv file: requires CSV and DataFrames
df = CSV.read("C:\\Users\\Andrea.Havron\\Documents\\github\\FIMSarch\\data\\gompertz.csv", DataFrame)
# . syntax extracts column from dataframe and casts as vector
gompertzDat = df.y
df = CSV.read("C:\\Users\\Andrea.Havron\\Documents\\github\\FIMSarch\\data\\logistic.csv", DataFrame)
# . syntax extracts column from dataframe and casts as vector
logisticDat = df.y


sampHMC = sample(gompertz(gompertzDat,1), NUTS(0.65), 1000)
describe(sampHMC)

sampNUTSlogistic = sample(logistic(logisticDat), NUTS(0.65), 1000)
