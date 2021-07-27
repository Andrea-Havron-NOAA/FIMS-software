using Turing
using StatsPlots

@model function gdemo(x, y)
  # Assumptions
  σ ~ InverseGamma(2,3)
  μ ~ Normal(0,sqrt(σ))
  # Observations
  x ~ Normal(μ, sqrt(σ))
  y ~ Normal(μ, sqrt(σ))
end
