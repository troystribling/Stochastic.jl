__precompile__(true)

module Stochastic

using Distributions

export
  RandomProcess,
  BrownianMotion

  include("Utils.jl")
  include("RandomProcess.jl")

end
