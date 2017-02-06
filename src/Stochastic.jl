__precompile__(true)

module Stochastic

using Distributions

import Distributions: rand

export

  # random process types
  RandomProcess,
  BrownianMotion,

  # methods
  params

  # source files
  include("Utils.jl")
  include("RandomProcess.jl")

end
