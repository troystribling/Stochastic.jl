__precompile__(true)

module Stochastic

using Distributions
using Gadfly

import Distributions: rand
import Gadfly: plot

export

  # random process types
  RandomProcess,
  BrownianMotion,

  # methods
  params

  # source files
  include("Plots.jl")
  include("BrownianMotion.jl")

end
