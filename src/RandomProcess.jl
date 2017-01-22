abstract  RandomProcess

immutable BrownianMotion{T <: Real} <: RandomProcess
  μ::T
  σ::T
  BrownianMotion(μ, σ) = σ < 0 ? error("σ < 0") : new(μ, σ)
end

BrownianMotion() = BrownianMotion(0.0, 1.0)
BrownianMotion(μ::Real) = BrownianMotion(μ, 1.0)
BrownianMotion(μ::Integer, σ::Integer) = BrownianMotion(Float64(μ), Float64(σ))
BrownianMotion(μ::Real, σ::Real) = BrownianMotion(promote(μ, σ)...)
BrownianMotion{T<:Real}(μ::T, σ::T) = BrownianMotion{T}(μ, σ)
