abstract  RandomProcess

doc"""
  SampleInterval(tmin, tmax, npts)

Time interval for evaluluation of randowm process where tmin is the start time, tmax the end time
and npts the number of points in the sample.
"""
immutable SampleInterval{T<:Integer, U<:Real}
  npts::T
  tmax::U
  tmin::U
  SampleInterval(npts::T, tmax::U, tmin::U) = tmax > tmin && npts > U(1) ? new(npts, tmax, tmin) : error("must have tmax > tmin && npts > 0")
end

# outer constructors
SampleInterval{T<:Integer, U<:Real}(npts::T, tmax::U, tmin::U) = SampleInterval{T, U}(npts, tmax, tmin)
SampleInterval(npts::Integer, tmax::Integer, tmin::Integer) = SampleInterval(UInt64(npts), Float64(tmax), Float64(tmin))
SampleInterval(npts::Real, tmax::Real, tmin::Real) = SampleInterval(UInt64(npts), promote(tmax, tmin)...)
SampleInterval(npts::Integer, tmax::Real, tmin::Real) = SampleInterval(npts, promote(tmax, tmin)...)
SampleInterval(npts::Integer, tmax::Real) = SampleInterval(npts, tmax, zero(tmax))

# conversions
convert{T<:Integer, U<:Real, V<:Integer, W<:Real}(::Type{SampleInterval{T, U}}, npts::V, tmax::W, tmin::W) = SampleInterval(T(npts), U(tmax), U(tmin))
convert{T<:Integer, U<:Real, V<:Integer, W<:Real}(::Type{SampleInterval{T, U}}, s::SampleInterval{V, W}) = SampleInterval(T(s.npts), U(s.tmax), U(s.tmin))

# parameters
params(sampleInterval::SampleInterval) = (sampleInterval.npts, sampleInterval.tmax, sampleInterval.tmin)
Δt(sampleInterval::SampleInterval) =  (sampleInterval.tmax - sampleInterval.tmin) / sampleInterval.npts

function samples(sampleInterval::SampleInterval)
    collect(linspace(sampleInterval.tmin, sampleInterval.tmax, sampleInterval.npts))
end

doc"""
BrownianMotion(μ, σ)

Brownian motion with mean μ and standard deviation σ.

https://en.wikipedia.org/wiki/Brownian_motion
"""
immutable BrownianMotion{T<:Real} <: RandomProcess
  μ::T
  σ::T
  BrownianMotion(μ, σ) = σ > zero(σ) ? new(μ, σ) : error("must have σ > 0")
end

# outer constructors
BrownianMotion{T<:Real}(μ::T, σ::T) = BrownianMotion{T}(μ, σ)
BrownianMotion(μ::Integer, σ::Integer) = BrownianMotion(Float64(μ), Float64(σ))
BrownianMotion(μ::Real, σ::Real) = BrownianMotion(promote(μ, σ)...)
BrownianMotion(μ::Real) = BrownianMotion(μ, 1.0)
BrownianMotion() = BrownianMotion(0.0, 1.0)

# conversions
convert{T<:Real, U<:Real}(::Type{BrownianMotion{T}}, μ::U, σ::U) = BrownianMotion(T(μ), T(σ))
convert{T<:Real, U<:Real}(::Type{BrownianMotion{T}}, bm::BrownianMotion{U}) = BrownianMotion(T(bm.μ), T(bm.σ))

# parameters
params(randomProcess::BrownianMotion) = (randomProcess.μ, randomProcess.σ)

# generation
function increments(randomProcess::BrownianMotion, interval::SampleInterval)
  μ = Δt(interval) * randomProcess.μ
  σ = sqrt(Δt(interval)) * randomProcess.σ
  bm = zeros(Float64, interval.npts-1)
  for n = 1:interval.npts-1
    bm[n] = rand(Normal(μ, σ))
  end
  return bm
end

function rand(randomProcess::BrownianMotion, interval::SampleInterval)
  bm = zeros(Float64, interval.npts)
  Δbm = increments(randomProcess, interval)
  for n = 2:interval.npts
    bm[n] = bm[n-1] + Δbm[n-1]
  end
  return bm
end

doc"""
GeometricBrownianMotion(μ, σ)

geometric Brownian motion with mean μ and standard deviation σ.

https://en.wikipedia.org/wiki/Geometric_Brownian_motion
"""

immutable GeometricBrownianMotion{T<:Real} <: RandomProcess
  brownianMotion::BrownianMotion{T}
  s0::T

  GeometricBrownianMotion(brownianMotion, s0) = s0 != zero(s0) ? new{T}(brownianMotion, s0) : error("must have s0 ! 0")
end

# outer constructors
GeometricBrownianMotion{T<:Real}(μ::T, σ::T, s0::T) = GeometricBrownianMotion{T}(BrownianMotion(μ, σ), s0)
GeometricBrownianMotion(μ::Integer, σ::Integer, s0::Integer) = GeometricBrownianMotion(Float64(μ), Float64(σ), Float64(s0))
GeometricBrownianMotion(μ::Real, σ::Real, s0::Real) = GeometricBrownianMotion(promote(μ, σ, s0)...)
GeometricBrownianMotion(μ::Real, σ::Real) = GeometricBrownianMotion(Float64(μ), Float64(σ), 1.0)
GeometricBrownianMotion(μ::Real) = GeometricBrownianMotion(Float64(μ), 1.0, 1.0)
GeometricBrownianMotion() = GeometricBrownianMotion(0.0, 1.0, 1.0)

# conversions
convert{T<:Real, U<:Real}(::Type{GeometricBrownianMotion{T}}, μ::U, σ::U, s0::T) = GeometricBrownianMotion(T(μ), T(σ, t(s0)))
convert{T<:Real, U<:Real}(::Type{GeometricBrownianMotion{T}}, gbm::GeometricBrownianMotion{U}) = GeometricBrownianMotion(T(gbm.μ), T(gbm.σ), T(gbm.s0))

# parameters
params(randomProcess::GeometricBrownianMotion) = (randomProcess.brownianMotion.μ, randomProcess.brownianMotion.σ, randomProcess.s0)

# generation
function increments(randomProcess::GeometricBrownianMotion, interval::SampleInterval)
  Δbm = increments(randomProcess.brownianMotion, interval)
  map(x -> exp(x), Δbm)
end

function rand(randomProcess::GeometricBrownianMotion, interval::SampleInterval)
  gbm = ones(Float64, interval.npts)
  Δgbm = increments(randomProcess, interval)
  gbm[1] = randomProcess.s0
  for n = 1:interval.npts-1
      gbm[n+1] = gbm[n] * Δgbm[n]
  end
  return gbm
end
