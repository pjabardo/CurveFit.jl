module CurveFit

# package code goes here
abstract Approx
abstract LeastSquares <: Approx

include("linfit.jl")
include("nonlinfit.jl")
include("king.jl")

end # module
