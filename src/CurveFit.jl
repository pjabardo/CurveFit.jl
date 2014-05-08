module CurveFit

export least_squares, linear_fit, log_fit, power_fit, exp_fit, poly_fit
export LinearFit, PolyFit, LogFit, PowerFit, ExpFit
export curve_fit, apply_fit
export nonlinear_fit, makeDerivFun
export linear_king_fit, king_fit
export LinearKingFit, KingFit

# package code goes here
abstract Approx
abstract LeastSquares <: Approx

include("linfit.jl")
include("nonlinfit.jl")
include("king.jl")

end # module
