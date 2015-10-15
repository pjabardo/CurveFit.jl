module CurveFit

export least_squares, linear_fit, log_fit, power_fit, exp_fit, poly_fit
export LinearFit, PolyFit, LogFit, PowerFit, ExpFit
export curve_fit, apply_fit
export nonlinear_fit, makeDerivFun
export linear_king_fit, king_fit
export LinearKingFit, KingFit
export fit, coef

# package code goes here
abstract ApproxFit
abstract LeastSquares <: ApproxFit


include("linfit.jl")
include("rationalfit.jl")
include("nonlinfit.jl")
include("king.jl")

end # module
