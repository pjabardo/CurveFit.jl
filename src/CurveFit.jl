
"""
# Simple Least Squares fitting

The `CurveFit` module provides functions that
implement a few least squares approximations.

It is simple in an engineering sense. It simply
returns the coefficients and does not do 
any error analysis. 

It does, however, provide a simple and common interface
to the routines.

The package also includes nonlinear least squares fitting
using a Newton type algorithm

The fitting algorithms include

 * Straight lines
 * Polynomial fitting
 * Power laws
 * log laws
 * exp laws
 * Rational polynopmial fitting
 * A generic non-linear fitting algorithm
 * King's law (used in hotwire anemometry)

"""
module CurveFit
using Polynomials

export linear_fit, log_fit, power_fit, exp_fit, poly_fit
export LinearFit, PolyFit, LogFit, PowerFit, ExpFit
export curve_fit, apply_fit
export nonlinear_fit
export linear_king_fit, king_fit
export LinearKingFit, KingFit
export fit, coef
export linear_rational_fit, RationalPoly, ratval, rational_fit
export Poly

# package code goes here
"Abstract base class for fitting data"
abstract ApproxFit

"Abstract class for least squares fitting of data"
abstract LeastSquares <: ApproxFit


include("linfit.jl")
include("rationalfit.jl")
include("nonlinfit.jl")
include("king.jl")

end # module
