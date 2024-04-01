
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

  - Straight lines
  - Polynomial fitting
  - Power laws
  - Log laws
  - Exp laws
  - Rational polynomial fitting
  - A generic non-linear fitting algorithm
  - King's law (used in hotwire anemometry)
"""
module CurveFit

using Markdown: @doc_str
using Polynomials

export linear_fit, log_fit, power_fit, exp_fit, poly_fit
export LinearFit, LogFit, PowerFit, ExpFit
export curve_fit, apply_fit
export nonlinear_fit, gauss_newton_fit, gauss_newton_generic_fit
export linear_king_fit, king_fit
export LinearKingFit, KingFit
export fit, coef
export linear_rational_fit, RationalPoly, ratval, rational_fit
export Polynomial
export ExpSumFit, expsum_fit, expsum_init

#"Abstract base class for fitting data"
abstract type AbstractApproxFit end

#"Abstract class for least squares fitting of data"
abstract type AbstractLeastSquares <: AbstractApproxFit end

include("linfit.jl")
include("rationalfit.jl")
include("nonlinfit.jl")
include("king.jl")
include("expsumfit.jl")

"""
# Uses the object created by `curve_fit` to estimate values

The `call` method is overloaded so that the fit object can
be used as a function:

## Example:

```julia
x = 1.0:10.0
@. y = 2*x + 1 + randn()

fit = curve_fit(LinearFit, x, y)

y1 = fit(5.1)
y2 = apply_fit(fit, 5.1)
```
"""
apply_fit(f::T, x) where {T <: AbstractLeastSquares} = f(x)

end # module
