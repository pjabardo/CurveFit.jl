



"Fits a straight line through a set of points, `y = a₁ + a₂ * x`"
linear_fit(x, y) = hcat(ones(x), x) \ y

"Fits a log function through a set of points: `y = a₁+ a₂*log(x)`"
log_fit(x, y) = linear_fit(log(x), y)

"Fits a power law through a set of points: `y = a₁*x^a₂`"
function power_fit(x, y)
    fit = linear_fit(log(x), log(y))
    [exp(fit[1]), fit[2]]
end

"Fits an `exp` through a set of points: `y = a₁*exp(a₂*x)`"
function exp_fit(x, y)
    fit = linear_fit(x, log(y))
    [exp(fit[1]), fit[2]]
end

"""
Fits a polynomial of degree `n` through a set of points.

Simple algorithm, doesn't use orthogonal polynomials or any such thing 
and therefore unconditioned matrices are possible. Use it only for low
degree polynomials. 

This function returns a the coefficients of the polynomial.
"""
function poly_fit(x, y, n)

    nx = length(x)
    A = zeros(eltype(x), nx, n+1)
    A[:,1] = 1.0
    for i in 1:n
        for k in 1:nx
            A[k,i+1] = A[k,i] * x[k]
        end
    end
    A\y
end


    


"""
High Level interface for fitting straight lines
"""
immutable LinearFit{T<:Number} <: LeastSquares
    coefs::Array{T,1}
    LinearFit(coefs) = new(copy(coefs))
end
LinearFit{T<:Number}(x::AbstractVector{T}, y::AbstractVector{T}) = LinearFit{T}(linear_fit(x, y))


"High Level interface for fitting log laws"
immutable LogFit{T<:Number} <: LeastSquares
    coefs::Array{T,1}
    LogFit(coefs) = new(copy(coefs))
end
LogFit{T<:Number}(x::AbstractVector{T}, y::AbstractVector{T}) = LogFit{T}(log_fit(x, y))

"High Level interface for fitting power laws"
immutable PowerFit{T<:Number} <: LeastSquares
    coefs::Array{T,1}
    PowerFit(coefs) = new(copy(coefs))
end
PowerFit{T<:Number}(x::AbstractVector{T}, y::AbstractVector{T}) = PowerFit{T}(power_fit(x, y))

"High Level interface for fitting exp laws"
immutable ExpFit{T<:Number} <: LeastSquares
    coefs::Array{T,1}
    ExpFit(coefs) = new(copy(coefs))
end
ExpFit{T<:Number}(x::AbstractVector{T}, y::AbstractVector{T}) = ExpFit{T}(exp_fit(x, y))



"""
# Generic interface for curve fitting.

The same function `curve_fit` can be used to fit the data depending on fit type, 
shich is specified in the first parameter. This function returns an object that
can be used to estimate the value of the fitting model using function `apply_fit`.

## A few examples:

 * `f = curve_fit(LinearFit, x, y)`
 * `f = curve_fit(Poly, x, y, n)`
"""
curve_fit{T<:LeastSquares}(::Type{T}, x, y) = T(x, y)
curve_fit{T<:LeastSquares}(::Type{T}, x, y, args...) = T(x, y, args...)
curve_fit(::Type{Poly}, x, y, n=1) = Poly(poly_fit(x, y, n))


"""
#Uses the object created by `curve_fit` to estimate values

The `call` method is overloaded so that the fit object can 
be used as a function:

## Example:
x = [linspace(1, 10, 10);]
y = 2*x + 1 + randn(10)

fit = curve_fit(LinearFit, x, y)

y1 = fit(5.1)
y2 = apply_fit(fit, 5.1)
"""
apply_fit(f::LinearFit, x) = f.coefs[1] .+ f.coefs[2].*x
apply_fit(f::PowerFit, x) = f.coefs[1] .* x .^ f.coefs[2]
apply_fit(f::LogFit, x) = f.coefs[1] .+ f.coefs[2] .* log(x)
apply_fit(f::ExpFit, x) = f.coefs[1] .* exp(f.coefs[2] .* x)

apply_fit(f::Poly, x) = polyval(f, x)

import Base.call

call{T<:LeastSquares}(f::T, x) = apply_fit(f, x)


