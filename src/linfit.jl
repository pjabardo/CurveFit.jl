using LinearAlgebra

sq(x) = x * x

"""
Fits a straight line through a set of points, `y = a₁ + a₂ * x`
"""
function linear_fit(x, y)
    sx = sum(x)
    sy = sum(y)

    m = length(x)

    sx2 = zero(sx .* sx)
    sy2 = zero(sy .* sy)
    sxy = zero(sx * sy)

    @inbounds for i in eachindex(x)
        sx2 += x[i] * x[i]
        sy2 += y[i] * y[i]
        sxy += x[i] * y[i]
    end

    a0 = (sx2 * sy - sxy * sx) / (m * sx2 - sx * sx)
    a1 = (m * sxy - sx * sy) / (m * sx2 - sx * sx)

    return (a0, a1)
end

"""
Fits a log function through a set of points: `y = a₁+ a₂*log(x)`
"""
log_fit(x, y) = linear_fit(log.(x), y)

"""
Fits a power law through a set of points: `y = a₁*x^a₂`
"""
function power_fit(x, y)
    fit = linear_fit(log.(x), log.(y))
    return (exp(fit[1]), fit[2])
end

"""
Fits an `exp` through a set of points: `y = a₁*exp(a₂*x)`
"""
function exp_fit(x, y)
    fit = linear_fit(x, log.(y))
    return (exp(fit[1]), fit[2])
end

"""
Create Vandermonde matrix for simple polynomial fit
"""
function vandermondepoly(x, y, n)
    m = length(x)
    A = zeros(eltype(y), m, n + 1)
    A[:, 1] .= 1.0

    @inbounds for i in 1:n, k in 1:m
        A[k, i + 1] = A[k, i] * x[k]
    end

    return A
end

"""
Calculates the coefficients of a linear model using Least Squares

Given a Vandermonde matrix, this function uses the QR factorization to compute the Least
Squares fit of a linear model.
"""
fit_linear_model(Av, y) = qr(Av) \ y

"""
Fits a polynomial of degree `n` through a set of points.

Simple algorithm, doesn't use orthogonal polynomials or any such thing
and therefore unconditioned matrices are possible. Use it only for low
degree polynomials.

This function returns a the coefficients of the polynomial.
"""
function poly_fit(x, y, n)
    A = vandermondepoly(x, y, n)
    coefs = fit_linear_model(A, y)
    return coefs
end

"""
High Level interface for fitting straight lines
"""
struct LinearFit{T <: Number} <: AbstractLeastSquares
    coefs::NTuple{2, T}

    LinearFit{T}(coefs) where {T <: Number} = new((coefs[1], coefs[2]))
    LinearFit{T}(c1, c2) where {T <: Number} = new((c1, c2))
end
function LinearFit(x::AbstractVector{T}, y::AbstractVector{T}) where {T <: Number}
    return LinearFit{T}(linear_fit(x, y))
end

"""
High Level interface for fitting log laws
"""
struct LogFit{T <: Number} <: AbstractLeastSquares
    coefs::NTuple{2, T}

    LogFit{T}(coefs) where {T <: Number} = new((coefs[1], coefs[2]))
    LogFit{T}(c1, c2) where {T <: Number} = new((c1, c2))
    #LogFit{T}(coefs) where {T<:Number} = new(copy(coefs))
end
function LogFit(x::AbstractVector{T}, y::AbstractVector{T}) where {T <: Number}
    return LogFit{T}(log_fit(x, y))
end

"""
High Level interface for fitting power laws
"""
struct PowerFit{T <: Number} <: AbstractLeastSquares
    coefs::NTuple{2, T}

    PowerFit{T}(coefs) where {T <: Number} = new((coefs[1], coefs[2]))
    PowerFit{T}(c1, c2) where {T <: Number} = new((c1, c2))
    #PowerFit{T}(coefs) where {T<:Number} = new(copy(coefs))
end
function PowerFit(x::AbstractVector{T}, y::AbstractVector{T}) where {T <: Number}
    return PowerFit{T}(power_fit(x, y))
end

"""
High Level interface for fitting exp laws
"""
struct ExpFit{T <: Number} <: AbstractLeastSquares
    coefs::NTuple{2, T}

    ExpFit{T}(coefs) where {T <: Number} = new((coefs[1], coefs[2]))
    ExpFit{T}(c1, c2) where {T <: Number} = new((c1, c2))
    #ExpFit{T}(coefs) where {T<:Number} = new(copy(coefs))
end
function ExpFit(x::AbstractVector{T}, y::AbstractVector{T}) where {T <: Number}
    return ExpFit{T}(exp_fit(x, y))
end

"""
# Generic interface for curve fitting.

The same function `curve_fit` can be used to fit the data depending on fit type,
shich is specified in the first parameter. This function returns an object that
can be used to estimate the value of the fitting model using function `apply_fit`.

## A few examples:

  - `f = curve_fit(LinearFit, x, y)`
  - `f = curve_fit(Polynomial, x, y, n)`
"""
curve_fit(::Type{T}, x, y) where {T <: AbstractLeastSquares} = T(x, y)
curve_fit(::Type{T}, x, y, args...) where {T <: AbstractLeastSquares} = T(x, y, args...)
curve_fit(::Type{Polynomial}, x, y, n = 1) = Polynomial(poly_fit(x, y, n))

(f::LinearFit)(x) = f.coefs[1] + f.coefs[2] * x
(f::PowerFit)(x) = f.coefs[1] * x^f.coefs[2]
(f::LogFit)(x) = f.coefs[1] + f.coefs[2] * log(x)
(f::ExpFit)(x) = f.coefs[1] * exp(f.coefs[2] * x)
