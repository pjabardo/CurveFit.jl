using LinearAlgebra
sq(x) = x*x

"Fits a straight line through a set of points, `y = a₁ + a₂ * x`"
function linear_fit(x, y)

    
    sx = sum(x)
    sy = sum(y)

    m = length(x)

    sx2 = zero(sx)
    sy2 = zero(sy)
    sxy = zero(sx*sy)

    for i = 1:m
        sx2 += x[i]*x[i]
        sy2 += y[i]*y[i]
        sxy += x[i]*y[i]
    end

    a0 = (sx2*sy - sxy*sx) / ( m*sx2 - sx*sx )
    a1 = (m*sxy - sx*sy) / (m*sx2 - sx*sx)

    return [a0, a1]
end

"Fits a log function through a set of points: `y = a₁+ a₂*log(x)`"
log_fit(x, y) = linear_fit(log.(x), y)

"Fits a power law through a set of points: `y = a₁*x^a₂`"
function power_fit(x, y)
    fit = linear_fit(log.(x), log.(y))
    [exp(fit[1]), fit[2]]
end

"Fits an `exp` through a set of points: `y = a₁*exp(a₂*x)`"
function exp_fit(x, y)
    fit = linear_fit(x, log.(y))
    [exp(fit[1]), fit[2]]
end

"""
Create Vandermonde matrix for simple polynomial fit
"""
function vandermondepoly(x, y, n)
    m = length(x)
    A = zeros(eltype(y), m, n+1)
    A[:,1] .= 1.0

    for i = 1:n
        for k in 1:m
            A[k, i+1] = A[k, i] * x[k]
        end
    end
    return A
end

"""
    Calculates the coefficients of a linear model using Least Squares

Given a Vandermonde matrix, this function uses the QR factorization to compute the Least Squares
fit of a linear model. The code probably could be optimized.

"""
function fit_linear_model(Av, y)
    m, n = size(Av)
    F = qr(Av)

    d = (F.Q * Matrix(I, m, n))' * y
    R = UpperTriangular(F.R)

    coefs = R\d

    return coefs
   
end

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
struct LinearFit{T<:Number} <: LeastSquares
    coefs::Array{T,1}
    LinearFit{T}(coefs) where {T<:Number} = new(copy(coefs))
end
LinearFit(x::AbstractVector{T}, y::AbstractVector{T}) where {T<:Number} = LinearFit{T}(linear_fit(x, y))


"High Level interface for fitting log laws"
struct LogFit{T<:Number} <: LeastSquares
    coefs::Array{T,1}
    LogFit{T}(coefs) where {T<:Number} = new(copy(coefs))
end
LogFit(x::AbstractVector{T}, y::AbstractVector{T}) where {T<:Number} = LogFit{T}(log_fit(x, y))

"High Level interface for fitting power laws"
struct PowerFit{T<:Number} <: LeastSquares
    coefs::Array{T,1}
    PowerFit{T}(coefs) where {T<:Number} = new(copy(coefs))
end
PowerFit(x::AbstractVector{T}, y::AbstractVector{T}) where {T<:Number} = PowerFit{T}(power_fit(x, y))

"High Level interface for fitting exp laws"
struct ExpFit{T<:Number} <: LeastSquares
    coefs::Array{T,1}
    ExpFit{T}(coefs) where {T<:Number} = new(copy(coefs))
end
ExpFit(x::AbstractVector{T}, y::AbstractVector{T}) where {T<:Number} = ExpFit{T}(exp_fit(x, y))



"""
# Generic interface for curve fitting.

The same function `curve_fit` can be used to fit the data depending on fit type, 
shich is specified in the first parameter. This function returns an object that
can be used to estimate the value of the fitting model using function `apply_fit`.

## A few examples:

 * `f = curve_fit(LinearFit, x, y)`
 * `f = curve_fit(Poly, x, y, n)`
"""
curve_fit(::Type{T}, x, y) where {T<:LeastSquares} = T(x, y)
curve_fit(::Type{T}, x, y, args...) where {T<:LeastSquares} = T(x, y, args...)
curve_fit(::Type{Poly}, x, y, n=1) = Poly(poly_fit(x, y, n))



(f::LinearFit)(x) = f.coefs[1] + f.coefs[2] *x
(f::PowerFit)(x) = f.coefs[1] * x ^ f.coefs[2]
(f::LogFit)(x) = f.coefs[1] + f.coefs[2] * log(x)
(f::ExpFit)(x) = f.coefs[1] * exp(f.coefs[2] * x)
