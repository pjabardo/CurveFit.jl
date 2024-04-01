# Rational polynomial interpolation

"""
Linear Rational LeastSquares

The following curvit is done:

`y = p(x) / q(x)`

where `p(x)` and `q(x)` are polynomials.

The linear case is solved by doing a least square fit on

`y*q(x) = p(x)`

where the zero order term o `q(x)` is assumed to be 1.
"""
function linear_rational_fit(
        x::AbstractVector{T}, y::AbstractVector{T}, p, q) where {T <: Number}
    n = size(x, 1)
    A = zeros(T, n, q + p + 1)
    for i in 1:n
        A[i, 1] = one(T)
        for k in 1:p
            A[i, k + 1] = x[i]^k
        end
        for k in 1:q
            A[i, p + 1 + k] = -y[i] * x[i]^k
        end
    end

    fit_linear_model(A, y)
    #    coefs[1:p+1], [1.0; coefs[p+2:end]]

end

"""
# Type defining a rational polynomial

A rational polynomial is the ratio of two polynomials
and it is often useful in approximating functions.
"""
struct RationalPoly{T <: Number} <: AbstractLeastSquares
    num::Vector{T}
    den::Vector{T}
end
function RationalPoly(p::Integer, q::Integer, ::Type{T} = Float64) where {T <: Number}
    RationalPoly(zeros(T, p + 1), zeros(T, q + 1))
end

function RationalPoly(coefs::AbstractVector{T}, p, q) where {T <: Number}
    RationalPoly(collect(coefs[1:(p + 1)]), [one(T); collect(coefs[(p + 2):end])])
end

"""
Evaluate a rational polynomial
"""
ratval(r::RationalPoly, x) = evalpoly(x, r.num) / evalpoly(x, r.den)

"""
`call` overload for calling directly `ratval`
"""
(r::RationalPoly)(x) = ratval(r, x)

"""
Auxiliary function used in nonlinear least squares
"""
function make_rat_fun(p, q)
    r = RationalPoly(p, q, Float64)

    function (x, a)
        for i in 0:p
            r.num[i + 1] = a[i + 1]
        end
        r.den[1] = 1
        for i in 1:q
            r.den[i + 1] = a[p + 1 + i]
        end
        evalpoly(x[1], r.num) / evalpoly(x[1], r.den) - x[2]
    end
end

"""
# Carry out a nonlinear least squares of rational polynomials

Find the polynomial coefficients that best approximate
the points given by `x` and `y`.
"""
function rational_fit(x, y, p, q, eps = 1e-8, maxiter = 200)
    coefs0 = linear_rational_fit(x, y, p, q)

    fun = make_rat_fun(p, q)

    coefs, converged, niter = nonlinear_fit(hcat(x, y), fun, coefs0, eps, maxiter)

    coefs
end

function curve_fit(::Type{RationalPoly}, x, y, p, q, eps = 1e-8, maxiter = 200)
    RationalPoly(rational_fit(x, y, p, q, eps, maxiter), p, q)
end
