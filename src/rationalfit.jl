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
function linear_rational_fit(x::AbstractVector{T}, y::AbstractVector{T}, p, q) where T<:Number
    n = size(x,1)
    A = zeros(T, n, q+p+1)
    for i = 1:n
        A[i,1] = one(T)
        for k = 1:p
            A[i,k+1] = x[i]^k
        end
        for k = 1:q
            A[i, p+1+k] = -y[i] * x[i]^k
        end
    end

    A \ y
#    coefs[1:p+1], [1.0; coefs[p+2:end]]

end

"""
# Type defining a rational polynomial

A rational polynomial is the ratio of two polynomials
and it is often useful in approximating functions.
"""
struct RationalPoly{T<:Number} <: LeastSquares
    num::Poly{T}
    den::Poly{T}
end
RationalPoly(a::AbstractVector{T}, b::AbstractVector{T}) where {T<:Number} = RationalPoly(Poly(a), Poly(b))
RationalPoly(p::Integer, q::Integer, ::Type{T}=Float64) where {T<:Number} = RationalPoly{T}(Poly(zeros(T,p+1)), Poly(zeros(T,q+1)))
RationalPoly(coefs::AbstractVector{T}, p, q) where {T<:Number} = RationalPoly(coefs[1:p+1],[1.0; coefs[p+2:end]])

"Evaluate a rational polynomial"
ratval(r::RationalPoly{T}, x) where {T<:Number} = polyval(r.num, x) ./ polyval(r.den, x)

"`call` overload for calling directly `ratval`"
(r::RationalPoly)(x) = ratval(r, x)

"Auxiliary function used in nonlinear least squares"
function make_rat_fun(p, q)
    r = RationalPoly(p, q, Float64)
    
    function(x, a)
        for i=0:p
            r.num[i] = a[i+1]
        end
        r.den[0] = 1
        for i = 1:q
            r.den[i] = a[p+1+i]
        end
        polyval(r.num, x[1]) / polyval(r.den, x[1]) - x[2]
    end
    
end

"""
# Carry out a nonlinear least squares of rational polynomials

Find the polynomial coefficients that best approximate
the points given by `x` and `y`.
"""
function rational_fit(x, y, p, q, eps=1e-8, maxiter=200)

    coefs0  = linear_rational_fit(x, y, p, q)

    fun = make_rat_fun(p, q)

    coefs, converged, niter = nonlinear_fit(hcat(x, y), fun, coefs0, eps, maxiter)
    
    coefs

end


function curve_fit(::Type{RationalPoly}, x, y, p, q, eps=1e-8, maxiter=200)
    RationalPoly(rational_fit(x, y, p, q, eps, maxiter), p, q)
end
    
    
