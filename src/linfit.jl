
using Polynomials

least_squares{T<:Number}(A::StridedVecOrMat{T}, y::Vector{T}) = A \ y


linear_fit(x, y) = least_squares(hcat(ones(x), x), y)


log_fit(x, y) = linear_fit(log(x), y)

function power_fit(x, y)
    fit = linear_fit(log(x), log(y))
    [exp(fit[1]), fit[2]]
end
function exp_fit(x, y)
    fit = linear_fit(x, log(y))
    [exp(fit[1]), fit[2]]
end

function poly_fit(x, y, n)

    nx = length(x)
    A = zeros(Float64, nx, n+1)
    A[:,1] = 1.0
    for i in 1:n
        for k in 1:nx
            A[k,i+1] = A[k,i] * x[k]
        end
    end
    Poly(least_squares(A,y))
end


    



type LinearFit{T} <: LeastSquares
    coefs::Array{T,1}
    LinearFit(coefs) = new(copy(coefs))
end
LinearFit{T}(x::AbstractVector{T}, y::AbstractVector{T}) = LinearFit{T}(linear_fit(x, y))



type LogFit{T} <: LeastSquares
    coefs::Array{T,1}
    LogFit(coefs) = new(copy(coefs))
end
LogFit{T}(x::AbstractVector{T}, y::AbstractVector{T}) = LogFit{T}(log_fit(x, y))

type PowerFit{T} <: LeastSquares
    coefs::Array{T,1}
    PowerFit(coefs) = new(copy(coefs))
end
PowerFit{T}(x::AbstractVector{T}, y::AbstractVector{T}) = PowerFit{T}(power_fit(x, y))

type ExpFit{T} <: LeastSquares
    coefs::Array{T,1}
    ExpFit(coefs) = new(copy(coefs))
end
ExpFit{T}(x::AbstractVector{T}, y::AbstractVector{T}) = ExpFit{T}(exp_fit(x, y))




curve_fit{T<:LeastSquares}(::Type{T}, x, y) = T(x, y)
curve_fit{T<:LeastSquares}(::Type{T}, x, y, args...) = T(x, y, args...)
curve_fit(::Type{Poly}, x, y, n) = poly_fit(x, y, n)
curve_fit(::Type{Poly}, x, y) = poly_fit(x, y, 1)


apply_fit(f::LinearFit, x) = f.coefs[1] .+ f.coefs[2].*x
apply_fit(f::PowerFit, x) = f.coefs[1] .* x .^ f.coefs[2]
apply_fit(f::LogFit, x) = f.coefs[1] .+ f.coefs[2] .* log(x)
apply_fit(f::ExpFit, x) = f.coefs[1] .* exp(f.coefs[2] .* x)

apply_fit(f::Poly, x) = polyval(f, x)

import Base.call

call{T<:LeastSquares}(f::T, x) = apply_fit(f, x)


