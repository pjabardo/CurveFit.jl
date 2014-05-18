    
least_squares{T<:Number}(A::StridedVecOrMat{T}, y::Vector{T}) = A \ y


linear_fit(x, y) = least_squares(hcat(ones(length(x)), x), y)


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
    least_squares(A,y)
end


    



type LinearFit <: LeastSquares
    coefs::Array{Float64,1}
    LinearFit(coefs) = new(copy(coefs))
end
LinearFit(x, y) = LinearFit(linear_fit(x, y))

type PolyFit <: LeastSquares
    coefs::Array{Float64,1}
    degree::Int
    PolyFit(coefs) = new(copy(coefs), length(coefs)-1)
end

PolyFit(x, y, n) = PolyFit(poly_fit(x, y, n))

type LogFit <: LeastSquares
    coefs::Array{Float64,1}
    LogFit(coefs) = new(copy(coefs))
end
LogFit(x, y) = LogFit(log_fit(x, y))

type PowerFit <: LeastSquares
    coefs::Array{Float64,1}
    PowerFit(coefs) = new(copy(coefs))
end
PowerFit(x, y) = PowerFit(power_fit(x, y))

type ExpFit <: LeastSquares
    coefs::Array{Float64,1}
    ExpFit(coefs) = new(copy(coefs))
end
ExpFit(x, y) = ExpFit(exp_fit(x, y))




curve_fit{T}(::Type{T}, x, y) = T(x, y)
curve_fit{T}(::Type{T}, x, y, args...) = T(x, y, args...)

fit{T<:LeastSquares}(::Type{T}, x, y) = T(x, y)
fit{T<:LeastSquares}(::Type{T}, x, y, args...) = T(x, y, args...)
coef{T<:LeastSquares}(f::T) = f.coefs

apply_fit(f::LinearFit, x) = f.coefs[1] .+ f.coefs[2].*x
apply_fit(f::PowerFit, x) = f.coefs[1] .* x .^ f.coefs[2]
apply_fit(f::LogFit, x) = f.coefs[1] .+ f.coefs[2] .* log(x)
apply_fit(f::ExpFit, x) = f.coefs[1] .* exp(f.coefs[2] .* x)


function apply_fit(f::PolyFit, x)
    a = f.coefs
    n = f.degree+1
    y = a[n] .+ 0.0.*x
    for i = (n-1):-1:1
        y = (a[i] .+ x .* y)
    end
    return y
end

