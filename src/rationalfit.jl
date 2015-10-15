# Rational polynomial interpolation


function linear_rational_fit{T<:Number}(x::AbstractArray{T}, y::AbstractVector{T}, p, q)
    A = zeros(T, q+p+1)
    n = size(x,1)
    
    for i = 1:n
        A[i,1] = one(T)
        for k = 1:p
            A[i,k+1] = A[i,k]*x[i]
        end
        A[i, p+2] = -y[i]*x[i]
        for k = 2:q
            A[i, p+1+k] = A[i, p+k] * x[i]
        end
    end

    A \ y

end



function rational_fit{T<:Number}(x::AbstractArray{T}, y::AbstractVector{T}, p, q)

    coefs0  = linear_rational_fit(x, y, p, q)

    

end

