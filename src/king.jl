


linear_king_fit(E, U) = linear_fit(sqrt(U), E .* E)
type LinearKingFit
    coefs::Array{Float64,1}
    LinearKingFit(coefs) = new(copy(coefs))
end
LinearKingFit(E,U) = LinearKingFit(linear_king_fit(E, U))

apply_fit(f::LinearKingFit, E) = ( (E.*E .- f.coefs[1]) ./ f.coefs[2] ) .^ 2



kingfun(x, a) = a[1] + a[2] * x[2] ^ a[3] - x[1]*x[1]
function deriv_kingfun(k, x, a)
    if k == 1
        return 1.0
    elseif k == 2
        return x[2]^a[3]
    else
        return a[2] * x[2]^a[3] * log(x[2])
    end
end


        
function king_fit(E, U, eps=1e-8, maxiter=200)

    f = linear_king_fit(E, U)
    a0 = [f[1], f[2], 0.5]
    #a0 = [1.0, 1.0, 0.4]

    nonlinear_fit(hcat(E,U), kingfun, deriv_kingfun, a0, eps, maxiter)

end

type KingFit <: LeastSquares
    coefs::Array{Float64,1}
    KingFit(coefs) = new(copy(coefs))
end
KingFit(E,U, eps=1e-8, maxiter=200) = KingFit(king_fit(E, U, eps, maxiter))

apply_fit(f::KingFit, E) = ( (E.*E .- f.coefs[1]) ./ f.coefs[2]) .^ (1./f.coefs[3])





