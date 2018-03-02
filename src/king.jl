

"""
Original Kings law (1910) represents the relation between voltage and velocity
in a hotwire anemometer. 
The law is given by:

`E^2 = A + B*U^0.5`

This function estimates `A` and `B`.
"""
linear_king_fit(E, U) = linear_fit(sqrt.(U), E .* E)

"Type that represents a Linear (original) King's law"
struct LinearKingFit <: LeastSquares
    coefs::Array{Float64,1}
    LinearKingFit(coefs) = new(copy(coefs))
end
LinearKingFit(E,U) = LinearKingFit(linear_king_fit(E, U))

(f::LinearKingFit)(E) = ( (E.*E .- f.coefs[1]) ./ f.coefs[2] ) .^ 2


"Equation that computes the error of the modified King's law "
kingfun(x, a) = a[1] + a[2] * x[2] ^ a[3] - x[1]*x[1]

"""
Uses nonlinear least squares to fit the modified King's law:

`E^2 = A + B*U^n`

The Original (linear) King's law is used to estimate `A` and `B` when `n=1/2`. 
This initial value is used as an initial guess for fitting the nonlinear modified King's law
using the function `nonlinear_fit`.
"""       
function king_fit(E, U, eps=1e-8, maxiter=200)

    f = linear_king_fit(E, U)
    a0 = [f[1], f[2], 0.5]

    coefs, converged, niter = nonlinear_fit(hcat(E,U), kingfun, a0, eps, maxiter)
    return coefs
end


"Type that represents the modified King's law "
struct KingFit <: LeastSquares
    coefs::Array{Float64,1}
    KingFit(coefs) = new(copy(coefs))
end
KingFit(E,U, eps=1e-8, maxiter=200) = KingFit(king_fit(E, U, eps, maxiter))

(f::KingFit)(E) = ( (E*E - f.coefs[1]) / f.coefs[2]) ^ (1 / f.coefs[3])





