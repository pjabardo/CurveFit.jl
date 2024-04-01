using CurveFit
using Test
using LinearAlgebra
using Polynomials

lspf(start, stop, len) = range(start, stop = stop, length = len)

@testset "CurveFit.jl" begin

    # Testing linear fit
    x = lspf(1, 10, 10)
    fun0(x) = 1.0 + 2.0 * x
    y = fun0.(x)
    f = linear_fit(x, y)

    @test f[1]≈1.0 atol=1.0e-7
    @test f[2]≈2.0 atol=1.0e-7

    f = curve_fit(LinearFit, x, y)
    @test f(1.5)≈fun0(1.5) atol=1.0e-7

    # power
    fun1(x) = 1.0 + 2.0 * log(x)
    y = fun1.(x)
    f = log_fit(x, y)
    @test f[1]≈1.0 atol=1.0e-7
    @test f[2]≈2.0 atol=1.0e-7
    f = curve_fit(LogFit, x, y)
    @test f(1.5)≈fun1(1.5) atol=1.0e-7

    # log
    fun2(x) = 2.0 * x .^ 0.8
    y = fun2(x)
    f = power_fit(x, y)
    @test f[1]≈2.0 atol=1.0e-7
    @test f[2]≈0.8 atol=1.0e-7
    f = curve_fit(PowerFit, x, y)
    @test f(1.5)≈fun2(1.5) atol=1.0e-7

    # Exp
    fun3(x) = 2.0 * exp.(0.8 * x)
    y = fun3(x)
    f = exp_fit(x, y)
    @test f[1]≈2.0 atol=1.0e-7
    @test f[2]≈0.8 atol=1.0e-7
    f = curve_fit(ExpFit, x, y)
    @test f(1.5)≈fun3(1.5) atol=1.0e-7

    # Poly
    fun4(x) = 1.0 + 2.0 * x + 3.0 * x^2 + 0.5 * x^3
    y = fun4.(x)
    f = poly_fit(x, y, 4)
    @test f[1]≈1.0 atol=1.0e-7
    @test f[2]≈2.0 atol=1.0e-7
    @test f[3]≈3.0 atol=1.0e-7
    @test f[4]≈0.5 atol=1.0e-7
    @test f[5]≈0.0 atol=1.0e-7

    f = curve_fit(Polynomial, x, y, 4)
    @test f(1.5)≈fun4(1.5) atol=1.0e-7

    # Polynomials with large numbers
    coefs = [80.0, -5e-18, -7e-20, -1e-36]
    P = Polynomial(coefs)
    x1 = 1e10 * (0:0.1:5)
    y1 = P.(x1)
    P2 = curve_fit(Polynomial, x1, y1, 3)
    c = coeffs(P2)
    @test coefs[1]≈c[1] rtol=1e-5
    @test coefs[2]≈c[2] rtol=1e-5
    @test coefs[3]≈c[3] rtol=1e-5
    @test coefs[4]≈c[4] rtol=1e-5

    # King's law
    U = [lspf(1, 20, 20);]
    A = 5.0
    B = 1.5
    n = 0.5
    E = sqrt.(A .+ B * U .^ n)
    fun5(E) = ((E .^ 2 - A) / B) .^ (1 ./ n)

    f = linear_king_fit(E, U)
    @test f[1]≈A atol=1.0e-7
    @test f[2]≈B atol=1.0e-7
    f = curve_fit(LinearKingFit, E, U)
    @test f(3.0)≈fun5(3.0) atol=1.0e-7

    # Modified King's law
    n = 0.42

    E = sqrt.(A .+ B .* U .^ n)
    fun6(E) = ((E^2 - A) / B)^(1 / n)

    f = king_fit(E, U)
    @test f[1]≈A atol=1.0e-7
    @test f[2]≈B atol=1.0e-7
    @test f[3]≈n atol=1.0e-7
    f = curve_fit(KingFit, E, U)
    @test f(3.0)≈fun6(3.0) atol=1.0e-5

    # Linear Rational fit
    r = RationalPoly([1.0, 0.0, -2.0], [1.0, 2.0, 3.0])
    y = r.(x)
    f = linear_rational_fit(x, y, 2, 3)
    @test f[1]≈1.0 atol=1.0e-8
    @test f[2]≈0.0 atol=1.0e-8
    @test f[3]≈-2.0 atol=1.0e-8
    @test f[4]≈2.0 atol=1.0e-8
    @test f[5]≈3.0 atol=1.0e-8
    @test f[6]≈0.0 atol=1.0e-8

    # Nonlinear Rational fit
    f = rational_fit(x, y, 2, 3)

    @test f[1]≈1.0 atol=1.0e-7
    @test f[2]≈0.0 atol=1.0e-7
    @test f[3]≈-2.0 atol=1.0e-7
    @test f[4]≈2.0 atol=1.0e-7
    @test f[5]≈3.0 atol=1.0e-7
    @test f[6]≈0.0 atol=1.0e-7

    f = curve_fit(RationalPoly, x, y, 2, 3)
    @test f(1.5)≈r(1.5) atol=1.0e-8
    @test f(4.5)≈r(4.5) atol=1.0e-8

    # Gauss-Newton curve fitting. Linear problem:
    x = 1.0:10.0
    a0 = [3.0, 2.0, 1.0]
    fun7(x, a) = a[1] + a[2] * x + a[3] * x^2
    y = fun7.(x, Ref(a0))

    function ∇fun7!(x, a, df)
        df[1] = 1.0
        df[2] = x
        df[3] = x^2
    end

    a = gauss_newton_fit(x, y, fun7, ∇fun7!, [0.5, 0.5, 0.5], 1e-8, 30)

    @test a[1]≈a0[1] rtol=1e-7
    @test a[2]≈a0[2] rtol=1e-7
    @test a[3]≈a0[3] rtol=1e-7

    # Gauss-Newton curve fitting. Nonlinear problem:
    x = 1.0:10.0
    fun8(x, a) = a[1] + a[2] * x^a[3]
    a0 = [3.0, 2.0, 0.7]
    y = fun8.(x, Ref(a0))
    function ∇fun8!(x, a, df)
        df[1] = 1.0
        df[2] = x^a[3]
        df[3] = a[2] * x^a[3] * log(x)
    end
    a = gauss_newton_fit(x, y, fun8, ∇fun8!, [0.5, 0.5, 0.5], 1e-8, 30)

    @test a[1]≈a0[1] rtol=1e-7
    @test a[2]≈a0[2] rtol=1e-7
    @test a[3]≈a0[3] rtol=1e-7

    # Gauss-Newton curve fitting (generic interface). Linear problem:
    x = 1.0:10.0
    a0 = [3.0, 2.0, 1.0]
    fun9(x, a) = a[1] + a[2] * x[1] + a[3] * x[1]^2 - x[2]
    P = length(x)
    X = zeros(P, 2)
    for i in 1:P
        X[i, 1] = x[i]
        X[i, 2] = a0[1] + a0[2] * x[i] + a0[3] * x[i]^2
    end

    function ∇fun9!(x, a, df)
        df[1] = 1.0
        df[2] = x[1]
        df[3] = x[1]^2
    end

    a = gauss_newton_generic_fit(X, fun9, ∇fun9!, [0.5, 0.5, 0.5], 1e-8, 30)

    @test a[1]≈a0[1] rtol=1e-7
    @test a[2]≈a0[2] rtol=1e-7
    @test a[3]≈a0[3] rtol=1e-7

    # Gauss-Newton curve fitting (generic interface). Nonlinear problem:
    x = 1.0:10.0
    funA(x, a) = a[1] + a[2] * x[1]^a[3] - x[2]
    a0 = [3.0, 2.0, 0.7]
    P = length(x)
    X = zeros(P, 2)
    for i in 1:P
        X[i, 1] = x[i]
        X[i, 2] = a0[1] + a0[2] * x[i]^a0[3]
    end
    function ∇funA!(x, a, df)
        df[1] = 1.0
        df[2] = x[1]^a[3]
        df[3] = a[2] * x[1]^a[3] * log(x[1])
    end
    a = gauss_newton_generic_fit(X, funA, ∇funA!, [0.5, 0.5, 0.5], 1e-8, 30)

    @test a[1]≈a0[1] rtol=1e-7
    @test a[2]≈a0[2] rtol=1e-7
    @test a[3]≈a0[3] rtol=1e-7

    # Gauss-Newton curve fitting (generic interface). Nonlinear problem:
    U = 0.5:0.5:10
    a0 = [2.0, 1.0, 0.35]
    E = @. sqrt(a0[1] + a0[2] * U^a0[3])

    X = hcat(E, U)
    funB(x, a) = a[1] + a[2] * x[2]^a[3] - x[1]^2
    function ∇funB!(x, a, df)
        df[1] = 1.0
        df[2] = x[2]^a[3]
        df[3] = a[2] * x[2]^a[3] * log(x[2])
    end

    a = gauss_newton_generic_fit(X, funB, ∇funB!, [0.5, 0.5, 0.5], 1e-8, 30)

    @test a[1]≈a0[1] rtol=1e-7
    @test a[2]≈a0[2] rtol=1e-7
    @test a[3]≈a0[3] rtol=1e-7

    # Secant method NLS curve fitting. Linear problem:
    x = 1.0:10.0
    a0 = [3.0, 2.0, 1.0]
    funC(x, a) = a[1] + a[2] * x + a[3] * x^2
    y = funC.(x, Ref(a0))

    a = CurveFit.secant_nls_fit(x, y, funC, [0.5, 0.5, 0.5], 1e-8, 30)

    @test a[1]≈a0[1] rtol=1e-7
    @test a[2]≈a0[2] rtol=1e-7
    @test a[3]≈a0[3] rtol=1e-7

    # Gauss-Newton curve fitting. Nonlinear problem:
    x = 1.0:10.0
    funD(x, a) = a[1] + a[2] * x^a[3]
    a0 = [3.0, 2.0, 0.7]
    y = funD.(x, Ref(a0))
    a = CurveFit.secant_nls_fit(x, y, funD, [0.5, 0.5, 0.5], 1e-8, 30)

    @test a[1]≈a0[1] rtol=1e-7
    @test a[2]≈a0[2] rtol=1e-7
    @test a[3]≈a0[3] rtol=1e-7

    @testset "ExpSumFit" begin
        include("expsumfit.jl")
    end
end
