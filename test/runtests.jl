using CurveFit
using Base.Test  
using Polynomials
# write your own tests here
@test 1 == 1


# Testing linear fit
x = [linspace(1,10,10);]
fun0(x) = 1.0 + 2.0.*x
y = fun0(x)
f = linear_fit(x, y)

@test_approx_eq_eps f[1] 1.0 1e-7
@test_approx_eq_eps f[2] 2.0 1e-7

f = curve_fit(LinearFit, x, y)
@test_approx_eq_eps f(1.5) fun0(1.5) 1e-7

# power
fun1(x) = 1.0 + 2.0*log(x)
y = fun1(x)
f = log_fit(x, y)
@test_approx_eq_eps f[1] 1.0 1e-7
@test_approx_eq_eps f[2] 2.0 1e-7
f = curve_fit(LogFit, x, y)
@test_approx_eq_eps f(1.5) fun1(1.5) 1e-7


# log
fun2(x) = 2.0*x.^0.8
y = fun2(x)
f = power_fit(x, y)
@test_approx_eq_eps f[1] 2.0 1e-7
@test_approx_eq_eps f[2] 0.8 1e-7
f = curve_fit(PowerFit, x, y)
@test_approx_eq_eps f(1.5) fun2(1.5) 1e-7


# Exp
fun3(x) = 2.0*exp(0.8*x)
y = fun3(x)
f = exp_fit(x, y)
@test_approx_eq_eps f[1] 2.0 1e-7
@test_approx_eq_eps f[2] 0.8 1e-7
f = curve_fit(ExpFit, x, y)
@test_approx_eq_eps f(1.5) fun3(1.5) 1e-7


# Poly
fun4(x) = 1.0 + 2.0*x + 3.0*x.^2 + 0.5*x.^3
y = fun4(x)
f = poly_fit(x, y, 4)
@test_approx_eq_eps f[1] 1.0 1e-7
@test_approx_eq_eps f[2] 2.0 1e-7
@test_approx_eq_eps f[3] 3.0 1e-7
@test_approx_eq_eps f[4] 0.5 1e-7
@test_approx_eq_eps f[5] 0.0 1e-7

f = curve_fit(Poly, x, y, 4)
@test_approx_eq_eps f(1.5) fun4(1.5) 1e-7

# King's law
U = [linspace(1, 20, 20);]
A = 5.0
B = 1.5
n = 0.5
E = sqrt(A + B*U.^n)
fun5(E) = ((E.^2 - A)/B).^(1./n)

f = linear_king_fit(E, U)
@test_approx_eq_eps f[1] A 1e-7
@test_approx_eq_eps f[2] B 1e-7
f= curve_fit(LinearKingFit, E, U)
@test_approx_eq_eps f(3.0) fun5(3.0) 1e-7

# Modified King's law
n = 0.42

E = sqrt(A + B*U.^n)
fun6(E) = ((E.^2 - A)/B).^(1./n)

f = king_fit(E, U)
@test_approx_eq_eps f[1] A 1e-7
@test_approx_eq_eps f[2] B 1e-7
@test_approx_eq_eps f[3] n 1e-7
f= curve_fit(KingFit, E, U)
@test_approx_eq_eps f(3.0) fun6(3.0) 1e-5



# Linear Rational fit
r =  RationalPoly([1.0, 0.0, -2.0], [1.0, 2.0, 3.0])
y = r(x)
f = linear_rational_fit(x, y, 2, 3)
@test_approx_eq_eps f[1] 1.0 1e-8
@test_approx_eq_eps f[2] 0.0 1e-8
@test_approx_eq_eps f[3] -2.0 1e-8
@test_approx_eq_eps f[4] 2.0 1e-8
@test_approx_eq_eps f[5] 3.0 1e-8
@test_approx_eq_eps f[6] 0.0 1e-8

# Nonlinear Rational fit
f = rational_fit(x, y, 2, 3)

@test_approx_eq_eps f[1] 1.0 1e-7
@test_approx_eq_eps f[2] 0.0 1e-7
@test_approx_eq_eps f[3] -2.0 1e-7
@test_approx_eq_eps f[4] 2.0 1e-7
@test_approx_eq_eps f[5] 3.0 1e-7
@test_approx_eq_eps f[6] 0.0 1e-7



f = curve_fit(RationalPoly, x, y, 2, 3)
@test_approx_eq_eps f(1.5) r(1.5) 1e-8
@test_approx_eq_eps f(4.5) r(4.5) 1e-8





