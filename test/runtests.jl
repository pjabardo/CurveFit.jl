using CurveFit
using Base.Test  

# write your own tests here
@test 1 == 1


# Testing linear fit
x = [linspace(1,10,10);]
fun(x) = 1.0 + 2.0.*x
y = fun(x)
f = linear_fit(x, y1)

@test_approx_eq_eps f[1] 1.0 1e-7
@test_approx_eq_eps f[2] 2.0 1e-7

f = curve_fit(LinearFit, x, y)
@test_approx_eq_eps f(1.5) fun(1.5) 1e-7

# power
fun(x) = 1.0 + 2.0*log(x)
y = fun(x)
f = log_fit(x, y)
@test_approx_eq_eps f[1] 1.0 1e-7
@test_approx_eq_eps f[2] 2.0 1e-7
f = curve_fit(LogFit, x, y)
@test_approx_eq_eps f(1.5) fun(1.5) 1e-7


# log
fun(x) = 2.0*x.^0.8
y = fun(x)
f = power_fit(x, y)
@test_approx_eq_eps f[1] 2.0 1e-7
@test_approx_eq_eps f[2] 0.8 1e-7
f = curve_fit(PowerFit, x, y)
@test_approx_eq_eps f(1.5) fun(1.5) 1e-7


# Exp
fun(x) = 2.0*exp(0.8*x)
y = fun(x)
f = exp_fit(x, y)
@test_approx_eq_eps f[1] 2.0 1e-7
@test_approx_eq_eps f[2] 0.8 1e-7
f = curve_fit(ExpFit, x, y)
@test_approx_eq_eps f(1.5) fun(1.5) 1e-7


# Poly
fun(x) = 1.0 + 2.0*x + 3.0*x.^2 + 0.5*x.^3
y = fun(x)
f = poly_fit(x, y, 4)
@test_approx_eq_eps f[0] 1.0 1e-7
@test_approx_eq_eps f[1] 2.0 1e-7
f = curve_fit(Poly, x, y, 4)
@test_approx_eq_eps f(1.5) fun(1.5) 1e-7

