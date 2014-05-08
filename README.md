# CurveFit

A package that implements a few curve fitting functions.

## Linear Least squares

Linear least square is commonly used technique to find approximation to a discrete
set of data. Given the sets of points `x[i]` and `y[i]` and a list of functions `f_i(x)`
the least squares method finds coefficients `a[i]` such that

```
a[1]*f_1(x) + a[2]*f_2(x) + ... + a[n]*f_n(x)
```
minimizes the squares of the errors in relation to `y[i]`.

The basic function is `least_squares`:
```
coefs = least_squares(A, y)
```
where `A[:,i] = f_i(x)`. While usually `x` is a single variable, in general if several
independent variables are used, the same procedure can be used: 
`A[:,i] = f_i(x1, x2, ..., xn)`.

Several typical cases are possible:
 * `linear_fit(x, y)` finds coeficients `a` and `b` for `y[i] = a + b*x[i]`
 * `power_fit(x, y)` finds coeficients `a` and `b` for `y[i] = a *x[i]^b`
 * `log_fit(x, y)` finds coeficients `a` and `b` for `y[i] = a + b*log(x[i])`
 * `exp_fit(x, y)` finds coeficients `a` and `b` for `y[i] = a*exp(b*x[i])`
 * `poly_fit(x, y, n)` finds coeficients `a[k]`  for 
   `y[i] = a[1] + a[2]*x[i] + a[3]*x[i]^2 + a[n+1]*x[i]^n`

### Example
```
using PyPlot
x = linspace(0, 2, 100)
y0 = 1 .+ x .+ x.*x .+ randn(100)/10
a = CurveFit.poly_fit(x, y0, 2)
y0b = a[1] .+ a[2] .* x .+ a[3] .* x.^2 
plot(x, y0, "o", x, y0b, "r-", linewidth=3)
```

## Nonlinear least squares

Sometimes the fitting function is not linear on the coefficients. In this case, givem
an approximation of the coefficients, the function is linearized around this 
approximation and linear least squares is used to calculate a better 
approximation of the coefficients. This iteration is repeated until convergence is 
reached. The fitting function has the following form:
```
f(x1, x2, x3, ..., xn, a1, a2, ... ap) = 0
```
where `xi` are the known data points and `ai` are the coefficients that 
should be calculated. 

The basic function that implements the nonlinear least squares has the following 
interface:
```
nonlinear_fit(x, fun::Function, dflst, a0, eps=1e-8, maxiter=200)
```
where 
 * `x` is a matrix where each column represents on data variable.
 * `fun` is a function that evaluates the fitting function. This function 
   has two arguments `x` and `a`, both arrays.
 * `dflst` is a function that calculates the derivatives of the fitting coefficients, 
  such that `dflst(k, x, a) = d / da_k f(x1, ..., xn, a1, ..., ak, ... ap)`.
 * `a0` Initial guess of the fitting coefficients.
 * `eps` convergence criteria.
 * `maxiter` maximum number of iterations to achieve convergence.

File `king.jl` that is an example using the function `nonlinear`.

Numerical derivatives are often enough and function `makeDerivFun` creates
a function that uses central differences to calculate the derivative.

## King's law

In hotwire anemometry, a simple expression for the calibration curve of the probe 
is known as King's law, expressed as:
```
E^2 = A + B*sqrt(U)
```
where `E` is voltage on the anemometer bridge, `U` is the flow velocity.
The coefficients A and B are obtained from a calibration. The function
`linear_king_fit` estimates coefficients `A` and `B`.

A better approximation for the calibration curve is known as modified
King's law:
```
E^2 = A + B*U^n
```
Now, this is a nonlinear curve fit. The linear fit (`linear_king_fit`) is usually
a very good first guess for the coefficients (where `n=0.5`). This curve fit is 
implemented in function `king_fit`.





## Generic interface

When different types of curve fits can be used, a common interface may be used. 
For each type of curve fitting, a corresponding `type` is defined. The generic function 
`curve_fit`  performs the curve fitting and the function `apply_fit` uses the object
returned by `curve_fit` to calculate the approximation.

The following types are defined in this package:
 * `LinearFit`
 * `LogFit`
 * `PowerFit`
 * `PolyFit`
 * `ExpFit`
 * `LinearKingFit`
 * `KingFit`


### Example
```
using PyPlot
U = linspace(1, 20, 20)
E = sqrt(2 .+ 1 .* U .^ 0.45) + randn(20)/60
e = linspace(minimum(E), maximum(E), 50)

f1 = CurveFit.curve_fit(KingFit, E, U)
U1 = CurveFit.apply_fit(f1, e)

f2 = CurveFit.curve_fit(PolyFit, E, U, 5)
U2 = CurveFit.apply_fit(f2, e)

plot(U, E, "o", U1, e, "r-", U2, e, "g-", linewidth=3)
```



[![Build Status](https://travis-ci.org/pjabardo/CurveFit.jl.svg)](https://travis-ci.org/pjabardo/CurveFit.jl)
