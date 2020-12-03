# Nonlinear curve fitting.

"""
# Nonlinear multi-variable least squares

Uses a Newton like algorithm to compute the least squares fit 
of a model

## Parameters:

 * `x` an array where each colum is a variable
 * `fun` The function that should be fitted (more later on)
 * `a0` An initial guess for each fitting parameter
 * `eps` Convergence criterion
 * `maxiter` Maximum number of iterations

## Specifying the model

The model is specified in argument `fun` that should be callable according
to the following signature:

`r = fun(x, a)`

where both `x` and `a` are vectors

The way the algorithm is implemented allows for the function to be implicit. This means
that we are trying to fit the relationship `fun(x, a) = 0`. 

## Initial guess

The initial guess is important and should be as close as possible to the expeted values. It should, at least, give an order of magnitude of each parameter. Zero values are not recommended since in this 
case no order of magnitude is available and it is assumed. The initial step is assumed to be `a0/10`.

## Convergence criterion

The initial guess `a0` will probably not fit the data very well. The largest value
of `maxerr = maxabs(fun(x, a0))` (maximum error) is used as a reference and 
every time a new approximation is computed, if the change in paramaters is compared to this
reference error. If it is small, `maxabs(da) < eps*maxerr`, the algorithm converged.

## Return values

A tuple containing:

 * The estimated coefficients
 * A `Bool` stating whether the algorithm converged
 * Number of iterations it took to converge.

"""
function nonlinear_fit(x, fun, a0, eps=1e-8, maxiter=200)

    na = length(a0)
    np = size(x, 1)
    nv = size(x, 2)
    for i in 1:na
        if a0[i] == 0
            a0[i] = 0.01
        end
    end

    da = a0 / 10
    a1 = zeros(na)
    a = zeros(na)
    for k = 1:na
        a1[k] = a0[k] - da[k]
    end
    xp = zeros(nv)
    A = zeros(np, na)
    r0 = zeros(np)
    r1 = zeros(np)
    for p = 1:np
        for k = 1:nv
            xp[k] = x[p,k]
        end
        r0[p] = fun(xp, a0)
        r1[p] = fun(xp, a1)
    end
    maxerr = maximum(abs, [maximum(abs, r0), maximum(abs, r1)])
    iter = 1
    convergence = false
    for i = 1:maxiter
        iter = i
        for p = 1:np
            for k = 1:nv
                xp[k] = x[p,k]
            end
            for k = 1:na
                aa = a1[k]
                a1[k] = a1[k] + da[k]
                r = fun(xp, a1)
                a1[k] = aa
                  
                A[p,k] = -(r1[p] - r) / da[k]
            end
        end
        da = A \ r1
        for k = 1:na
            a1[k] = a1[k] - da[k]
        end
        for p = 1:np
            r0[p] = r1[p]
            for k = 1:nv
                xp[k] = x[p,k]
            end
            r1[p] = fun(xp, a1)
        end
        if maximum(abs, r1) < eps * maxerr
            convergence = true
            break
        end
    end

    return a1, convergence, iter
end


"""
   a = gauss_newton_fit(x, y, fun, ∇fun!, a0[[, eps,] maxiter])

Gauss-Newton nonlinear least squares. Given vectors `x` and `y`, the tries to fit parameters `a` to 
a function `f` using least squares approximation:

 ``y = f(x, a₁, ..., aₙ)``

For more general approximations, see [`gauss_newton_fit`](@ref).

### Arguments:

 * `x` Vector with x values
 * `y` Vector with y values
 * `fun` a function that is called as `fun(x, a)` where `a` is a vector of parameters.
 * `∇fun!` A function that calculares the derivatives with respect to parameters `a`
 * `a0` Vector with the initial guesses of the parameters
 * `eps` Maximum residual for convergence
 * `maxiter` Maximum number of iterations for convergence

## Return value

A vector with the convrged array. If no convergence is achieved, the function throws an error.

## Specification of the fitting function

The function that should be fitted shoud be specified by Julia funcion with the following signature:

```julia
fun(x::T, a::AbstractVector{T}) where {T<:Number}
```

The derivatives with respect to each fitting parameter `a[i]` should have the following signature:

```julia
∇fun!(x::T, a::AbstractVector{T}, df::AbstractVector{T}) where {T<:Number}
```

No return value is expected and the derivatives are returned in argument `df`.

## Initial approximation (guess)

If the initial approximation is not good enough, divergence is possible. 

**Careful** with parameters close to 0. The initial guess should never be 0.0 because the initial
value of the parameter is used as reference value for computing resiuduals.

## Convergence criteria

The argumento `maxiter` specifies the maximum number of iterations that should be carried out. 
At each iteration, 

``aₖⁿ⁺¹ = aₖⁿ + δₖ``

Convergence is achieved when

``|δᵏ / aₖ⁰| < ε``

## Example
```julia
x = 1.0:10.0
a = [3.0, 2.0, 1.0]
y = a[1] + a[2]*x + a[3]*x^2
fun(x, a) = a[1] + a[2]*x + a[3]*x^2

function ∇fun!(x, a, df) 
    df[1] = 1.0
    df[2] = x
    df[3] = x^2
end

a = gauss_newton_fit(x, y, fun, ∇fun!, [0.5, 0.5, 0.5], 1e-8, 30)
```
"""
function gauss_newton_fit(x::AbstractVector{T}, y::AbstractVector{T}, fun, ∇fun!, a0::AbstractVector{T},
                          eps=1e-8, maxiter=200) where {T<:Number}
    P = length(x) # Number of points
    N = length(a0) # Number of parameters
    
    xi = zero(T)
    df = zeros(T, N)
    a = zeros(T, N)
    for i in 1:N
        a[i] = a0[i]
        if a[i] == 0
            a[i] = 0.01
        end
    end

    A = zeros(T, N, N)
    b = zeros(T, N)

    δref = abs.(a)
    maxerr = zero(T)
    for iter = 1:maxiter
        A .= zero(T)
        b .= zero(T)
        for i in 1:P
            xi = x[i]
            yi = y[i]
            f = fun(xi, a) - yi
            ∇fun!(xi, a, df)

            # Assemble LHS
            for k in 1:N
                for j in 1:N
                    A[j,k] += df[j] * df[k]
                end
            end
            # Assemble RHS
            for j in 1:N
                b[j] -= f * df[j]
            end
        end

        δ = A\b
        a .+= δ

        # Verify convergence:
        maxerr = maximum(abs, δ./δref)
        if maxerr < eps
            return(a)
        end
        
    end

    error("gauss_newton_fit failed to converge in $maxiter iterations with relative residual of $maxerr !")
    
    return(a)
end

"""
   a = gauss_newton_generic_fit(x, y, fun, ∇fun!, a0[[, eps,] maxiter])

Gauss-Newton nonlinear least squares. Given matrix `x` the function tries to fit parameters `a` to 
an implicit function `f` using least squares approximation:

``f(x₁,..., xₘ, a₁, ..., aₙ) = 0``

This function is very similar to [`gauss_newton_fit`](@ref) but more generic. It doesn't limit the number of 
independent variables and doesn't require a `y` LHS.


### Arguments:

 * `x` Matrix where each column is a variable
 * `
 * `fun` a function that is called as `fun(x, a)` where `a` is a vector of parameters.
 * `∇fun!` A function that calculares the derivatives with respect to parameters `a`
 * `a0` Vector with the initial guesses of the parameters
 * `eps` Maximum residual for convergence
 * `maxiter` Maximum number of iterations for convergence

## Return value

A vector with the convrged array. If no convergence is achieved, the function throws an error.

## Specification of the fitting function

The function that should be fitted shoud be specified by Julia funcion with the following signature:

```julia
fun(x::T, a::AbstractVector{T}) where {T<:Number}
```

The derivatives with respect to each fitting parameter `a[i]` should have the following signature:

```julia
∇fun!(x::T, a::AbstractVector{T}, df::AbstractVector{T}) where {T<:Number}
```

No return value is expected and the derivatives are returned in argument `df`.

## Initial approximation (guess)

If the initial approximation is not good enough, divergence is possible. 

**Careful** with parameters close to 0. The initial guess should never be 0.0 because the initial
value of the parameter is used as reference value for computing resiuduals.

## Convergence criteria

The argumento `maxiter` specifies the maximum number of iterations that should be carried out. 
At each iteration, 

``aₖⁿ⁺¹ = aₖⁿ + δₖ``

Convergence is achieved when

``|δᵏ / aₖ⁰| < ε``

## Example
```julia
x = 1.0:10.0
a = [3.0, 2.0, 1.0]
y = a[1] + a[2]*x + a[3]*x^2
fun(x, a) = a[1] + a[2]*x + a[3]*x^2

function ∇fun!(x, a, df) 
    df[1] = 1.0
    df[2] = x
    df[3] = x^2
end

a = gauss_newton_fit(x, y, fun, ∇fun!, [0.5, 0.5, 0.5], 1e-8, 30)
```

"""
function gauss_newton_generic_fit(x::AbstractMatrix{T}, fun, ∇fun!, a0::AbstractVector{T},
                          eps=1e-8, maxiter=200) where {T<:Number}
    P = size(x, 1) # Number of points
    M = size(x, 2) # Number of columns (variables)
    N = length(a0) # Number of parameters
    
    xi = zeros(T, M)
    df = zeros(T, N)
    a = zeros(T, N)
    for i in 1:N
        a[i] = a0[i]
        if a[i] == 0
            a[i] = 0.01
        end
    end

    A = zeros(T, N, N)
    b = zeros(T, N)

    δref = abs.(a)
    maxerr = zero(T)
    for iter = 1:maxiter
        A .= zero(T)
        b .= zero(T)
        for i in 1:P
            for k in 1:M
                xi[k] = x[i,k]
            end
            f = fun(xi, a)
            
            ∇fun!(xi, a, df)
            # Assemble LHS
            for k in 1:N
                for j in 1:N
                    A[j,k] += df[j] * df[k]
                end
            end
            # Assemble RHS
            for j in 1:N
                b[j] -= f * df[j]
            end
        end

        δ = A\b
        a .+= δ

        # Verify convergence:
        maxerr = maximum(abs, δ./δref)
        if maxerr < eps
            return(a)
        end
        
    end

    error("gauss_newton_fit failed to converge in $maxiter iterations with relative residual of $maxerr !")
    
    return(a)
end



        
        
    
"""
   a = secant_nls_fit(x, y, fun, ∇fun!, a0[[, eps,] maxiter])

Secant/Gauss-Newton nonlinear least squares. DOESN'T NEED A DERIVATIVE FUNCTION. Given vectors `x` and `y`, the tries to fit parameters `a` to 
a function `f` using least squares approximation:

 ``y = f(x, a₁, ..., aₙ)``

For more general approximations, see [`gauss_newton_fit`](@ref).

### Arguments:

 * `x` Vector with x values
 * `y` Vector with y values
 * `fun` a function that is called as `fun(x, a)` where `a` is a vector of parameters.
 * `∇fun!` A function that calculares the derivatives with respect to parameters `a`
 * `a0` Vector with the initial guesses of the parameters
 * `eps` Maximum residual for convergence
 * `maxiter` Maximum number of iterations for convergence

## Return value

A vector with the convrged array. If no convergence is achieved, the function throws an error.

## Specification of the fitting function

The function that should be fitted shoud be specified by Julia funcion with the following signature:

```julia
fun(x::T, a::AbstractVector{T}) where {T<:Number}
```

The derivatives with respect to each fitting parameter `a[i]` should have the following signature:

```julia
∇fun!(x::T, a::AbstractVector{T}, df::AbstractVector{T}) where {T<:Number}
```

No return value is expected and the derivatives are returned in argument `df`.

## Initial approximation (guess)

If the initial approximation is not good enough, divergence is possible. 

**Careful** with parameters close to 0. The initial guess should never be 0.0 because the initial
value of the parameter is used as reference value for computing resiuduals.

## Convergence criteria

The argumento `maxiter` specifies the maximum number of iterations that should be carried out. 
At each iteration, 

``aₖⁿ⁺¹ = aₖⁿ + δₖ``

Convergence is achieved when

``|δᵏ / aₖ⁰| < ε``

## Example
```julia
x = 1.0:10.0
a = [3.0, 2.0, 1.0]
y = a[1] + a[2]*x + a[3]*x^2
fun(x, a) = a[1] + a[2]*x + a[3]*x^2

a = secant_nls_fit(x, y, fun, ∇fun!, [0.5, 0.5, 0.5], 1e-8, 30)
```
"""

function secant_nls_fit(x::AbstractVector{T}, y::AbstractVector{T}, fun, aguess::AbstractVector{T},
                          eps=1e-8, maxiter=200) where {T<:Number}
    P = length(x) # Number of points
    N = length(aguess) # Number of parameters
    
    xi = zero(T)
    df = zeros(T, N)
    a = zeros(T, N)
    for i in 1:N
        a[i] = aguess[i]
        if a[i] == 0
            a[i] = 0.01
        end
    end

    δ = a .* (one(T)/20)
    f1 = zeros(T, P)
    a .+= δ
    A = zeros(T, N, N)
    b = zeros(T, N)

    δref = abs.(a)
    maxerr = zero(T)
    for iter = 1:maxiter

        A .= zero(T)
        b .= zero(T)
        f1 .= fun.(x, Ref(a))
        for i in 1:P
            xi = x[i]
            yi = y[i]
            f = f1[i] - yi
            for k in 1:N
                a[k] -= δ[k]
                df[k] = (f1[i] - fun(xi, a)) / δ[k]
                a[k] += δ[k]
            end
            # Assemble LHS
            for k in 1:N
                for j in 1:N
                    A[j,k] += df[j] * df[k]
                end
            end
            # Assemble RHS
            for j in 1:N
                b[j] -= f * df[j]
            end
        end
        δ = A\b
        a .+= δ
        # Verify convergence:
        maxerr = maximum(abs, δ./δref)
        if maxerr < eps
            return(a)
        end
        
    end

    error("gauss_newton_fit failed to converge in $maxiter iterations with relative residual of $maxerr !")
    
    return(a)
end


        
