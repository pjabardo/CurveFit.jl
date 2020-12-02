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

function nonlinear_fit0(x, fun, dflst, a0, eps=1e-8, maxiter=200)

    na = length(a0)
    np = size(x, 1)
    for i in 1:na
        if a0[i] == 0
            a0[i] = 0.01
        end
    end
    aref = eps .* abs(a0)
    
    A = zeros(Float64, np, na)
    r = zeros(Float64, np)
    iter = 1
    convergence = false
    for i = 1:maxiter
        iter = i
        for p = 1:np
            xp = x[p,:]
            r[p] = -fun(xp, a0)
            for k = 1:na
                A[p,k] = dflst(k, xp, a0)
            end
        end
        da = least_squares(A, r)
        a0 = a0 + da
        if maxabs(da/aref) < eps
            convergence = true
            break
        end
    end

    return a0, convergence, iter
end


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

myfun0(x, a) = a[1] + a[2]*x[1]^a[3] - x[2]^2
function mydfun0!(x, a, df)
    df[1] = 1.0
    df[2] = x[1]^a[3]
    df[3] = a[2] * x[1]^a[3] * log(x[1])
end

myfun1(x, a) = a[1] + a[2]*x[1] + a[3]*x[1]^2 - x[2]
function mydfun1!(x, a, df)
    df[1] = 1.0
    df[2] = x[1]
    df[3] = x[1]^2
end


        
        
    


        
