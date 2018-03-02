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


function makeDerivFun(fun, na, da=1e-7)
    

    
    if length(da)==1
        da = [da for i in 1:na]
    end

    dfun(k, x, a) = begin
        a[k] = a[k] + da[k]
        x1 = fun(x, a)
        a[k] = a[k] - 2*da[k]
        x2 = fun(x, a)
        a[k] = a[k] + da[k]
        (x1 - x2) / (2*da[k])
    end

    return dfun
end




        
        
    


        
