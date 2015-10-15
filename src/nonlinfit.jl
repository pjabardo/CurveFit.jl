# Nonlinear curve fitting.

using Debug

function nonlinear_fit(x, fun, a0, eps=1e-8, maxiter=200)

    na = length(a0)
    np = size(x, 1)
    nv = size(x, 2)
    for i in 1:na
        if a0[i] == 0
            a0[i] = 0.01
        end
    end
    aref = eps .* abs(a0)

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
    maxerr = max(maxabs(r0), maxabs(r1))
    iter = 1
    convergence = false
    for iter = 1:maxiter
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
        if maxabs(r1) < eps * maxerr
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
    for iter = 1:maxiter
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




        
        
    


        
