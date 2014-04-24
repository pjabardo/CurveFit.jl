# Nonlinear curve fitting.



function nonlinear_fit(x, fun::Function, dflst, a0, eps=1e-8, maxiter=200)

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
        if max(abs(da/aref)) < eps
            convergence = true
            break
        end
    end
    if !convergence
        println("NAO CONVERGIU")
    end
    return a0
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




        
        
    


        