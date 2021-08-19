"""
    E = kepler(M, e; atol=1e-12, maxit=1000)

Solves [Kepler's equation](https://en.wikipedia.org/wiki/Kepler%27s_equation) for eccentric anomaly or
hyperbolic anomaly using Newton's method.
"""
function kepler(M::Real, e::Real; atol = 1e-12, maxit = 1000)
    if e == 1           # Parabolic case
        throw(ArgumentError("Parabolic case (e=1) has not been implemented"))
    elseif e < 1        # Elliptical case
        # # from Seppo Mikkola (slower, but better initial guess)
        # α = (1-e)/(4*e+0.5)
        # β = 0.5*M/(4*e+0.5)
        # z = copysign((abs(β) + √(β^2 + α^3))^(1/3), β)
        # s = z - α/z
        # s = s - 0.078*s^5/(1+e)
        # E = M + e*(3*s - 4*s^3)
        E = M
        # iterate until convergence
        err = 1+atol
        it = 1
        while (err > atol) && (it < maxit)
            it = it + 1
            Ep = E
            errp = err
            E = E - (E - e*sin(E) - M)/(1-e*cos(E))       # Newton's method
            # f = E - e*sin(E) - M
            # df = 1 - e*cos(E)
            # ddf = e*sin(E)
            # E = E - 2*f*df/(2*df^2 - f*ddf)               # Halley's method (slower, no accuracy improvement)
            err = abs(E-Ep)
            # if error gets worse, use bisection (makes up for poor initial guess)
            if err > errp
                E = (E+Ep)/2
            end
        end
        if it == maxit
            @warn "Kepler solver failed to converge: err = $err"
        end
        return E
    else                # Hyperbolic case
        if abs(M) > 4π 
            H = sign(M) * 4π
        else
            H = M
        end
        err = 1+atol
        it = 1
        while (err > atol) && (it < maxit)
            it = it + 1
            Hp = H
            errp = err
            H = H + (M - e*sinh(H) + H)/(e*cosh(H) - 1)     # Newton's method
            err = abs(H-Hp)
            # if error gets worse, use bisection (makes up for poor initial guess)
            if err > errp
                H = (H+Hp)/2
            end
        end
        if it == maxit
            @warn "Kepler solver failed to converge: err = $err"
        end
        return H
    end
end
