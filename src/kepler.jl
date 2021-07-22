function kepler(M::Real, e::Real; tol = 1e-12, maxit = 1000)
    if e == 1           # Parabolic case
        throw(ArgumentError("Parabolic case (e=1) has not been implemented"))
    elseif e < 1        # Elliptical case
        if e < 0.08
            E = M
        else
            E = π
        end
        Ep = E + 2*tol
        it = 1
        while (abs(E-Ep) > tol) && (it < maxit)
            it = it + 1
            Ep = E
            E = E - (E - e*sin(E) - M)/(1-e*cos(E))
        end
        return E
    else                # Hyperbolic case
        if abs(M) > 4π 
            H = sign(M) * 4π
        else
            H = M
        end
        Hp = H + 2*tol
        it = 1
        while (abs(H-Hp) > tol) && (it < maxit)
            it = it + 1
            Hp = H
            H = H + (M - e*sinh(H) + H)/(e*cosh(H) - 1)
            if it == maxit
                @show it, M, e
            end
        end
        return H
    end
end
