function bound_angle(θ::Real, min::Real = 0)
    if min <= θ < (min + 2π)
        return θ
    else
        return θ - 2π*floor((θ-min)/2π)
    end
end

function isapproxangle(θ₁,θ₂; atol=1e-12, rtol=1e-12)
    diff = θ₁ - θ₂
    rem = (diff%2π)/(2π)
    return isapprox(rem,0,atol=atol) || isapprox(rem,1,atol=atol)
end
