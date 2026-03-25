### Functions and derivatives for f₅
@inline function sinc_t(x)
    if x isa Interval
        return (inf(x) <= 0.0 <= sup(x)) ? interval(-1.0, 1.0) : sin(x) / x
    else
        return Base.sinc(x / pi)
    end
end


# function Λ(κ, μ)
#     κ, μ = promote(κ, μ)
#     r = sqrt(κ)
#     q1 = q₁(κ, μ)
#     q2 = q₂(κ, μ)

#     s2 = sin(r * q2)

#     num = r * (q1 / μ) * sinc_t(r * q1)
#     return cbrt_t(num / s2)
# end

function f₅(x, κ, μ)
    x, κ, μ = promote(x, κ, μ)

    r = sqrt(κ)
    q1 = q₁(κ, μ)
    q2 = q₂(κ, μ)
    G = Γ(κ, μ)
    L = Λ(κ, μ)

    s1 = sin(r * q1)
    s2 = sin(r * q2)

    sx = sin(x)
    cx = cos(x)

    inner = one(sx) - (L * L) * (sx * sx)
    rt = sqrt(inner)

    bracket = s1 * cx + s2 * rt

    pref = (L * L * L) * s2 * (sx * sx * sx)
    tail = G * (sin(r) * sin(r))

    return pref * bracket - tail
end

function f₅θ(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)

    r = sqrt(κ)
    q1 = q₁(κ, μ)
    q2 = q₂(κ, μ)
    G = Γ(κ, μ)
    L = Λ(κ, μ)

    s1 = sin(r * q1)
    s2 = sin(r * q2)

    sθ = sin(θ)
    cθ = cos(θ)

    A = one(sθ) - (L * L) * (sθ * sθ)
    rt = sqrt(A)

    bracket = s1 * cθ + s2 * rt
    pref = (L * L * L) * s2 * (sθ * sθ * sθ)

    prefθ = (L * L * L) * s2 * (three(one(sθ)) * (sθ * sθ * cθ))

    rtθ = -((L * L) * (sθ * cθ)) / rt

    bracketθ = -(s1 * sθ) + s2 * rtθ

    return prefθ * bracket + pref * bracketθ
end

function Intθ_5(κ, μ)
    κ, μ = promote(κ, μ)
    r = sqrt(κ)

    L = Λ(κ, μ)

    y = sin(r)
    c = cos(r)

    xL = L + c
    xR = L - c

    a = atan(y, xL)
    b = pi - atan(y, xR)

    return (a, b)
end