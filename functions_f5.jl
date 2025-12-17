                

# Λ(κ, μ) = ( sin(√κ q1) / ( μ sin(√κ q2) ) )^(1/3)
function Λ(κ, μ)
    κ, μ = promote(κ, μ)
    r  = sqrt(κ)
    q1 = q₁(κ, μ)
    q2 = q₂(κ, μ)

    s1 = sin(r * q1)
    s2 = sin(r * q2)

    return cbrt_t(s1 / (μ * s2))
end

function f₅(x, κ, μ)
    x, κ, μ = promote(x, κ, μ)

    r  = sqrt(κ)
    q1 = q₁(κ, μ)
    q2 = q₂(κ, μ)
    G  = Γ(κ, μ)
    L  = Λ(κ, μ)

    s1 = sin(r * q1)
    s2 = sin(r * q2)

    sx = sin(x)
    cx = cos(x)

    # sqrt(1 - Λ^2 sin^2(x))
    inner = one(sx) - (L * L) * (sx * sx)
    rt    = sqrt(inner)

    bracket = s1 * cx + s2 * rt

    pref = (L * L * L) * s2 * (sx * sx * sx)    # Λ^3 sin(√κ q2) sin^3(x)
    tail = G * (sin(r) * sin(r))                # Γ sin^2(√κ)

    return pref * bracket - tail
end

@inline twice(x)  = x + x
@inline three(x)  = x + x + x
@inline inv3(x)   = inv(three(one(x)))
@inline cbrt_t(x) = x ^ inv3(x)

function f₅θ(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)

    r  = sqrt(κ)
    q1 = q₁(κ, μ)
    q2 = q₂(κ, μ)
    G  = Γ(κ, μ)
    L  = Λ(κ, μ)

    s1 = sin(r * q1)
    s2 = sin(r * q2)

    sθ = sin(θ)
    cθ = cos(θ)

    # rt = sqrt(1 - L^2 sin^2(θ))
    A  = one(sθ) - (L * L) * (sθ * sθ)
    rt = sqrt(A)

    bracket  = s1 * cθ + s2 * rt
    pref     = (L * L * L) * s2 * (sθ * sθ * sθ)   # Λ^3 sin(√κ q2) sin^3 θ
    # tail = G * sin(r)^2  (independent of θ) → derivative 0

    # d/dθ pref = Λ^3 sin(√κ q2) * 3 sin^2 θ cos θ
    prefθ = (L * L * L) * s2 * (three(one(sθ)) * (sθ * sθ * cθ))

    # d/dθ rt = -(L^2 sin θ cos θ)/rt
    rtθ = -((L * L) * (sθ * cθ)) / rt

    # d/dθ bracket = s1*(-sin θ) + s2*rtθ
    bracketθ = -(s1 * sθ) + s2 * rtθ

    return prefθ * bracket + pref * bracketθ
end

function Intx_5(κ, μ)
    κ, μ = promote(κ, μ)
    r = sqrt(κ)

    L = Λ(κ, μ)

    y = sin(r)
    c = cos(r)

    xL = L + c
    xR = L - c

    a = atan(y, xL)
    b = atan(y, xR)

    return (a, b)
end