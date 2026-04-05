### Functions and derivatives for f₁ with κ < 0
function f₁neg(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    r = sqrt(-κ)
    a = θ + r * q₁neg(κ, μ)
    b = θ - r * q₂neg(κ, μ)
    term = inv(sqr(sinh(a))) - μ * inv(sqr(sinh(b)))
    return sinh(θ) * cosh(θ) - Γneg(κ, μ) * term
end

function f₁negθ(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    r = sqrt(-κ)
    a = θ + r * q₁neg(κ, μ)
    b = θ - r * q₂neg(κ, μ)
    Γv = Γneg(κ, μ)

    d_main = cosh(twice(θ))

    sa = sinh(a)
    ca = cosh(a)
    sb = sinh(b)
    cb = cosh(b)

    inv_sa3 = inv(sa * sa * sa)
    inv_sb3 = inv(sb * sb * sb)
    d_term = (-twice(ca) * inv_sa3) + (twice(μ) * cb * inv_sb3)

    return d_main - Γv * d_term
end
