### Functions and derivatives for f₃ with κ < 0
function f₃neg(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    r = sqrt(-κ)
    a = θ + r * q₁neg(κ, μ)
    b = θ - r * q₂neg(κ, μ)
    term = inv(sqr(sinh(a))) + μ * inv(sqr(sinh(b)))
    return sinh(θ) * cosh(θ) + Γneg(κ, μ) * term
end

function f₃negθ(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    r = sqrt(-κ)
    q1 = q₁neg(κ, μ)
    q2 = q₂neg(κ, μ)
    Γv = Γneg(κ, μ)

    a = θ + r * q1
    b = θ - r * q2

    d_main = cosh(twice(θ))

    sa = sinh(a)
    ca = cosh(a)
    sb = sinh(b)
    cb = cosh(b)

    inv_sa3 = inv(sa * sa * sa)
    inv_sb3 = inv(sb * sb * sb)

    d_term = (-twice(ca) * inv_sa3) + (-twice(μ) * cb * inv_sb3)

    return d_main + Γv * d_term
end

function f₃negθθ(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    oneT = one(θ)
    two = twice(oneT)

    r = sqrt(-κ)
    q1 = q₁neg(κ, μ)
    q2 = q₂neg(κ, μ)
    Γv = Γneg(κ, μ)

    a = θ + r * q1
    b = θ - r * q2

    d2_main = two * sinh(twice(θ))

    sa = sinh(a)
    ca = cosh(a)
    sb = sinh(b)
    cb = cosh(b)

    inv_sa4 = inv(sa * sa * sa * sa)
    inv_sb4 = inv(sb * sb * sb * sb)

    d2_term =
        (two * (oneT + two * ca * ca) * inv_sa4) +
        μ * (two * (oneT + two * cb * cb) * inv_sb4)

    return d2_main + Γv * d2_term
end
