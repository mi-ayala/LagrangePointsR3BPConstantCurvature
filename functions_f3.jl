function f₃(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    r = sqrt(κ)
    a = θ + r * q₁(κ, μ)
    b = θ - r * q₂(κ, μ)
    term = inv(sqr(sin(a))) + μ * inv(sqr(sin(b)))
    return sin(θ) * cos(θ) + Γ(κ, μ) * term
end

function f₃θ(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    r  = sqrt(κ)
    q1 = q₁(κ, μ)
    q2 = q₂(κ, μ)
    Γv = Γ(κ, μ)

    a = θ + r * q1
    b = θ - r * q2

    # derivative of sinθ cosθ
    d_main = cos(twice(θ))      # = cos(2θ)

    # pieces for the Γ-term
    sa = sin(a); ca = cos(a)
    sb = sin(b); cb = cos(b)

    inv_sa3 = inv(sa * sa * sa) # 1 / sin(a)^3
    inv_sb3 = inv(sb * sb * sb) # 1 / sin(b)^3

    # d/dθ[ sin(a)^(-2) ] = -2 cos(a) / sin(a)^3
    # d/dθ[ μ sin(b)^(-2) ] = -2 μ cos(b) / sin(b)^3
    d_term = (-twice(ca) * inv_sa3) + (-twice(μ) * cb * inv_sb3)

    return d_main + Γv * d_term
end

function f₃θθ(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    r  = sqrt(κ)
    q1 = q₁(κ, μ)
    q2 = q₂(κ, μ)
    Γv = Γ(κ, μ)

    a = θ + r * q1
    b = θ - r * q2

    # d²/dθ² [sinθ cosθ] = d/dθ[cos(2θ)] = -2 sin(2θ)
    d2_main = -twice(sin(twice(θ)))

    sa = sin(a); ca = cos(a)
    sb = sin(b); cb = cos(b)

    inv_sa2 = inv(sa * sa)          # 1/sin(a)^2
    inv_sb2 = inv(sb * sb)          # 1/sin(b)^2
    inv_sa4 = inv_sa2 * inv_sa2     # 1/sin(a)^4
    inv_sb4 = inv_sb2 * inv_sb2     # 1/sin(b)^4

    two   = twice(one(θ))
    three = one(θ) + two
    six   = two * three

    # d²/dθ²[ sin(a)^(-2) ] = 2/sin(a)^2 + 6 cos(a)^2 / sin(a)^4
    # d²/dθ²[ μ sin(b)^(-2) ] = μ*(2/sin(b)^2 + 6 cos(b)^2 / sin(b)^4)
    d2_term =
        (two * inv_sa2 + six * (ca * ca) * inv_sa4) +
        μ * (two * inv_sb2 + six * (cb * cb) * inv_sb4)

    return d2_main + Γv * d2_term
end

function Intθ_3(κ, μ)
    return (sqrt(κ) * q₂(κ, μ) - pi, -sqrt(κ) * q₁(κ, μ))
end

function brackets_f₃(κ, μ; n::Integer=400, margin::Real=0.0)
    θmin = Intθ_3(κ, μ)[1]
    θmax = Intθ_3(κ, μ)[2]
    δ = margin * (θmax - θmin)
    a, b = θmin + δ, θmax - δ
    brackets_signchanges(θ -> f₃(θ, κ, μ), a, b; n=n)
end