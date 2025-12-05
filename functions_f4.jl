function f₄(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    r = sqrt(κ)
    a = θ + r * q₁(κ, μ)
    b = θ - r * q₂(κ, μ)
    term = inv(sqr(sin(a))) - μ * inv(sqr(sin(b)))
    return sin(θ) * cos(θ) + Γ(κ, μ) * term
end

function f₄θ(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)

    r = sqrt(κ)
    q1 = q₁(κ, μ)
    q2 = one(q1) - q1
    Γv = Γ(κ, μ)

    a = θ + r * q1
    b = θ - r * q2

    d_main = cos(twice(θ))

    sa = sin(a)
    ca = cos(a)
    sb = sin(b)
    cb = cos(b)

    inv_sa3 = inv(sa * sa * sa)
    inv_sb3 = inv(sb * sb * sb)

    d_term = (-twice(ca) * inv_sa3) + (twice(μ) * cb * inv_sb3)

    return d_main + Γv * d_term   # note the + (for f₁θ it was -)
end



function Intθ_4(κ, μ)
    return (pi - sqrt(κ) * q₁(κ, μ), pi + sqrt(κ) * q₂(κ, μ))
end

function θ4_min(κ, μ)
    κ, μ = promote(κ, μ)
    r = sqrt(κ)
    q1 = q₁(κ, μ)
    q2 = one(q1) - q1
    G = Γ(κ, μ)

    πT = π * one(κ)                    

    oneG = one(G)
    two = oneG + oneG
    four = two + two
    denom = two * κ * (q2^2) + (πT^2) * G * μ
    rad = (four * G) / denom
    s = sqrt(rad)

    ones = one(s)                      
    return πT - (ones - s) * r * q1
end

function θ4_max(κ, μ)
    κ, μ = promote(κ, μ)
    r = sqrt(κ)
    q1 = q₁(κ, μ)
    q2 = one(q1) - q1
    G = Γ(κ, μ)

    πT = π * one(κ)

    oneG = one(G)
    two = oneG + oneG
    four = two + two
    denom = two * κ * q1^2 + (πT^2) * G
    rad = (four * G * μ) / denom
    s = sqrt(rad)

    ones = one(s)
    return πT + (ones - (q1 / q2) * s) * r * q2
end

function brackets_f₄(κ, μ; n::Integer=400, margin::Real=0.0)
    θmin = Intθ_4(κ, μ)[1]
    θmax = Intθ_4(κ, μ)[2]
    δ = margin * (θmax - θmin)
    a, b = θmin + δ, θmax - δ
    brackets_signchanges(θ -> f₄(θ, κ, μ), a, b; n=n)
end