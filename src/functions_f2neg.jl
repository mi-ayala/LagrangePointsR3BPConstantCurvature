### Functions and derivatives for f₂ with κ < 0
function f₂neg(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    r = sqrt(-κ)
    a = θ + r * q₁neg(κ, μ)
    b = θ - r * q₂neg(κ, μ)
    term = inv(sqr(sinh(a))) + μ * inv(sqr(sinh(b)))
    return sinh(θ) * cosh(θ) - Γneg(κ, μ) * term
end

function f₂negθ(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    r = sqrt(-κ)
    q1 = q₁neg(κ, μ)
    q2 = one(q1) - q1
    Γv = Γneg(κ, μ)

    a = θ + r * q1
    b = θ - r * q2

    sa, ca = sinh(a), cosh(a)
    sb, cb = sinh(b), cosh(b)

    inv_sa3 = inv(sa * sa * sa)
    inv_sb3 = inv(sb * sb * sb)

    d_main = cosh(twice(θ))
    d_term = (twice(ca) * inv_sa3) + (twice(μ) * cb * inv_sb3)

    return d_main + Γv * d_term
end

function f₂negθθ(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    oneT = one(θ)
    two = twice(oneT)

    r = sqrt(-κ)
    q1 = q₁neg(κ, μ)
    q2 = q₂neg(κ, μ)
    Γv = Γneg(κ, μ)

    A = θ + r * q1
    B = θ - r * q2

    sA = sinh(A)
    cA = cosh(A)
    sB = sinh(B)
    cB = cosh(B)

    main = two * sinh(twice(θ))

    inv_sA4 = inv(sA * sA * sA * sA)
    inv_sB4 = inv(sB * sB * sB * sB)

    two_cA2 = two * cA * cA
    two_cB2 = two * cB * cB

    FA′ = -(oneT + two_cA2) * inv_sA4
    FB′ = -(oneT + two_cB2) * inv_sB4

    return main + two * Γv * (FA′ + μ * FB′)
end

function f₂negκ(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    oneT = one(κ)
    twoT = twice(oneT)
    fourT = twoT + twoT

    r = sqrt(-κ)
    @assert inf(r) > zero(r) "(-κ) interval must be > 0 for ∂f₂neg/∂κ (1/(2r) present)."

    q1 = q₁neg(κ, μ)
    q2 = one(q1) - q1
    G = Γneg(κ, μ)

    inv_2r = one(r) / (twoT * r)

    s2r = sinh(twice(r))
    c2r = cosh(twice(r))
    y = μ * s2r
    x = one(μ) + μ * c2r

    yκ = -μ * twice(c2r) * inv_2r
    xκ = -μ * twice(s2r) * inv_2r
    denom = x * x - y * y
    atanhκ = (x * yκ - y * xκ) / denom

    g = inv_2r
    g_over = one(r) / (twoT * r * r)
    dq1_dκ = g_over * q1 + g * atanhκ

    sr, cr = sinh(r), cosh(r)
    num = sr^3 * cr
    den = sqrt(μ * μ + twoT * (μ * c2r) + one(μ))

    dnum = -((oneT + oneT + oneT) * sr^2 * cr^2 + sr^4) * inv_2r
    dden = ((-fourT * μ * s2r) * inv_2r) / (twoT * den)

    dG = (dnum * den - num * dden) / (den * den)

    a = θ + r * q1
    b = θ - r * q2

    sa, ca = sinh(a), cosh(a)
    sb, cb = sinh(b), cosh(b)

    inv_sa2 = inv(sa * sa)
    inv_sb2 = inv(sb * sb)
    inv_sa3 = inv(sa * sa * sa)
    inv_sb3 = inv(sb * sb * sb)

    aκ = -inv_2r * q1 + r * dq1_dκ
    bκ = inv_2r * q2 + r * dq1_dκ

    d_term = (-twoT * ca * inv_sa3) * aκ + (-twoT * μ * cb * inv_sb3) * bκ
    term = inv_sa2 + μ * inv_sb2

    return -dG * term - G * d_term
end

function f₂negμ(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)

    oneT = one(μ)
    twoT = twice(oneT)

    r = sqrt(-κ)
    q1 = q₁neg(κ, μ)
    q2 = one(q1) - q1
    G = Γneg(κ, μ)

    a = θ + r * q1
    b = θ - r * q2

    sa, ca = sinh(a), cosh(a)
    sb, cb = sinh(b), cosh(b)

    inv_sa2 = inv(sa * sa)
    inv_sb2 = inv(sb * sb)
    inv_sa3 = inv(sa * sa * sa)
    inv_sb3 = inv(sb * sb * sb)

    s2r = sinh(twice(r))
    c2r = cosh(twice(r))
    den = sqrt(μ * μ + twoT * (μ * c2r) + one(μ))

    Gμ = -G * (μ + c2r) / (den * den)

    y = μ * s2r
    x = one(μ) + μ * c2r
    yμ = s2r
    xμ = c2r
    atanhμ = (x * yμ - y * xμ) / (x * x - y * y)
    inv_2r = one(r) / (twoT * r)
    q1μ = inv_2r * atanhμ

    aμ = r * q1μ
    bμ = r * q1μ

    d_csch2_a = -twoT * ca * inv_sa3 * aμ
    d_mu_csch2_b = inv_sb2 + (-twoT * μ * cb * inv_sb3 * bμ)

    term = inv_sa2 + μ * inv_sb2

    return -Gμ * term - G * (d_csch2_a + d_mu_csch2_b)
end
