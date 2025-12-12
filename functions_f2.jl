### Functions and derivatives for  f₂
function f₂(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    r = sqrt(κ)
    a = θ + r * q₁(κ, μ)
    b = θ - r * q₂(κ, μ)
    term = inv(sqr(sin(a))) + μ * inv(sqr(sin(b)))
    return sin(θ) * cos(θ) - Γ(κ, μ) * term
end

function f₂θ(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    r  = sqrt(κ)
    q1 = q₁(κ, μ)
    q2 = one(q1) - q1
    Γv = Γ(κ, μ)

    a = θ + r * q1
    b = θ - r * q2

    sa, ca = sin(a), cos(a)
    sb, cb = sin(b), cos(b)

    inv_sa3 = inv(sa * sa * sa)
    inv_sb3 = inv(sb * sb * sb)

    d_main = cos(twice(θ))
    d_term = (twice(ca) * inv_sa3) + (twice(μ) * cb * inv_sb3)

    return d_main + Γv * d_term
end


function f₂θθ(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    oneT = one(θ)
    two = twice(oneT)

    r = sqrt(κ)
    q1 = q₁(κ, μ)
    q2 = q₂(κ, μ)
    Γv = Γ(κ, μ)

    A = θ + r * q1
    B = θ - r * q2

    sA = sin(A)
    cA = cos(A)
    sB = sin(B)
    cB = cos(B)

    main = -two * sin(twice(θ))

    inv_sA4 = inv(sA * sA * sA * sA)
    inv_sB4 = inv(sB * sB * sB * sB)

    two_cA2 = two * cA * cA
    two_cB2 = two * cB * cB

    FA′ = -(oneT + two_cA2) * inv_sA4
    FB′ = -(oneT + two_cB2) * inv_sB4

    return main + two * Γv * (FA′ + μ * FB′)
end

function f₂κ(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    oneT = one(κ); twoT = twice(oneT); fourT = twoT + twoT

    r  = sqrt(κ)
    @assert inf(r) > zero(r) "κ interval must be > 0 for ∂f₂/∂κ (1/(2r) present)."

    q1 = q₁(κ, μ)
    q2 = one(q1) - q1
    G  = Γ(κ, μ)

    inv_2r = one(r) / (twoT * r)          

    s2r = sin(twice(r))
    c2r = cos(twice(r))
    y   = μ * s2r
    x   = one(μ) + μ * c2r

    y′   = μ * twice(c2r) * inv_2r       
    x′   = -μ * twice(s2r) * inv_2r
    denom = x*x + y*y
    atan′ = (x*y′ - y*x′) / denom

    g      = inv_2r                        
    g_over = - one(r) / (twoT * r * r)      
    dq1_dκ = g_over * q1 + g * atan′


    sr, cr = sin(r), cos(r)
    num = sr^3 * cr
    den = sqrt( μ*μ + twoT*(μ * c2r) + one(μ) )

    dnum = ( (oneT+oneT+oneT)*sr^2*cr^2 - sr^4 ) * inv_2r  
    dden = ( (-fourT * μ * s2r) * inv_2r ) / (twoT * den)  

    dG = (dnum*den - num*dden) / (den*den)

    a = θ + r*q1
    b = θ - r*q2

    sa, ca = sin(a), cos(a)
    sb, cb = sin(b), cos(b)

    inv_sa2 = inv(sa*sa)
    inv_sb2 = inv(sb*sb)
    inv_sa3 = inv(sa*sa*sa)
    inv_sb3 = inv(sb*sb*sb)

    aκ = inv_2r*q1 + r*dq1_dκ
    bκ = -inv_2r*q2 + r*dq1_dκ

    d_term = (twoT*ca*inv_sa3)*aκ + (twoT*μ*cb*inv_sb3)*bκ
    term   = inv_sa2 + μ*inv_sb2

    return -dG*term + G*d_term
end

function f₂μ(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)

    oneT = one(μ)
    twoT = twice(oneT)

    r  = sqrt(κ)
    q1 = q₁(κ, μ)
    q2 = one(q1) - q1
    G  = Γ(κ, μ)

    a = θ + r*q1
    b = θ - r*q2

    sa, ca = sin(a), cos(a)
    sb, cb = sin(b), cos(b)

    inv_sa2 = inv(sa*sa)
    inv_sb2 = inv(sb*sb)
    inv_sa3 = inv(sa*sa*sa)
    inv_sb3 = inv(sb*sb*sb)

    # Γ_μ: Γ = num/den, num = sin(r)^3 cos(r) (no μ), den = sqrt(μ^2 + 2μ cos(2r) + 1)
    s2r = sin(twice(r))
    c2r = cos(twice(r))
    sr, cr = sin(r), cos(r)
    num = sr^3 * cr
    den = sqrt( μ*μ + twoT*(μ * c2r) + one(μ) )

    # den_μ = (μ + cos(2r)) / den  ⇒  Γ_μ = - Γ * (μ + cos(2r)) / den^2
    Gμ = - G * (μ + c2r) / (den*den)

    # q₁_μ: q1 = (1/(2r)) * atan(y,x), with y = μ sin(2r), x = 1 + μ cos(2r)
    y  = μ * s2r
    x  = one(μ) + μ * c2r
    yμ = s2r
    xμ = c2r
    atanμ = (x*yμ - y*xμ) / (x*x + y*y)
    inv_2r = one(r) / (twoT * r)
    q1μ = inv_2r * atanμ

    # a_μ = r*q1_μ,  b_μ = r*q1_μ  (since q₂ = 1 - q₁)
    aμ = r * q1μ
    bμ = r * q1μ

    # d/dμ csc²(a) = -2 cos(a)/sin^3(a) * a_μ  ;  d/dμ( μ csc²(b) ) = csc²(b) + μ * (-2 cos(b)/sin^3(b)) * b_μ
    d_csc2_a = - twoT * ca * inv_sa3 * aμ
    d_mu_csc2_b = inv_sb2 + ( - twoT * μ * cb * inv_sb3 * bμ )

    term = inv_sa2 + μ*inv_sb2

    # assemble ∂f₂/∂μ
    return - Gμ * term - G * ( d_csc2_a + d_mu_csc2_b )
end


### Intervals and brackets for f₁
function Intθ_2(κ, μ)
    return (sqrt(κ) * q₂(κ, μ), pi - sqrt(κ) * q₁(κ, μ))
end


function brackets_f₂(κ, μ; n::Integer=400, margin::Real=0.0)
    θmin = Intθ_2(κ, μ)[1]
    θmax = Intθ_2(κ, μ)[2]
    δ = margin * (θmax - θmin)
    a, b = θmin + δ, θmax - δ
    brackets_signchanges(θ -> f₂(θ, κ, μ), a, b; n=n)
end