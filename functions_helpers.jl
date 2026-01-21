#### Helper functions that we use in multiple places

@inline twice(x) = x + x
@inline csc(x) = one(x) / sin(x)
@inline powabs3(x) = abs(x)^3
@inline sqr(x) = x * x
@inline three(x) = x + x + x
@inline inv3(x) = inv(three(one(x)))
@inline cbrt_t(x) = x^inv3(x)

function q₁(κ, μ)
    κ, μ = promote(κ, μ)
    r = sqrt(κ)
    y = μ * sin(twice(r))
    x = one(μ) + μ * cos(twice(r))
    return (one(r) / twice(r)) * atan(y, x)
end

function q₂(κ, μ)
    κ, μ = promote(κ, μ)
    return one(κ) - q₁(κ, μ)
end

function Γ(κ, μ)
    κ, μ = promote(κ, μ)
    r = sqrt(κ)
    num = sqr(sin(r)) * sin(r) * cos(r)
    den = sqrt(μ * μ + twice(μ * cos(twice(r))) + one(μ))
    return num / den
end

function λ₁(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    r = sqrt(κ)
    q1 = q₁(κ, μ)
    q2 = one(q1) - q1
    G = Γ(κ, μ)

    s1 = sin(θ + r * q1)
    s2 = sin(θ - r * q2)
    t1 = cos(θ + r * q1) / powabs3(s1)
    t2 = cos(θ - r * q2) / powabs3(s2)

    return -(cos(twice(θ)) + twice(G) * (t1 + μ * t2))
end

function λ₂(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    r = sqrt(κ)
    q1 = q₁(κ, μ)
    q2 = one(q1) - q1
    G = Γ(κ, μ)

    s1 = sin(θ + r * q1)
    s2 = sin(θ - r * q2)
    sκ1 = sin(r * q1)
    sκ2 = sin(r * q2)

    return -G * sin(θ) * (sκ1 / powabs3(s1) - μ * (sκ2 / powabs3(s2)))
end

detλ(θ, κ, μ) = λ₁(θ, κ, μ) * λ₂(θ, κ, μ)

function c₂(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    r = sqrt(κ)
    q1 = q₁(κ, μ)
    q2 = one(q1) - q1
    G = Γ(κ, μ)

    s1 = sin(θ + r * q1)
    s2 = sin(θ - r * q2)
    a = sin(twice(θ) + r * q1) / powabs3(s1)
    b = sin(twice(θ) - r * q2) / powabs3(s2)

    two = one(θ) + one(θ)
    return two + cos(twice(θ)) - G * csc(θ) * (a - μ * b)
end


function Δ(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    c2 = c₂(θ, κ, μ)
    four = 4 * one(c2)
    sθ = sin(θ)
    two = one(θ) + one(θ)
    four = two + two
    return c2^2 - (four * detλ(θ, κ, μ)) / (sθ^2)
end


function q₁μ(κ, μ)
    κ, μ = promote(κ, μ)
    r = sqrt(κ)

    oneT = one(r)
    two = twice(oneT)

    s2r = sin(twice(r))
    c2r = cos(twice(r))

    y = μ * s2r
    x = one(μ) + μ * c2r

    yμ = s2r
    xμ = c2r

    denom = x * x + y * y
    atanμ = (x * yμ - y * xμ) / denom

    inv_2r = one(r) / (two * r)

    return inv_2r * atanμ
end

function q₁κ(κ, μ)
    κ, μ = promote(κ, μ)
    r = sqrt(κ)

    oneT = one(r)
    two = twice(oneT)

    q1 = q₁(κ, μ)

    dr_dκ = one(r) / (two * r)   # 1/(2r)

    s2r = sin(twice(r))
    c2r = cos(twice(r))

    y = μ * s2r
    x = one(μ) + μ * c2r

    yκ = μ * twice(c2r) * dr_dκ
    xκ = -μ * twice(s2r) * dr_dκ

    denom = x * x + y * y
    atanκ = (x * yκ - y * xκ) / denom

    g = one(r) / (two * r)
    g_over = -one(r) / (two * r * r)     # g'/g = -1/(2r²)

    return g_over * q1 + g * atanκ
end


function Γμ(κ, μ)
    κ, μ = promote(κ, μ)
    r = sqrt(κ)

    oneT = one(κ)
    two = twice(oneT)

    s2r = sin(twice(r))
    c2r = cos(twice(r))

    sr, cr = sin(r), cos(r)

    num = sr * sr * sr * cr
    den = sqrt(μ * μ + two * (μ * c2r) + one(μ))

    G = num / den
    return -G * (μ + c2r) / (den * den)
end

function Γκ(κ, μ)
    κ, μ = promote(κ, μ)
    r = sqrt(κ)

    oneT = one(κ)
    two = twice(oneT)
    four = two + two

    dr_dκ = one(r) / (two * r)

    s2r = sin(twice(r))
    c2r = cos(twice(r))

    sr, cr = sin(r), cos(r)

    num = sr * sr * sr * cr
    den = sqrt(μ * μ + two * (μ * c2r) + one(μ))

    three = two + oneT

    dnum = (three * sr * sr * cr * cr - sr * sr * sr * sr) * dr_dκ

    uκ = -four * μ * s2r * dr_dκ
    dden = uκ / (two * den)

    return (dnum * den - num * dden) / (den * den)
end

function Λ(κ, μ)
    κ, μ = promote(κ, μ)
    r = sqrt(κ)
    q1 = q₁(κ, μ)
    q2 = q₂(κ, μ)

    s1 = cos(r * q1)
    s2 = cos(r * q2)

    return cbrt_t(s2 / s1)
end

function Λκ(κ, μ)
    κ, μ = promote(κ, μ)

    r = sqrt(κ)
    q1 = q₁(κ, μ)
    q2 = q₂(κ, μ)

    q1κv = q₁κ(κ, μ)
    q2κv = -q1κv

    oneκ = one(κ)
    two = oneκ + oneκ

    rκ = inv(two * r)   # dr/dκ

    A = cos(r * q2)
    B = cos(r * q1)

    Aκ = -sin(r * q2) * (rκ * q2 + r * q2κv)
    Bκ = -sin(r * q1) * (rκ * q1 + r * q1κv)

    U = A / B
    Uκ = (Aκ * B - A * Bκ) / (B * B)

    return inv3(U) * inv(cbrt_t(U)^2) * Uκ
end

function Λμ(κ, μ)
    κ, μ = promote(κ, μ)

    r = sqrt(κ)
    q1 = q₁(κ, μ)
    q2 = q₂(κ, μ)

    q1μv = q₁μ(κ, μ)
    q2μv = -q1μv

    A = cos(r * q2)
    B = cos(r * q1)

    Aμ = -sin(r * q2) * (r * q2μv)
    Bμ = -sin(r * q1) * (r * q1μv)

    U = A / B
    Uμ = (Aμ * B - A * Bμ) / (B * B)

    return inv3(U) * inv(cbrt_t(U)^2) * Uμ
end

function brackets_signchanges(f, a::Real, b::Real; n::Integer=200)
    xs = collect(range(a, b; length=n + 1))
    br = Tuple{Float64,Float64}[]

    i = 1
    fa = f(xs[i])
    while i <= n && !isfinite(fa)
        i += 1
        fa = f(xs[i])
    end
    i > n && return br

    sgn(x) = x > 0 ? 1 : x < 0 ? -1 : 0
    sa = sgn(fa)

    for j in i:n
        x1, x2 = xs[j], xs[j+1]
        fx = f(x2)
        if !isfinite(fx)
            continue
        end
        sb = sgn(fx)

        if sb == 0
            push!(br, (x2, x2))              # exact root at node
        elseif sa == 0
            push!(br, (x1, x1))              # exact root at previous node
        elseif sa != sb
            push!(br, (x1, x2))              # sign flip bracket
        end

        fa, sa = fx, sb                      # advance
    end
    return br
end


