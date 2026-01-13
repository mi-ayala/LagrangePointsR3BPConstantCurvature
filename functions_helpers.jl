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


# λ₁(θ) = -[ cos(2θ) + 2Γ(  cos(θ+√κ q₁)/|sin(θ+√κ q₁)|^3  + μ cos(θ-√κ q₂)/|sin(θ-√κ q₂)|^3 ) ]
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

# λ₂(θ) = -Γ sinθ [ sin(√κ)/|sin(θ+√κ q₁)|^3  - μ sin(√κ)/|sin(θ-√κ q₂)|^3 ]
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

# det(θ) = λ₁(θ) * λ₂(θ)
detλ(θ, κ, μ) = λ₁(θ, κ, μ) * λ₂(θ, κ, μ)


# c₂(θ) = 2 + cos(2θ) - Γ csc(θ) [ sin(2θ+√κ q₁)/|sin(θ+√κ q₁)|^3 - μ sin(2θ-√κ q₂)/|sin(θ-√κ q₂)|^3 ]
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

# Δ(θ) = c₂(θ)^2 - 4*det(θ)/sin^2(θ)
function Δ(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    c2 = c₂(θ, κ, μ)
    four = 4 * one(c2)
    sθ = sin(θ)
    two = one(θ) + one(θ)
    four = two + two
    return c2^2 - (four * detλ(θ, κ, μ)) / (sθ^2)
end


function brackets_signchanges(f, a::Real, b::Real; n::Integer=200)
    xs = collect(range(a, b; length=n + 1))
    br = Tuple{Float64,Float64}[]

    # evaluate first finite value
    i = 1
    fa = f(xs[i])
    while i <= n && !isfinite(fa)
        i += 1
        fa = f(xs[i])
    end
    i > n && return br  # nothing finite

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


