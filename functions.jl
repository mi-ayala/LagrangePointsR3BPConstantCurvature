# helpers (safe for Reals and Intervals)
@inline twice(x) = x + x
@inline csc(x) = one(x) / sin(x)
@inline powabs3(x) = abs(x)^3
@inline sqr(x) = x * x


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

### Functions and interval
function f₁(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    r = sqrt(κ)
    a = θ + r * q₁(κ, μ)
    b = θ - r * q₂(κ, μ)
    term = inv(sqr(sin(a))) - μ * inv(sqr(sin(b)))
    return sin(θ) * cos(θ) - Γ(κ, μ) * term
end

function Intθ_1(κ, μ)
    return (-sqrt(κ) * q₁(κ, μ), sqrt(κ) * q₂(κ, μ))
end

function θ1_max(κ, μ)
    κ, μ = promote(κ, μ)
    r = sqrt(κ)
    y = tan(r * q₂(κ, μ))
    x = one(μ) + μ
    return atan(y, x)
end


# left endpoint of the above interval
function θ4_min(κ, μ)
    κ, μ = promote(κ, μ)
    r = sqrt(κ)
    q1 = q₁(κ, μ)
    q2 = one(q1) - q1
    G = Γ(κ, μ)

    πT = π * one(κ)                    #

    # sqrt( 4Γ / (2κ q2^2 + π^2 Γ μ) )
    oneG = one(G)
    two = oneG + oneG
    four = two + two
    denom = two * κ * (q2^2) + (πT^2) * G * μ
    rad = (four * G) / denom
    s = sqrt(rad)

    ones = one(s)                      
    return πT - (ones - s) * r * q1
end

# right endpoint of the interval
function θ4_max(κ, μ)
    κ, μ = promote(κ, μ)
    r = sqrt(κ)
    q1 = q₁(κ, μ)
    q2 = one(q1) - q1
    G = Γ(κ, μ)

    πT = π * one(κ)

    # sqrt( 4Γ μ / (2κ q1^2 + π^2 Γ) )
    oneG = one(G)
    two = oneG + oneG
    four = two + two
    denom = two * κ * q1^2 + (πT^2) * G
    rad = (four * G * μ) / denom
    s = sqrt(rad)

    ones = one(s)
    return πT + (ones - (q1 / q2) * s) * r * q2
end


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

function f₂κ(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    oneT = one(κ); twoT = twice(oneT); fourT = twoT + twoT

    r  = sqrt(κ)
    @assert inf(r) > zero(r) "κ interval must be > 0 for ∂f₂/∂κ (1/(2r) present)."

    q1 = q₁(κ, μ)
    q2 = one(q1) - q1
    G  = Γ(κ, μ)

    inv_2r = one(r) / (twoT * r)           # 1/(2r)

    # --- dq1/dκ using q1 (avoid fresh atan) ---
    s2r = sin(twice(r))
    c2r = cos(twice(r))
    y   = μ * s2r
    x   = one(μ) + μ * c2r

    y′   = μ * twice(c2r) * inv_2r         # μ*(2 cos(2r))*(1/(2r))
    x′   = -μ * twice(s2r) * inv_2r
    denom = x*x + y*y
    atan′ = (x*y′ - y*x′) / denom

    g      = inv_2r                         # 1/(2r)
    g_over = - one(r) / (twoT * r * r)      # (g′/g) = -1/(2 r^2)
    dq1_dκ = g_over * q1 + g * atan′

    # --- dΓ/dκ ---
    sr, cr = sin(r), cos(r)
    num = sr^3 * cr
    den = sqrt( μ*μ + twoT*(μ * c2r) + one(μ) )

    dnum = ( (oneT+oneT+oneT)*sr^2*cr^2 - sr^4 ) * inv_2r  # (3 s^2 c^2 - s^4)/(2r)
    dden = ( (-fourT * μ * s2r) * inv_2r ) / (twoT * den)  # (du/dr * dr/dκ)/(2√u)

    dG = (dnum*den - num*dden) / (den*den)

    # --- angles, terms ---
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

    # typed “constants”
    oneT = one(μ)
    twoT = twice(oneT)

    r  = sqrt(κ)
    q1 = q₁(κ, μ)
    q2 = one(q1) - q1
    G  = Γ(κ, μ)

    # angles
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


function Intθ_2(κ, μ)
    return (sqrt(κ) * q₂(κ, μ), pi - sqrt(κ) * q₁(κ, μ))
end


function f₃(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    r = sqrt(κ)
    a = θ + r * q₁(κ, μ)
    b = θ - r * q₂(κ, μ)
    term = inv(sqr(sin(a))) + μ * inv(sqr(sin(b)))
    return sin(θ) * cos(θ) + Γ(κ, μ) * term
end

function Intθ_3(κ, μ)
    return (sqrt(κ) * q₂(κ, μ) - pi, -sqrt(κ) * q₁(κ, μ))
end

function f₄(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    r = sqrt(κ)
    a = θ + r * q₁(κ, μ)
    b = θ - r * q₂(κ, μ)
    term = inv(sqr(sin(a))) - μ * inv(sqr(sin(b)))
    return sin(θ) * cos(θ) + Γ(κ, μ) * term
end

function Intθ_4(κ, μ)
    return (pi - sqrt(κ) * q₁(κ, μ), pi + sqrt(κ) * q₂(κ, μ))
end

### Derivatives
function f₁θ(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    r = sqrt(κ)
    a = θ + r * q₁(κ, μ)
    b = θ - r * q₂(κ, μ)
    Γv = Γ(κ, μ)

    d_main = cos(twice(θ))

    sa = sin(a)
    ca = cos(a)
    sb = sin(b)
    cb = cos(b)

    inv_sa3 = inv(sa * sa * sa)
    inv_sb3 = inv(sb * sb * sb)
    d_term = (-twice(ca) * inv_sa3) + (twice(μ) * cb * inv_sb3)

    return d_main - Γv * d_term
end


# f₄θ(θ, κ, μ) = ∂f₄/∂θ
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


function brackets_f₁(κ, μ; n::Integer=400, margin::Real=0.0)
    θmin = -sqrt(κ) * q₁(κ, μ)
    θmax = sqrt(κ) * q₂(κ, μ)
    δ = margin * (θmax - θmin)
    a, b = θmin + δ, θmax - δ
    brackets_signchanges(θ -> f₁(θ, κ, μ), a, b; n=n)
end

function brackets_f₂(κ, μ; n::Integer=400, margin::Real=0.0)
    θmin = Intθ_2(κ, μ)[1]
    θmax = Intθ_2(κ, μ)[2]
    δ = margin * (θmax - θmin)
    a, b = θmin + δ, θmax - δ
    brackets_signchanges(θ -> f₂(θ, κ, μ), a, b; n=n)
end

function brackets_f₃(κ, μ; n::Integer=400, margin::Real=0.0)
    θmin = Intθ_3(κ, μ)[1]
    θmax = Intθ_3(κ, μ)[2]
    δ = margin * (θmax - θmin)
    a, b = θmin + δ, θmax - δ
    brackets_signchanges(θ -> f₃(θ, κ, μ), a, b; n=n)
end

function brackets_f₄(κ, μ; n::Integer=400, margin::Real=0.0)
    θmin = Intθ_4(κ, μ)[1]
    θmax = Intθ_4(κ, μ)[2]
    δ = margin * (θmax - θmin)
    a, b = θmin + δ, θmax - δ
    brackets_signchanges(θ -> f₄(θ, κ, μ), a, b; n=n)
end

function brackets_signchanges(f, a::Real, b::Real; n::Integer=400)
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
            continue  # skip nonfinite; keep previous finite fa/sa
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


# ---- Notable functions ----
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
    sκ = sin(r)

    return -G * sin(θ) * (sκ / powabs3(s1) - μ * (sκ / powabs3(s2)))
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

    two = one(θ) + one(θ)      # “2” with the same type as θ (Float64, Interval, BigFloat,…)
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


"""
    classify_equilibrium_fp(θ, κ, μ; verbose=false)

Floating-point version of the interval classifier.

Implements:

1) detλ(θ) < 0                → :saddle_hyperbolic
2) λ₂(θ)  > 0                 → :lyapunov_stable
3) Δ(θ)   > 0                 → :elliptic
4) Δ(θ)   < 0                 → :linearly_unstable
5) otherwise → :unknown

Returns `(label::Symbol, det_val, lam2_val, Δ_val)`.
"""
function classify_equilibrium_fp(θ, κ, μ; verbose::Bool=false)
    det_val = detλ(θ, κ, μ)

    if verbose
        println("det(θ) = ", det_val)
    end
    if det_val < 0
        return :saddle_hyperbolic, det_val, NaN, NaN
    end

    lam2_val = λ₂(θ, κ, μ)
    if verbose
        println("λ₂(θ) = ", lam2_val)
    end
    if lam2_val > 0
        return :lyapunov_stable, det_val, lam2_val, NaN
    end

    Δ_val = Δ(θ, κ, μ)
    if verbose
        println("Δ(θ) = ", Δ_val)
    end
    if Δ_val > 0
        return :elliptic, det_val, lam2_val, Δ_val
    elseif Δ_val < 0
        return :linearly_unstable, det_val, lam2_val, Δ_val
    else
        return :unknown, det_val, lam2_val, Δ_val
    end
end