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
    sκ1 = sin(r*q1)
    sκ2 = sin(r*q2)

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