using RadiiPolynomial

include("functions_helpers.jl")

function F_appendix_D(x)
    θ, κ, μ = x
    θ, κ, μ = promote(θ, κ, μ)

    r = sqrt(κ)
    q1 = q₁(κ, μ)
    q2 = q₂(κ, μ)
    L = Λ(κ, μ)
    G = Γ(κ, μ)

    oneθ = one(θ)
    two = oneθ + oneθ

    s = sin(θ + r * q1)

    F1 = s - L * sin(θ - r * q2)

    F2 = sin(two * θ) * (s * s) + (two * G) * (oneθ + μ * (L * L))

    F3 = cos(two * θ) * (s * s * s) * sin(r * q2) - (two * G) * cos(θ) * sin(r)

    return [F1, F2, F3]
end

function J_equilibrium(x)
    θ, κ, μ = x
    θ, κ, μ = promote(θ, κ, μ)
    T = typeof(θ)

    # core scalars
    r = sqrt(κ)
    q1 = q₁(κ, μ)
    q2 = q₂(κ, μ)
    L = Λ(κ, μ)
    G = Γ(κ, μ)

    # interval-safe constants
    oneθ = one(θ)
    two = oneθ + oneθ
    three = two + oneθ

    # derivatives of r
    rκ = inv(two * r)      # dr/dκ = 1/(2√κ)

    # derivatives of q1,q2, L, G
    q1κv = q₁κ(κ, μ)
    q1μv = q₁μ(κ, μ)
    q2κv = -q1κv
    q2μv = -q1μv

    Lκv = Λκ(κ, μ)
    Lμv = Λμ(κ, μ)
    Gκv = Γκ(κ, μ)
    Gμv = Γμ(κ, μ)

    # phase variables
    a = θ + r * q1           # θ + √κ q1
    b = θ - r * q2           # θ - √κ q2
    c = r * q2               # √κ q2
    d = r                  # √κ

    s = sin(a)
    ca = cos(a)
    sb = sin(b)
    cb = cos(b)
    sc = sin(c)
    cc = cos(c)
    sd = sin(d)
    cd = cos(d)

    s2θ = sin(two * θ)
    c2θ = cos(two * θ)

    # chain pieces
    aθ = oneθ
    bθ = oneθ

    aκ = rκ * q1 + r * q1κv
    bκ = -(rκ * q2 + r * q2κv)
    cκ = (rκ * q2 + r * q2κv)
    dκ = rκ

    aμ = r * q1μv
    bμ = -(r * q2μv)
    cμ = r * q2μv
    dμ = zero(T)

    # -------------------------
    # F1 = sin(a) - L*sin(b)
    # -------------------------
    F1θ = ca * aθ - L * cb * bθ
    F1κ = ca * aκ - (Lκv * sb + L * cb * bκ)
    F1μ = ca * aμ - (Lμv * sb + L * cb * bμ)

    # -----------------------------------------------
    # F2 = sin(2θ) * sin(a)^2 + 2G*(1 + μ L^2)
    # -----------------------------------------------
    s_sq = s * s
    ds_dθ = ca * aθ
    ds_dκ = ca * aκ
    ds_dμ = ca * aμ

    # d/dθ sin(2θ) = 2 cos(2θ)
    F2θ = (two * c2θ) * s_sq + s2θ * (two * s * ds_dθ)

    # κ-derivative
    term_allowed = (oneθ + μ * (L * L))
    F2κ = s2θ * (two * s * ds_dκ) + (two * Gκv) * term_allowed + (two * G) * (μ * (two * L * Lκv))

    # μ-derivative
    F2μ = s2θ * (two * s * ds_dμ) + (two * Gμv) * term_allowed + (two * G) * ((L * L) + μ * (two * L * Lμv))

    # -------------------------------------------------------
    # F3 = cos(2θ)*sin(a)^3*sin(c) - 2G*cos(θ)*sin(√κ)
    # -------------------------------------------------------
    s_cu = s * s * s

    # θ-derivative first term:
    # d/dθ cos(2θ) = -2 sin(2θ)
    # d/dθ sin(a)^3 = 3 sin(a)^2 cos(a) * aθ
    F3θ = (-(two * s2θ)) * s_cu * sc + c2θ * (three * (s * s) * ca * aθ) * sc + (two * G) * sin(θ) * sd

    # κ-derivative first term:
    # d/dκ sin(a)^3 = 3 sin(a)^2 cos(a) * aκ
    # d/dκ sin(c)   = cos(c)*cκ
    F3κ = c2θ * (three * (s * s) * ca * aκ) * sc + c2θ * s_cu * (cc * cκ) - (two * Gκv) * cos(θ) * sd - (two * G) * cos(θ) * (cd * dκ)

    # μ-derivative (sin(√κ) has no μ dependence)
    F3μ = c2θ * (three * (s * s) * ca * aμ) * sc + c2θ * s_cu * (cc * cμ) - (two * Gμv) * cos(θ) * sd

    J = Matrix{T}(undef, 3, 3)
    J[1, 1] = F1θ
    J[1, 2] = F1κ
    J[1, 3] = F1μ
    J[2, 1] = F2θ
    J[2, 2] = F2κ
    J[2, 3] = F2μ
    J[3, 1] = F3θ
    J[3, 2] = F3κ
    J[3, 3] = F3μ
    return J
end


### Approximate zero
x0 = [-1.0471975511965979, 0.7444610429129391, 0.08867648509140362]
# F_appendix_D(x)

### Check derivative with ForwardDiff
# using ForwardDiff
# J_fd = ForwardDiff.jacobian(F_appendix_D, x0')
# J_eq = J_equilibrium(x0)


### Proof 
DF_interval = mid.(J_equilibrium(interval.(x0)))
A = inv(DF_interval)

x0_interval = interval.(x0)
F_interval = F_appendix_D(x0_interval)
Y = norm(interval.(A) * F_interval, Inf)

r0 = 1e-6
x0r_interval = interval.(x0, r0; format=:midpoint)
DF_interval = J_equilibrium(x0r_interval)
Z₁ = opnorm(LinearOperator(interval.(A) * DF_interval - UniformScaling(interval(1))), Inf)

interval_of_existence(Y, Z₁, r0; verbose=true)
