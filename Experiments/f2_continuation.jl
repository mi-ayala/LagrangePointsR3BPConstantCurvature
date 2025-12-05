using RadiiPolynomial
using LinearAlgebra
# import LinearAlgebra: ⋅


include("functions.jl")

# --------------------------------------------------------
# q₁μ
# --------------------------------------------------------
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

# --------------------------------------------------------
# q₁κ
# --------------------------------------------------------
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

# --------------------------------------------------------
# Γμ
# --------------------------------------------------------
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

# --------------------------------------------------------
# Γκ
# --------------------------------------------------------
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

# --------------------------------------------------------
# f₂
# --------------------------------------------------------
function f₂(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    r = sqrt(κ)
    q1 = q₁(κ, μ)
    q2 = q₂(κ, μ)

    A = θ + r * q1
    B = θ - r * q2

    term = inv(sqr(sin(A))) + μ * inv(sqr(sin(B)))
    return sin(θ) * cos(θ) - Γ(κ, μ) * term
end

# --------------------------------------------------------
# f₂θ
# --------------------------------------------------------
function f₂θ(θ, κ, μ)
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

    FA = cA / (sA * sA * sA)
    FB = cB / (sB * sB * sB)

    return cos(twice(θ)) + two * Γv * (FA + μ * FB)
end

# --------------------------------------------------------
# f₂κ
# --------------------------------------------------------
function f₂κ(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)

    oneT = one(θ)
    two = twice(oneT)

    r = sqrt(κ)
    dr = one(r) / (two * r)

    q1 = q₁(κ, μ)
    q2 = q₂(κ, μ)
    dq1 = q₁κ(κ, μ)
    dq2 = -dq1

    Γv = Γ(κ, μ)
    Γκv = Γκ(κ, μ)

    A = θ + r * q1
    B = θ - r * q2

    sA = sin(A)
    cA = cos(A)
    sB = sin(B)
    cB = cos(B)

    inv_sA2 = inv(sA * sA)
    inv_sB2 = inv(sB * sB)

    T0 = inv_sA2 + μ * inv_sB2

    Aκ = dr * q1 + r * dq1
    Bκ = -(dr * q2 + r * dq2)

    FAκ = -(two * cA / (sA * sA * sA)) * Aκ
    FBκ = μ * (-(two * cB / (sB * sB * sB)) * Bκ)

    return -Γκv * T0 - Γv * (FAκ + FBκ)
end

# --------------------------------------------------------
# f₂μ
# --------------------------------------------------------
function f₂μ(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)

    oneT = one(θ)
    two = twice(oneT)

    r = sqrt(κ)

    q1 = q₁(κ, μ)
    q2 = q₂(κ, μ)
    dq1 = q₁μ(κ, μ)
    dq2 = -dq1

    Γv = Γ(κ, μ)
    Γμv = Γμ(κ, μ)

    A = θ + r * q1
    B = θ - r * q2

    sA = sin(A)
    cA = cos(A)
    sB = sin(B)
    cB = cos(B)

    inv_sA2 = inv(sA * sA)
    inv_sB2 = inv(sB * sB)

    T0 = inv_sA2 + μ * inv_sB2

    Aμ = r * dq1
    Bμ = r * dq1

    FAμ = -(two * cA / (sA * sA * sA)) * Aμ
    FBμ = μ * (-(two * cB / (sB * sB * sB)) * Bμ) + inv_sB2

    return -Γμv * T0 - Γv * (FAμ + FBμ)
end

# --------------------------------------------------------
# f₂θθ
# --------------------------------------------------------
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

# --------------------------------------------------------
# f₂θκ
# --------------------------------------------------------
function f₂θκ(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    oneT = one(θ)
    two = twice(oneT)

    r = sqrt(κ)
    dr = one(r) / (two * r)

    q1 = q₁(κ, μ)
    q2 = q₂(κ, μ)
    dq1 = q₁κ(κ, μ)
    dq2 = -dq1

    Γv = Γ(κ, μ)
    Γκv = Γκ(κ, μ)

    A = θ + r * q1
    B = θ - r * q2

    sA = sin(A)
    cA = cos(A)
    sB = sin(B)
    cB = cos(B)

    inv_sA3 = inv(sA * sA * sA)
    inv_sB3 = inv(sB * sB * sB)

    FA = cA * inv_sA3
    FB = cB * inv_sB3

    inv_sA4 = inv(sA * sA * sA * sA)
    inv_sB4 = inv(sB * sB * sB * sB)

    two_cA2 = two * cA * cA
    two_cB2 = two * cB * cB

    FA′ = -(oneT + two_cA2) * inv_sA4
    FB′ = -(oneT + two_cB2) * inv_sB4

    Aκ = dr * q1 + r * dq1
    Bκ = -(dr * q2 + r * dq2)

    FAκ = FA′ * Aκ
    FBκ = FB′ * Bκ

    return two * (Γκv * (FA + μ * FB) + Γv * (FAκ + μ * FBκ))
end

# --------------------------------------------------------
# f₂θμ
# --------------------------------------------------------
function f₂θμ(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    oneT = one(θ)
    two = twice(oneT)

    r = sqrt(κ)

    q1 = q₁(κ, μ)
    q2 = q₂(κ, μ)
    dq1 = q₁μ(κ, μ)
    dq2 = -dq1

    Γv = Γ(κ, μ)
    Γμv = Γμ(κ, μ)

    A = θ + r * q1
    B = θ - r * q2

    sA = sin(A)
    cA = cos(A)
    sB = sin(B)
    cB = cos(B)

    inv_sA3 = inv(sA * sA * sA)
    inv_sB3 = inv(sB * sB * sB)

    FA = cA * inv_sA3
    FB = cB * inv_sB3

    inv_sA4 = inv(sA * sA * sA * sA)
    inv_sB4 = inv(sB * sB * sB * sB)

    two_cA2 = two * cA * cA
    two_cB2 = two * cB * cB

    FA′ = -(oneT + two_cA2) * inv_sA4
    FB′ = -(oneT + two_cB2) * inv_sB4

    Aμ = r * dq1
    Bμ = r * dq1

    FAμ = FA′ * Aμ
    FBμ = FB′ * Bμ

    return two * (Γμv * (FA + μ * FB) + Γv * (FAμ + μ * FBμ + FB))
end


function J_f₂(θ, κ, μ)
    θ, κ, μ = promote(θ, κ, μ)
    T = typeof(one(θ))
    J = Matrix{T}(undef, 2, 3)

    J[1, 1] = f₂θ(θ, κ, μ)
    J[1, 2] = f₂κ(θ, κ, μ)
    J[1, 3] = f₂μ(θ, κ, μ)

    J[2, 1] = f₂θθ(θ, κ, μ)
    J[2, 2] = f₂θκ(θ, κ, μ)
    J[2, 3] = f₂θμ(θ, κ, μ)

    return J
end


### u is (theta, kappa, mu)

mu = 0.02
kappa = 0.5
theta = 1.2

g = (x) -> [f₂(x[1], x[2], x[3]); f₂θ(x[1], x[2], x[3])]
Dᵤg = (x) -> J_f₂(x[1], x[2], x[3])


# μ= 0.05


# f₂(u, γ, μ)
# # Dᵤf=(u, γ) -> f₂θ(u, γ, μ)
# # ∂γf=(u, γ) -> f₂μ(u, κ, γ) 
# f = (u, γ) -> [f₂(u, κ, μ); f₂θ(u, κ, μ)]



f = (u, γ) -> [f₂(u[1], u[2], mu); f₂θ(u[1], u[2], mu)]
Dᵤf = (u, γ) -> J_f₂(u[1], u[2], mu)[1:2, 1:2]


F(x, v, w) = [f(x[1:2], x[3]); (x - w) ⋅ v]
DF(x, v) = [J_f₂(x[1], x[2], x[3]); transpose(v)]



# # initial point on the branch of equilibria
# γ_init = 0.02
# u_init = [theta, kappa]

# f(u_init, γ_init)
# Dᵤf(u_init, γ_init)

# u_init, success = newton(u -> (f(u, γ_init), Dᵤf(u, γ_init)), u_init;verbose=true)

# # next point on the branch of equilibria
# x_init = [u_init ; γ_init]

# v = vec(nullspace(J_f₂(x_init[1], x_init[2], x_init[3]))) # initial tangent vector

# δ = 5e-6 # step size

# w = x_init + δ * v # predictor

# w

#         # --- UPDATE nonlinear system for this γ ---
#         f  = make_f(w[3])
#         Df = make_Du_f(w[3])
#         Ffun  = make_F(f, v, w)
#         DFfun = make_DF(w[3], v)

# x_final, success = newton(x -> (Ffun(x), DFfun(x)), w; verbose=true, maxiter=50)


# f(u, γ) = system in (θ, κ)
make_f(γ) = (u -> @SVector [f₂(u[1], u[2], γ),
    f₂θ(u[1], u[2], γ)])

# Dᵤ f(u, γ) = 2×2 Jacobian wrt u = (θ, κ)
make_Du_f(γ) = (u -> J_f₂(u[1], u[2], γ)[1:2, 1:2])

# Full pseudo-arclength system
make_F(f, v, w) = (x -> vcat(f(x[1:2]), dot(x - w, v)))

# Full Jacobian of F
make_DF(γ, v) = (x -> vcat(J_f₂(x[1], x[2], γ),
    transpose(v)))

function continue_branch(x_init;
    γ_target::Real,
    δ::Real=5e-4,
    max_steps::Int=100,
    newton_kwargs=(; verbose=false, maxiter=50, tol=1e-8),
    verbose::Bool=true)

    x = copy(x_init)                 # [θ, κ, γ]
    path = [copy(x)]
    γ_idx = 3                        # index of parameter γ=μ

    # direction of continuation in γ
    γ0 = x[γ_idx]
    Δγ = γ_target - γ0
    sgn = sign(Δγ)
    sgn == 0 && error("Already at γ_target")

    # Build f and Du f at initial γ
    f = make_f(x[γ_idx])
    Df = make_Du_f(x[γ_idx])

    # Initial tangent from nullspace
    J = J_f₂(x[1], x[2], x[3])
    v = vec(nullspace(J))
    if sign(v[γ_idx]) != sgn
        v .= -v
    end

    for step in 1:max_steps
        γ = x[γ_idx]
        if sgn * (γ_target - γ) ≤ 0
            verbose && println("Reached γ_target at γ=$(x[γ_idx])")
            break
        end

        # --- PREDICTOR ---
        w = x .+ δ * v

        # --- UPDATE nonlinear system for this γ ---
        f = make_f(w[γ_idx])
        Df = make_Du_f(w[γ_idx])
        Ffun = make_F(f, v, w)
        DFfun = make_DF(w[γ_idx], v)

        # --- CORRECTOR ---
        x_new, success = newton((x -> (Ffun(x), DFfun(x))), w; newton_kwargs...)
        if !success
            verbose && println("Newton failed at step $step")
            break
        end

        push!(path, copy(x_new))

        # Update tangent
        J_new = J_f₂(x_new[1], x_new[2], x_new[3])
        v_new = vec(nullspace(J_new))
        if dot(v_new, v) < 0
            v_new .= -v_new
        end

        x .= x_new
        v .= v_new
    end

    return path
end

# 1. You already did:
γ_init = 0.02
u_init = [theta, kappa]

# refine initial equilibrium in u:
u_init, success = newton(u -> (f(u, γ_init), Dᵤf(u, γ_init)), u_init; verbose=true)

# build initial x = [θ, κ, γ]
x_init = [u_init; γ_init]

# 2. Run continuation up to some γ_target:
path = continue_branch(x_init; γ_target=0.0886765, δ=5e-6, max_steps=30000)

# # 3. Extract coordinates for plotting:
θs = [p[1] for p in path]
κs = [p[2] for p in path]
γs = [p[3] for p in path]

# # Example Plots.jl usage:
using Plots
Plots.plot(γs, κs; xlabel="γ", ylabel="κ", title="Continuation branch (κ vs γ)")