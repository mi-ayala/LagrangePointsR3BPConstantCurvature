using Roots
using CairoMakie
using LaTeXStrings

# include("functions_helpers.jl")
# include("functions_validation.jl")
# include("functions_f5.jl")
include("functions_plots.jl")
include("functions_numerics.jl")

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


@inline twice(x) = x + x
@inline csc(x) = one(x) / sin(x)
@inline powabs3(x) = abs(x)^3
@inline sqr(x) = x * x
@inline three(x) = x + x + x
@inline inv3(x) = inv(three(one(x)))
@inline cbrt_t(x) = x^inv3(x)

function q₁(κ, μ)
    κ, μ = promote(κ, μ)
    r = sqrt(-κ)
    y = μ * sinh(twice(r))
    x = one(μ) + μ * cosh(twice(r))
    return (one(r) / twice(r)) * atanh(y / x)
end

function q₂(κ, μ)
    κ, μ = promote(κ, μ)
    return one(κ) - q₁(κ, μ)
end

function Γ(κ, μ)
    κ, μ = promote(κ, μ)
    r = sqrt(-κ)
    num = sqr(sinh(r)) * sinh(r) * cosh(r)
    den = sqrt(μ * μ + twice(μ * cosh(twice(r))) + one(μ))
    return num / den
end



# Λ(κ, μ) = ( sinh(-√κ q1) / ( μ sinh(-√κ q2) ) )^(1/3)
function Λ(κ, μ)
    κ, μ = promote(κ, μ)
    r = sqrt(-κ)
    q1 = q₁(κ, μ)
    q2 = q₂(κ, μ)

    s1 = sinh((-r) * q1)
    s2 = sinh((-r) * q2)

    return cbrt_t(s1 / (μ * s2))
end

# f(x, κ, μ) from the screenshot
function f₅(x, κ, μ)
    x, κ, μ = promote(x, κ, μ)

    r = sqrt(-κ)
    q1 = q₁(κ, μ)
    q2 = q₂(κ, μ)
    G = Γ(κ, μ)
    L = Λ(κ, μ)

    s1 = sinh((-r) * q1)
    s2 = sinh((-r) * q2)

    shx = sinh(x)
    chx = cosh(x)

    # sqrt(1 + Λ^2 sinh^2(x))
    inner = one(shx) + (L * L) * (shx * shx)
    rt = sqrt(inner)

    bracket = s1 * chx + s2 * rt

    pref = (L * L * L) * s2 * (shx * shx * shx)   # Λ^3 sinh(-r q2) sinh^3(x)
    tail = G * (sinh(-r) * sinh(-r))              # Γ sinh^2(-r)

    return pref * bracket - tail
end

# Numerics-only classification using the HYPERBOLIC test in your screenshot
# Same name + same inputs as your current function.
#
# Returns: label::Symbol
#   :center_saddle_unstable, :lyapunov_stable, :elliptic, :complex_saddle_unstable, or :unknown
#
function classify_center_equilibrium_num(x0::Real, κ::Real, μ::Real; verbose::Bool=false)
    x0 = float(x0)
    κ = float(κ)
    μ = float(μ)

    r = sqrt(-κ)
    q1 = q₁(κ, μ)
    q2 = q₂(κ, μ)

    # Λ = ( cosh(-r q2) / cosh(-r q1) )^(1/3)
    # (equivalently (sinh(-r q1)/(μ*sinh(-r q2)))^(1/3) depending on your definition,
    # but the screenshot shows the cosh ratio too.)
    Λ = cbrt(cosh((-r) * q2) / cosh((-r) * q1))

    # d1 = arcsinh( Λ sinh(x0) ),  d2 = x0
    d1 = asinh(Λ * sinh(x0))
    d2 = x0

    # Common pieces
    sh_r = sinh(-r)
    sh_r2 = sh_r * sh_r
    sh_r4 = sh_r2 * sh_r2

    sh_rq1 = sinh((-r) * q1)
    sh_rq2 = sinh((-r) * q2)

    ch_r = cosh(-r)

    ch_d1 = cosh(d1)
    ch_d2 = cosh(d2)

    # det(d1,d2) =
    #  [3 sinh(-r q1) sinh(-r q2) / sinh^4(-r)]
    #  * [sinh(-r q2) cosh(d1) + sinh(-r q1) cosh(d2)]
    #  * [(1+2 cosh(2 d2)) sinh(-r q1) cosh(d1) + (1+2 cosh(2 d1)) sinh(-r q2) cosh(d2)]
    pref = (3.0 * sh_rq1 * sh_rq2) / sh_r4

    B1 = sh_rq2 * ch_d1 + sh_rq1 * ch_d2
    B2 = (1.0 + 2.0 * cosh(2.0 * d2)) * sh_rq1 * ch_d1 +
         (1.0 + 2.0 * cosh(2.0 * d1)) * sh_rq2 * ch_d2

    det = pref * B1 * B2

    # tr(d1,d2) =
    #  1/sinh^2(-r) * [ (1-4 cosh^2(d2)) sinh^2(-r q1)
    #                  -6 cosh(d1)cosh(d2) sinh(-r q1)sinh(-r q2)
    #                  + (1-4 cosh^2(d1)) sinh^2(-r q2) ]
    tr = (1.0 / sh_r2) * (
        (1.0 - 4.0 * (ch_d2 * ch_d2)) * (sh_rq1 * sh_rq1)
        -
        6.0 * ch_d1 * ch_d2 * sh_rq1 * sh_rq2
        +
        (1.0 - 4.0 * (ch_d1 * ch_d1)) * (sh_rq2 * sh_rq2)
    )

    # m(d1,d2) = - sinh^2(d1) sinh^2(d2) / [(cosh(-r)-cosh(d1+d2)) (cosh(-r)-cosh(d1-d2))]
    sh_d1 = sinh(d1)
    sh_d2 = sinh(d2)
    num_m = (sh_d1 * sh_d1) * (sh_d2 * sh_d2)
    den_m = (ch_r - cosh(d1 + d2)) * (ch_r - cosh(d1 - d2))
    m = -(num_m / den_m)

    # C(d1,d2) = m(d1,d2) - 4 det(d1,d2)
    C = m - 4.0 * det

    verbose && println("Λ=$Λ  d1=$d1  d2=$d2  det=$det  tr=$tr  m=$m  C=$C")

    # Classification (as in screenshot)
    if !(isfinite(det) && isfinite(tr) && isfinite(C))
        return :unknown
    end

    if det < 0
        return :center_saddle_unstable
    elseif tr > 0
        return :lyapunov_stable
    elseif C > 0
        return :elliptic
    else
        return :complex_saddle_unstable
    end
end



function Intθ_5(κ, μ)
    κ, μ = promote(κ, μ)
    r = sqrt(-κ)
    L = Λ(κ, μ)

    num = sinh(-r)
    den = L + cosh(-r)

    return (atanh(num / den), 20.0)   # 20 is just a large number to ensure we cover the whole range
end


function cover_rectangle_grid_rect(κ_lo::Float64, κ_hi::Float64,
    μ_lo::Float64, μ_hi::Float64;
    nx::Int, ny::Int,
    verbose::Bool=true)

    @assert κ_hi > κ_lo "κ_hi must be > κ_lo"
    @assert μ_hi > μ_lo "μ_hi must be > μ_lo"
    @assert nx > 0 && ny > 0 "nx, ny must be positive"

    Lκ = κ_hi - κ_lo
    Lμ = μ_hi - μ_lo
    Δκ = Lκ / nx
    Δμ = Lμ / ny


    boxes = Tuple{Symbol,Float64,Float64,Float64,Float64,Symbol,Symbol}[]

    # iterate cell centers
    for i in 1:nx, j in 1:ny
        κc = κ_lo + (i - 0.5) * Δκ
        μc = μ_lo + (j - 0.5) * Δμ

        θmin = Intθ_5(κc, μc)[1]
        θmax = Intθ_5(κc, μc)[2]

        br = brackets_signchanges(θ -> f₅(θ, κc, μc), θmin, θmax)

        num_roots = length(br)

        if num_roots == 0
            push!(boxes, (:zero_roots, κc, μc, NaN, NaN, :zero_roots, :zero_roots))

        elseif num_roots == 1

            root = find_zero(θ -> f₅(θ, κc, μc), br[1]; atol=1e-12)
            stability_root = classify_center_equilibrium_num(root, κc, μc; verbose=false)
            push!(boxes, (:one_root, κc, μc, root, NaN, stability_root, :zero_roots))

        elseif num_roots == 2

            root1 = find_zero(θ -> f₅(θ, κc, μc), br[1]; atol=1e-12)
            stability_root1 = classify_center_equilibrium_num(root1, κc, μc; verbose=false)


            root2 = find_zero(θ -> f₅(θ, κc, μc), br[2]; atol=1e-12)
            stability_root2 = classify_center_equilibrium_num(root2, κc, μc; verbose=false)

            push!(boxes, (:two_roots, κc, μc, root1, root2, stability_root1, stability_root2))
        else
            push!(boxes, (:multiple_roots, κc, μc, NaN, NaN, :unknown, :unknown))


        end


    end


    return boxes
end


κ_lo, κ_hi = -0.5, -0.1
μ_lo, μ_hi = 0.0, 0.0886765

nx, ny = 1000, 1000


boxes = cover_rectangle_grid_rect(κ_lo, κ_hi, μ_lo, μ_hi; nx=nx, ny=ny, verbose=true)

# --- Root count plot ---
fig1, ax1 = plot_regions_all_labels_from_boxes_numerics(boxes;
    label_type=:root,
    ms=4,
)

save("figures/numerics/f5_root_regions_negative_numerics.png", fig1)

# --- Stability plot (root 1) ---
fig2, ax2 = plot_regions_all_labels_from_boxes_numerics(boxes;
    label_type=:stability,
    root_index=1,
    ms=4
)

save("figures/numerics/f5_stability_root1_negative_numerics.png", fig2)

# --- Stability plot (root 2) ---
fig3, ax3 = plot_regions_all_labels_from_boxes_numerics(boxes;
    label_type=:stability,
    root_index=2,
    ms=4
)

save("figures/numerics/f5_stability_root2_negative_numerics.png", fig3)


