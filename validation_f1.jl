using IntervalArithmetic
using Printf
using Roots

include("functions_helpers.jl")
include("functions_stability.jl")
include("functions_validation.jl")
include("functions_f1.jl")


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
    rκ = Δκ / 2
    rμ = Δμ / 2

    boxes = Tuple{Bool,Bool,Float64,Float64,Float64,Float64,Symbol,Symbol}[]
    total = nx * ny
    seeded = 0
    validated = 0

    # iterate cell centers
    for i in 1:nx, j in 1:ny
        κc = κ_lo + (i - 0.5) * Δκ
        μc = μ_lo + (j - 0.5) * Δμ

        (θmin, θmax) = Intθ_1(κc, μc)

        br = brackets_signchanges(θ -> f₁(θ, κc, μc), θmin, θmax)

        if isempty(br)
            continue
        end
        seeded += 1 # Seeded means we found a root 

        # 2) refine the root 
        root = find_zero(θ -> f₁(θ, κc, μc), br[1]; atol=1e-12)

        Iκ = interval(κ_lo, κ_hi)
        Iμ = interval(μ_lo, μ_hi)

        θmax_refined = sup(θ1_max(Iκ, Iμ))

        ok_root, data = validate_root_sign_mono(
            f₁,
            f₁θ,
            root,
            interval(0.0, θmax_refined),
            Iκ,
            Iμ;
            ntiles=16,
            min_width=1e-12,
            min_radius=1e-12,
            verbose=false,
        )

        keep_invalid = true

        if ok_root

            _, root_label = finalize_root_label(complements_nonzero(f₁, interval(0.0, θmax_refined), data.J, Iκ, Iμ))

            _, stability_label = classify_equilibrium(data.Iexist, Iκ, Iμ)

            push!(boxes, (true, true, κc, μc, rκ, rμ, root_label, stability_label))
            validated += 1
            
        elseif keep_invalid
            push!(boxes, (false, false, κc, μc, rκ, rμ, :unknown, :unknown))
        end

    end

    if verbose
        println("cover_rectangle_grid_rect:")
        println("  Δκ=$(Δκ), Δμ=$(Δμ); rκ=$(rκ), rμ=$(rμ)")
        println("  cells total=$(total), seeded=$(seeded), validated=$(validated)")
    end

    stats = (Δκ=Δκ, Δμ=Δμ, rκ=rκ, rμ=rμ, total=total, seeded=seeded, validated=validated)
    return boxes, stats
end



κ_lo, κ_hi = 0.5, (0.5 * pi)^2
μ_lo, μ_hi = 0.0, 0.0886765

nx, ny = 2464, 89
# nx, ny = 2500, 100
boxes, stats = cover_rectangle_grid_rect(κ_lo, κ_hi, μ_lo, μ_hi; nx=nx, ny=ny, verbose=true)


# ### Plotting

# # ============================================================
# # Color maps for CairoMakie / Plots
# # ============================================================

# # ---------- Existence / root-count colors ----------
# @inline function rootcount_color(lbl::Symbol)
#     return (
#         lbl === :two_roots  ? "#2C5AA0" :   # strong blue
#         lbl === :one_root   ? "#8FBCE6" :   # light blue
#         lbl === :zero_roots ? "#D0D0D0" :   # light gray
#         lbl === :unknown    ? "#ED1B2F" :   # red (failed / unknown)
#                               "#000000"    # fallback
#     )
# end


# # ---------- Stability colors (stable → unstable) ----------
# @inline function stability_color(lbl::Symbol)
#     return (
#         lbl === :lyapunov_stable   ? "#1B9E77" :  # green
#         lbl === :elliptic          ? "#A6D854" :  # yellow-green
#         lbl === :linearly_unstable ? "#E67E22" :  # orange
#         lbl === :saddle_hyperbolic ? "#B35806" :  # dark orange / brown
#         lbl === :unknown           ? "#ED1B2F" :  # red
#                                     "#000000"
#     )
# end

# using CairoMakie