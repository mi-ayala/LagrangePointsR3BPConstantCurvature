using IntervalArithmetic
using Printf
using Roots
using CairoMakie

include("functions_helpers.jl")
include("functions_validation.jl")
include("functions_f4.jl")
include("functions_plots.jl")


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

    boxes = Tuple{Symbol,Float64,Float64,Float64,Float64,Bool,Bool,Symbol,Symbol}[]
    total = nx * ny
    seeded = 0
    validated = 0

    # iterate cell centers
    for i in 1:nx, j in 1:ny
        κc = κ_lo + (i - 0.5) * Δκ
        μc = μ_lo + (j - 0.5) * Δμ

        (θmin, θmax) = Intθ_4(κc, μc)

        br = brackets_signchanges(θ -> f₄(θ, κc, μc), θmin, θmax)

        if isempty(br)
            continue
        end
        seeded += 1 # Seeded means we found a root 

        # 2) refine the root 
        root = find_zero(θ -> f₄(θ, κc, μc), br[1]; atol=1e-12)

        Iκ = interval(κc, rκ; format=:midpoint)
        Iμ = interval(μc, rμ; format=:midpoint)

        θmax_refined = sup(θ4_max(Iκ, Iμ))
        θmin_refined = inf(θ4_min(Iκ, Iμ))


        ok_root, data = validate_root_sign_mono(
            f₄,
            f₄θ,
            root,
            interval(θmin_refined, θmax_refined),
            Iκ,
            Iμ;
            ntiles=30,
            min_width=1e-12,
            min_radius=1e-12,
            verbose=false,
        )

        keep_invalid = true

        if ok_root

            _, root_label = finalize_root_label(complements_nonzero(f₄, interval(θmin_refined, θmax_refined), data.J, Iκ, Iμ))

            _, stability_label = classify_equilibrium(data.Iexist, Iκ, Iμ)

            if root_label === :unknown

                push!(boxes, (root_label, κc, μc, rκ, rμ, ok_root, false, :unknown, :zero_roots))

            else
                push!(boxes, (root_label, κc, μc, rκ, rμ, ok_root, false, stability_label, :zero_roots))
                validated += 1
            end


        elseif keep_invalid
            push!(boxes, (:unknown, κc, μc, rκ, rμ, ok_root, false, :unknown, :zero_roots))
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


κ_lo, κ_hi = 0.0001, (0.5 * pi)^2
μ_lo, μ_hi = 0.0001, 0.0886765

nx, ny =  2*2464, 2*89
boxes, stats = cover_rectangle_grid_rect(κ_lo, κ_hi, μ_lo, μ_hi; nx=nx, ny=ny, verbose=true)


# --- Root count plot ---
fig1, ax1 = plot_regions_all_labels_from_boxes_quick(boxes;
    label_type=:root,
    # title="Root-count",
)

save("/Users/miguelayala/Research/R3BP-on-spaces-with-constant-curvature/figures Miguel/validated/f4_root_regions.png", fig1)

# --- Stability plot (root 1) ---
fig2, ax2 = plot_regions_all_labels_from_boxes_quick(boxes;
    label_type=:stability,
    # title="Stability ",
)

save("/Users/miguelayala/Research/R3BP-on-spaces-with-constant-curvature/figures Miguel/validated/f4_stability.png", fig2)



