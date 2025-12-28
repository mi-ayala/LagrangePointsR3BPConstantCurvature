using IntervalArithmetic
using Printf
using Roots
using CairoMakie

include("functions_helpers.jl")
include("functions_stability.jl")
include("functions_validation.jl")
include("functions_f2.jl")


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

        Iκ = interval(κc, rκ; format=:midpoint)
        Iμ = interval(μc, rμ; format=:midpoint)

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

            if root_label === :unknown
                push!(boxes, (true, true, κc, μc, rκ, rμ, root_label, :unknown))

            else
                push!(boxes, (true, true, κc, μc, rκ, rμ, root_label, stability_label))
                validated += 1
            end


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


κ_lo, κ_hi = 0.0, (0.5 * pi)^2
μ_lo, μ_hi = 0.0, 0.0886765

nx, ny = 2*2464, 2*89
boxes, stats = cover_rectangle_grid_rect(κ_lo, κ_hi, μ_lo, μ_hi; nx=nx, ny=ny, verbose=true)


fig1, ax = plot_regions_all_labels_from_boxes(
    boxes;
    label_type = :root,
    κlims = (0.0, (0.5*pi)^2),
    μlims = (0.0, 0.0886765),
    resolution = 2000,        
    alpha = 0.7,
    show_boundary = true,
    boundary_lw = 2.0,
    title = "Root-count Regions",
)

fig2, ax = plot_regions_all_labels_from_boxes(
    boxes;
    label_type = :stability,
    κlims = (0.0, (0.5*pi)^2),
    μlims = (0.0, 0.0886765),
    resolution = 2000,        
    alpha = 0.7,
    show_boundary = true,
    boundary_lw = 2.0,
    title = "Stability Regions",
)

CairoMakie.save("figures/f1_root_regions.png", fig1)  
CairoMakie.save("figures/f1_stab_regions.png", fig2)  




