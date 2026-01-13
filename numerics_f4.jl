using Roots
using CairoMakie


include("functions_helpers.jl")
include("functions_validation.jl")
include("functions_f4.jl")
include("functions_plots.jl")
include("functions_numerics.jl")



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


        θmin = θ4_min(κc, μc)
        θmax = θ4_max(κc, μc)


        br = brackets_signchanges(θ -> f₄(θ, κc, μc), θmin, θmax)

        num_roots = length(br)

        if num_roots == 0
            push!(boxes, (:zero_roots, κc, μc, NaN, NaN, :zero_roots, :zero_roots))

        elseif num_roots == 1

            root = find_zero(θ -> f₄(θ, κc, μc), br[1]; atol=1e-12)
            stability_root = classify_equilibrium_f(root, κc, μc; verbose=false)
            push!(boxes, (:one_root, κc, μc, root, NaN, stability_root, :zero_roots))

        else
            push!(boxes, (:multiple_roots, κc, μc, NaN, NaN, :unknown, :unknown))


        end



    end


    return boxes
end


κ_lo, κ_hi = 0.0, (0.5 * pi)^2
μ_lo, μ_hi = 0.0, 0.0886765

nx, ny = 100, 100
boxes = cover_rectangle_grid_rect(κ_lo, κ_hi, μ_lo, μ_hi; nx=nx, ny=ny, verbose=true)

# --- Root count plot ---
fig1, ax1 = plot_regions_all_labels_from_boxes_numerics(boxes;
    label_type=:root,
    ms=4
)

save("figures/f4_root_regions_numerics.png", fig1)

# --- Stability plot (root 1) ---
fig2, ax2 = plot_regions_all_labels_from_boxes_numerics(boxes;
    label_type=:stability,
    root_index=1,
    ms=4
)

save("figures/f4_stability_numerics.png", fig2)

