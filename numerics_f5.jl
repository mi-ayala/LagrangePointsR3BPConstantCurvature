using Roots
using CairoMakie


include("functions_helpers.jl")
include("functions_validation.jl")
include("functions_f5.jl")
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


κ_lo, κ_hi = 0.0, (0.5 * pi)^2
μ_lo, μ_hi = 0.0, 0.0886765

nx, ny = 1000, 1000


boxes = cover_rectangle_grid_rect(κ_lo, κ_hi, μ_lo, μ_hi; nx=nx, ny=ny, verbose=true)

# --- Root count plot ---
fig1, ax1 = plot_regions_all_labels_from_boxes_numerics(boxes;
    label_type = :root,
    ms = 4,
)

save("figures/numerics/f5_root_regions_numerics.png", fig1)

# --- Stability plot (root 1) ---
fig2, ax2 = plot_regions_all_labels_from_boxes_numerics(boxes;
    label_type = :stability,
    root_index = 1,
        ms = 4
)

save("figures/numerics/f5_stability_root1_numerics.png", fig2)

# --- Stability plot (root 2) ---
fig3, ax3 = plot_regions_all_labels_from_boxes_numerics(boxes;
    label_type = :stability,
    root_index = 2,
        ms = 4
)

save("figures/numerics/f5_stability_root2_numerics.png", fig3)




