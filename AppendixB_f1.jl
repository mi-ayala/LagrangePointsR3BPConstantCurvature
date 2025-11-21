using Roots
using Plots
using PlotlyJS
plotlyjs()

include("functions.jl")

# bounds for the box
μT = 0.0886765
κ_lo, κ_hi = 0.0, (0.5 * pi)^2
μ_lo, μ_hi = 0.0, μT

# make an interior grid (with margin)
interior_linspace(lo, hi, n; margin=0.001) = begin
    a = lo + margin * (hi - lo)
    b = hi - margin * (hi - lo)
    range(a, b; length=n)
end

# # make all grid pairs as a vector of tuples (κ, μ)
# nk, nm = 40, 40
# κs = collect(interior_linspace(κ_lo, κ_hi, nk))  # κ-grid
# μs = collect(interior_linspace(μ_lo, μ_hi, nm))  # μ-grid

# grid = [(κ, μ) for κ in κs, μ in μs]  # nk×nm matrix of pairs
# roots_grid = zeros(nk * nm)  # storage for roots at each (κ, μ)
# pairs = vec(grid)  # length = nk*nm

# function pairs_to_roots(pairs; n=400)
#     roots = zeros(length(pairs))
#     for i in 1:length(pairs)
#         κ, μ = pairs[i]
#         br = brackets_f₁(κ, μ; n=n)
#         println("i=$i, (κ,μ)=($κ,$μ), brackets = ", br)
#         if length(br) == 0
#             roots[i] = NaN  # no root found
#         else
#             rts = find_zero(θ -> f₁(θ, κ, μ), br[1])
#             roots[i] = rts[1]  # take first root
#         end
#     end
#     return roots
# end

# roots = pairs_to_roots(pairs; n=400)

# ### Plotting f1
# xs = first.(pairs)           
# ys = last.(pairs)            

# # collect sorted unique axes
# κs = sort(unique(first.(pairs)))
# μs = sort(unique(last.(pairs)))
# nk, nm = length(κs), length(μs)

# # fill θ*(κ_i, μ_j) matrix
# θgrid = fill(NaN, nk, nm)
# for (i, (κ, μ)) in enumerate(pairs)
#     ii = searchsortedfirst(κs, κ)
#     jj = searchsortedfirst(μs, μ)
#     θgrid[ii, jj] = roots[i]
# end

# Plots.surface(μs, κs, θgrid;
#     xlabel="μ", ylabel="κ", zlabel="θ*",
#     title="Root surface θ*(κ, μ)",
#     legend=false,
#     color=:cork, cbar_title="θ*",
#     linewidth=0,
#     fillalpha=0.95,
#     camera=(30, 60),
#     grid=:on, framestyle=:box)



# GLOBAL OR EXTERNAL STABILITY LOG
const STABILITY_LOG = Vector{NamedTuple}()
"""
    classify_equilibrium_log(θ, κ, μ)

Same as classify_equilibrium_fp, but also logs:
- detλ
- λ₂
- Δ
- label

Returns: (label, det, lam2, Δ)
"""
function classify_equilibrium_log(θ, κ, μ; verbose=false)
    # compute values
    det_val = detλ(θ, κ, μ)
    lam2_val = λ₂(θ, κ, μ)
    Δ_val   = Δ(θ, κ, μ)

    # classify
    label =
        det_val < 0        ? :saddle_hyperbolic :
        lam2_val > 0       ? :lyapunov_stable :
        Δ_val > 0          ? :elliptic :
        Δ_val < 0          ? :linearly_unstable :
                             :unknown

    # ---- LOG ENTRY ----
    push!(STABILITY_LOG, (
        θ      = θ,
        κ      = κ,
        μ      = μ,
        det    = det_val,
        lam2   = lam2_val,
        Δ      = Δ_val,
        label  = label
    ))

    return label, det_val, lam2_val, Δ_val
end


# create grid
nk, nm = 40, 40
κs = collect(interior_linspace(κ_lo, κ_hi, nk))
μs = collect(interior_linspace(μ_lo, μ_hi, nm))
pairs = [(κ, μ) for κ in κs, μ in μs] |> vec

# storage
results = Vector{Tuple{Float64,Float64,Float64,Symbol}}(undef, length(pairs))

function evaluate_grid(pairs; n=400)
    out = Vector{Tuple{Float64,Float64,Float64,Symbol}}(undef, length(pairs))

    for i in 1:length(pairs)
        κ, μ = pairs[i]

        br = brackets_f₁(κ, μ; n=n)
        println("i=$i, (κ,μ)=($κ,$μ), brackets = ", br)

        if isempty(br)
            # no root exists → store NaN + :no_root
            out[i] = (κ, μ, NaN, :no_root)
            continue
        end

        # numerical root
        θstar = find_zero(θ -> f₁(θ, κ, μ), br[1])

        # floating-point stability classification
        # stab, _, _, _ = classify_equilibrium_fp(θstar, κ, μ)
        stab, _, _, _ = classify_equilibrium_log(θstar, κ, μ)

        # store everything
        out[i] = (κ, μ, θstar, stab)
    end

    return out
end

results = evaluate_grid(pairs; n=400)

θgrid   = fill(NaN, nk, nm)
stabgrid = fill(Symbol(:none), nk, nm)

for (κ, μ, θstar, stab) in results
    ii = searchsortedfirst(κs, κ)
    jj = searchsortedfirst(μs, μ)
    θgrid[ii, jj] = θstar
    stabgrid[ii, jj] = stab
end

plt_surface = Plots.surface(
    μs, κs, θgrid;
    xlabel="μ", ylabel="κ", zlabel="θ*",
    title="Root surface θ*(κ, μ)",
    legend=false,
    color=:cork,
    cbar_title="θ*",
    linewidth=0,
    fillalpha=0.95,
    camera=(30, 60),
    grid=:on, framestyle=:box
)

using Plots

"""
    stability_color(sym::Symbol)

Map stability labels to colors for plotting.
"""
function stability_color(sym::Symbol)
    if sym === :saddle_hyperbolic
        return :red
    elseif sym === :lyapunov_stable
        return :green
    elseif sym === :elliptic
        return :blue
    elseif sym === :linearly_unstable
        return :orange
    elseif sym === :unknown
        return :gray
    elseif sym === :no_root
        return :black
    else
        return :black
    end
end

"""
    plot_stability_centers(results;
        ms=5,
        title="Stability in (κ, μ) plane",
        show_bounding_box=true)

`results` is a Vector of (κ, μ, θstar, stab::Symbol).

Makes a 2D scatter plot in the (κ, μ) plane with colors by stability.
"""
function plot_stability_centers(results;
        ms::Real=5,
        title::AbstractString="Stability in (κ, μ) plane",
        show_bounding_box::Bool=true)

    κ_vals = getindex.(results, 1)
    μ_vals = getindex.(results, 2)
    θ_vals = getindex.(results, 3)
    stabs  = getindex.(results, 4)

    # parameter bounds (adjust if needed)
    μT   = 0.0886765
    κ_lo, κ_hi = 0.0, (0.5 * π)^2
    μ_lo, μ_hi = 0.0, μT

    # Group by stability
    groups = Dict{Symbol, Vector{Tuple{Float64,Float64}}}()
    for (κ, μ, θ, stab) in results
        # skip points with NaN root if you like
        if isnan(θ)
            continue
        end
        push!(get!(groups, stab, Tuple{Float64,Float64}[]), (κ, μ))
    end

    plt = Plots.plot(; xlabel="κ", ylabel="μ",
                     title=title, legend=true,
                     xlims=(κ_lo, κ_hi), ylims=(μ_lo, μ_hi))

    for (stab, pts) in groups
        xs = first.(pts)
        ys = last.(pts)
        col = stability_color(stab)
        Plots.scatter!(plt, xs, ys;
            ms=ms, color=col, label=String(stab))
    end

    if show_bounding_box
        Xb = [κ_lo, κ_hi, κ_hi, κ_lo, κ_lo]
        Yb = [μ_lo, μ_lo, μ_hi, μ_hi, μ_lo]
        Plots.plot!(plt, Xb, Yb; color=:transparent, linecolor=:black, lw=2, label=false)
    end

    Plots.display(plt)
    return plt
end

κs = sort(unique(first.(pairs)))
μs = sort(unique(last.(pairs)))
θgrid = fill(NaN, nk, nm)
stabgrid = fill(Symbol(:none), nk, nm)

# filling θgrid and stabgrid:
for (κ, μ, θstar, stab) in results
    ii = searchsortedfirst(κs, κ)
    jj = searchsortedfirst(μs, μ)
    θgrid[ii, jj] = θstar
    stabgrid[ii, jj] = stab
end


"""
    plot_stability_surface(κs, μs, θgrid, stabgrid;
        surface_kwargs...)

Plot θ*(κ, μ) as a surface and overlay points colored by stability.

- `κs`, `μs` are 1D grids (length nk, nm)
- `θgrid` is nk×nm matrix of roots
- `stabgrid` is nk×nm matrix of Symbols (same layout)
"""
function plot_stability_surface(κs, μs, θgrid, stabgrid;
        title="Root surface with stability",
        surface_kwargs...)

    nk, nm = size(θgrid)
    @assert size(stabgrid) == size(θgrid)

    # Base surface of θ*
    plt = Plots.surface(
        μs, κs, θgrid;
        xlabel="μ", ylabel="κ", zlabel="θ*",
        title=title,
        legend=false,
        color=:cork,
        cbar_title="θ*",
        linewidth=0,
        fillalpha=0.95,
        camera=(30, 60),
        grid=:on, framestyle=:box,
        surface_kwargs...
    )

    # Overlay colored markers at each (κ, μ, θ*) with stable root
    for i in 1:nk, j in 1:nm
        θ = θgrid[i, j]
        stab = stabgrid[i, j]
        if isnan(θ) || stab == :no_root || stab == :none
            continue
        end
        κ = κs[i]
        μ = μs[j]
        col = stability_color(stab)
        Plots.scatter!(plt, [μ], [κ], [θ];
            ms=4, color=col)
    end

    Plots.display(plt)
    return plt
end

# # 2D κ–μ plot colored by stability
# plot_stability_centers(results)

# 3D surface θ*(κ, μ) with colored stability points
plot_stability_surface(κs, μs, θgrid, stabgrid)


κvals = getfield.(STABILITY_LOG, :κ)
μvals = getfield.(STABILITY_LOG, :μ)
Δvals = getfield.(STABILITY_LOG, :det)

Plots.scatter(κvals, μvals; zcolor=Δvals, ms=5, xlabel="κ", ylabel="μ", title="Δ values")