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
interior_linspace(lo, hi, n; margin=0.007) = begin
    a = lo + margin * (hi - lo)
    b = hi - margin * (hi - lo)
    range(a, b; length=n)
end

# Example: make all grid pairs as a vector of tuples (κ, μ)
nk, nm = 90, 90
κs = collect(interior_linspace(κ_lo, κ_hi, nk))  # κ-grid
μs = collect(interior_linspace(μ_lo, μ_hi, nm))  # μ-grid

grid = [(κ, μ) for κ in κs, μ in μs]  # nk×nm matrix of pairs
roots_grid = zeros(nk * nm)  # storage for roots at each (κ, μ)
pairs = vec(grid)  # length = nk*nm


function pairs_to_roots(pairs; n=400)
    roots_1 = zeros(length(pairs))
    roots_2 = zeros(length(pairs))


    for i in 1:length(pairs)
        κ, μ = pairs[i]
        br = brackets_f₂(κ, μ; n=n)
        println("i=$i, (κ,μ)=($κ,$μ), brackets = ", br)
        if length(br) == 0
            roots_1[i] = NaN  # no root found
            roots_2[i] = NaN
        else
            rts = find_zero(θ -> f₂(θ, κ, μ), br[1])
            roots_1[i] = rts[1]  # take first root
            rts = find_zero(θ -> f₂(θ, κ, μ), br[2])
            roots_2[i] = rts[1]  # take first root
        end
    end
    return (roots_1, roots_2)
end
roots_1, roots_2 = pairs_to_roots(pairs; n=400)

### Plotting f1
xs = first.(pairs)
ys = last.(pairs)            #

# collect sorted unique axes
κs = sort(unique(first.(pairs)))
μs = sort(unique(last.(pairs)))
nk, nm = length(κs), length(μs)

# fill θ*(κ_i, μ_j) matrix
θgrid = fill(NaN, nk, nm)
for (i, (κ, μ)) in enumerate(pairs)
    ii = searchsortedfirst(κs, κ)
    jj = searchsortedfirst(μs, μ)
    θgrid[ii, jj] = roots_1[i]
end

θgrid2 = fill(NaN, nk, nm)
for (i, (κ, μ)) in enumerate(pairs)
    ii = searchsortedfirst(κs, κ)
    jj = searchsortedfirst(μs, μ)
    θgrid2[ii, jj] = roots_2[i]
end

Plots.surface(μs, κs, θgrid;
    xlabel="μ", ylabel="κ", zlabel="θ*",
    title="Root surface θ*(κ, μ)",
    legend=false,
    color=:cork, cbar_title="θ*",
    linewidth=0,
    fillalpha=0.95,
    camera=(30, 60),
    grid=:on, framestyle=:box)


Plots.surface!(μs, κs, θgrid2;
    color=:cork, cbar_title="θ*",
    linewidth=0,
    fillalpha=0.95,
    camera=(30, 60),
    grid=:on, framestyle=:box)

"""
    evaluate_grid_f2(pairs; n=400)

Compute up to two roots of f₂ at each (κ, μ) point and classify stability
using `classify_equilibrium_fp`.

Returns a vector of tuples:
    (κ, μ, θ1, stab1, θ2, stab2)

with the convention that, when both exist,
    θ1 ≤ θ2   (i.e. root 1 is always the smaller root).
"""
function evaluate_grid_f2(pairs; n=400)
    out = Vector{Tuple{Float64,Float64,Float64,Symbol,Float64,Symbol}}(undef, length(pairs))

    for i in 1:length(pairs)
        κ, μ = pairs[i]

        br = brackets_f₂(κ, μ; n=n)
        println("i=$i, (κ,μ)=($κ,$μ), brackets = ", br)

        # default: no roots
        θ1 = NaN
        stab1 = :no_root
        θ2 = NaN
        stab2 = :no_root

        if isempty(br)
            # keep defaults
            out[i] = (κ, μ, θ1, stab1, θ2, stab2)
            continue
        end

        # always take the first bracket
        θa = find_zero(θ -> f₂(θ, κ, μ), br[1])
        staba, _, _, _ = classify_equilibrium_fp(θa, κ, μ)

        if length(br) == 1
            # only one root we trust
            θ1, stab1 = θa, staba
            out[i] = (κ, μ, θ1, stab1, θ2, stab2)
            continue
        end

        # second root from second bracket
        θb = find_zero(θ -> f₂(θ, κ, μ), br[2])
        stabb, _, _, _ = classify_equilibrium_fp(θb, κ, μ)

        # reorder so that θ1 ≤ θ2
        if θa ≤ θb
            θ1, stab1 = θa, staba
            θ2, stab2 = θb, stabb
        else
            θ1, stab1 = θb, stabb
            θ2, stab2 = θa, staba
        end

        out[i] = (κ, μ, θ1, stab1, θ2, stab2)
    end

    return out
end

results₂ = evaluate_grid_f2(pairs; n=400)

κs = sort(unique(first.(pairs)))
μs = sort(unique(last.(pairs)))
nk, nm = length(κs), length(μs)

θgrid1 = fill(NaN, nk, nm)
θgrid2 = fill(NaN, nk, nm)

stabgrid1 = fill(Symbol(:none), nk, nm)
stabgrid2 = fill(Symbol(:none), nk, nm)

for (κ, μ, θ₁, stab₁, θ₂, stab₂) in results₂
    ii = searchsortedfirst(κs, κ)
    jj = searchsortedfirst(μs, μ)
    θgrid1[ii, jj] = θ₁
    θgrid2[ii, jj] = θ₂
    stabgrid1[ii, jj] = stab₁
    stabgrid2[ii, jj] = stab₂
end



function plot_stability_centers_f2(results; root_index=1,
    ms=4, title="Stability (f₂)")
    κvals = Float64[]
    μvals = Float64[]
    cols = Symbol[]

    for (κ, μ, θ₁, stab₁, θ₂, stab₂) in results
        stab = (root_index == 1) ? stab₁ : stab₂
        if !isnan((root_index == 1 ? θ₁ : θ₂))
            push!(κvals, κ)
            push!(μvals, μ)
            push!(cols, stability_color(stab))
        end
    end

    plt = Plots.scatter(κvals, μvals;
        ms=ms, color=cols,
        xlabel="κ", ylabel="μ",
        title=title, legend=false)

    return plt
end
plot_stability_centers_f2(results₂; root_index=1,
    title="Stability for f₂ — first branch")

plot_stability_centers_f2(results₂; root_index=2,
    title="Stability for f₂ — second branch")