using ForwardDiff
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

# Example: make all grid pairs as a vector of tuples (κ, μ)
nk, nm = 20, 20
κs = collect(interior_linspace(κ_lo, κ_hi, nk))  # κ-grid
μs = collect(interior_linspace(μ_lo, μ_hi, nm))  # μ-grid

grid = [(κ, μ) for κ in κs, μ in μs]  # nk×nm matrix of pairs
roots_grid = zeros(nk * nm)  # storage for roots at each (κ, μ)
pairs = vec(grid)  # length = nk*nm

function pairs_to_roots(pairs; n=400)
    roots = zeros(length(pairs))
    for i in 1:length(pairs)
        κ, μ = pairs[i]
        br = brackets_f₁(κ, μ; n=n)
        println("i=$i, (κ,μ)=($κ,$μ), brackets = ", br)
        if length(br) == 0
            roots[i] = NaN  # no root found
        else
            rts = find_zero(θ -> f₁(θ, κ, μ), br[1])
            roots[i] = rts[1]  # take first root
        end
    end
    return roots
end

roots = pairs_to_roots(pairs; n=400)

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
    θgrid[ii, jj] = roots[i]
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




