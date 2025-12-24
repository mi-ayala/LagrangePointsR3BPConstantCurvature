using Roots
using Plots
using PlotlyJS
plotlyjs()

include("functions_helpers.jl")
include("functions_f5.jl")
include("functions_stability.jl")
include("functions_plots.jl")
include("functions_grid.jl")

# bounds for the box
μT = 0.0886765
κ_lo, κ_hi = 0.0, (0.5 * pi)^2
# κ_lo, κ_hi = 0.8,1.0
μ_lo, μ_hi = 0.0, μT


# Example: make all grid pairs as a vector of tuples (κ, μ)
nk, nm = 250, 250
κs = collect(interior_linspace(κ_lo, κ_hi, nk))  # κ-grid
μs = collect(interior_linspace(μ_lo, μ_hi, nm))  # μ-grid

grid = [(κ, μ) for κ in κs, μ in μs]  # nk×nm matrix of pairs
roots_grid = zeros(nk * nm)  # storage for roots at each (κ, μ)
pairs = vec(grid)  # length = nk*nm


function pairs_to_roots(pairs; n=400, atol=1e-12, verbose=true)
    roots_1 = fill(NaN, length(pairs))
    roots_2 = fill(NaN, length(pairs))

    for i in eachindex(pairs)
        κ, μ = pairs[i]

        θlo, θhi = Intθ_5(κ, μ)

        br = brackets_signchanges(x -> f₅(x, κ, μ), θlo, θhi; n = 600)
        verbose && println("i=$i, (κ,μ)=($κ,$μ), brackets = ", br)

        nb = length(br)
        if nb >= 1
            roots_1[i] = find_zero(θ -> f₅(θ, κ, μ), br[1]; atol=atol)
        end
        if nb >= 2
            # println("Finding second root for i=$i, (κ,μ)=($κ,$μ)")
            roots_2[i] = find_zero(θ -> f₅(θ, κ, μ), br[2]; atol=atol)
        end
    end

    return roots_1, roots_2
end
roots_1, roots_2 = pairs_to_roots(pairs; n=400)

### Plotting f1
xs = first.(pairs)
ys = last.(pairs)            #

# collect sorted unique axes
κs = sort(unique(first.(pairs)))
μs = sort(unique(last.(pairs)))
nk, nm = length(κs), length(μs)


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


    # how many actual roots you got in roots_2
n2 = count(!isnan, roots_2)

# indices where roots_2 has a value
idx2 = findall(!isnan, roots_2)

# the actual values (no NaNs)
vals2 = roots_2[idx2]