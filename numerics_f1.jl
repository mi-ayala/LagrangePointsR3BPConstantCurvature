using Roots
using Plots


include("functions.jl")
include("functions_f1.jl")
include("stability_functions.jl")
include("plotting_functions_floats.jl")

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

# make all grid pairs as a vector of tuples (κ, μ)
nk, nm = 40, 40
κs = collect(interior_linspace(κ_lo, κ_hi, nk))  # κ-grid
μs = collect(interior_linspace(μ_lo, μ_hi, nm))  # μ-grid

grid = [(κ, μ) for κ in κs, μ in μs] |> vec

"""
    classify_equilibrium_fp(θ, κ, μ; verbose=false)

Floating-point version of the stability classifier.

Returns:
    (label::Symbol, det_val, lam2_val, Δ_val)
"""
function classify_equilibrium_f(θ, κ, μ; verbose::Bool=false)
    det_val = detλ(θ, κ, μ)
    lam2_val = λ₂(θ, κ, μ)
    Δ_val = Δ(θ, κ, μ)

    label =
        det_val < 0 ? :saddle_hyperbolic :
        lam2_val > 0 ? :lyapunov_stable :
        Δ_val > 0 ? :elliptic :
        Δ_val < 0 ? :linearly_unstable :
        :unknown

    if verbose
        println("θ=$θ κ=$κ μ=$μ → label=$label, det=$det_val, λ₂=$lam2_val, Δ=$Δ_val")
    end

    return label, det_val, lam2_val, Δ_val
end

"""
    evaluate_grid_f1(pairs; n=400)

For each (κ, μ) in `pairs`:
  - try to bracket a root of f₁,
  - if found, compute θ* with `find_zero`,
  - classify stability with `classify_equilibrium_fp`.

Returns:
    results, log

where

    results :: Vector{Tuple{Float64,Float64,Float64,Symbol}}
        (κ, μ, θ*, label)

    log :: Vector{NamedTuple}
        each entry: (θ, κ, μ, det, lam2, Δ, label)
"""
function evaluate_grid_f1(pairs; n::Int=400)
    N = length(pairs)
    results = Vector{Tuple{Float64,Float64,Float64,Symbol}}(undef, N)
    log = NamedTuple[]

    for i in 1:N
        κ, μ = pairs[i]

        br = brackets_f₁(κ, μ; n=n)
        println("i=$i, (κ,μ)=($κ,$μ), brackets = ", br)

        if isempty(br)
            # no root: store NaN + :no_root and skip logging (or log if you like)
            results[i] = (κ, μ, NaN, :no_root)
            push!(log, (θ=NaN, κ=κ, μ=μ,
                det=NaN, lam2=NaN, Δ=NaN, label=:no_root))
            continue
        end

        # numerical root
        θstar = find_zero(θ -> f₁(θ, κ, μ), br[1])

        # stability at that root
        label, det_val, lam2_val, Δ_val = classify_equilibrium_f(θstar, κ, μ)

        # store compact result for plots
        results[i] = (κ, μ, θstar, label)

        # store full diagnostic info in the log
        push!(log, (
            θ=θstar,
            κ=κ,
            μ=μ,
            det=det_val,
            lam2=lam2_val,
            Δ=Δ_val,
            label=label
        ))
    end

    return results, log
end


results, STABILITY_LOG = evaluate_grid_f1(grid; n=400)

θgrid = fill(NaN, nk, nm)
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

plot_stability_centers(results)

