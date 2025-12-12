using Roots
using Plots


include("functions.jl")
include("functions_f2.jl")


function brackets_f₂θ(κ, μ; n::Integer=400, margin::Real=0.0)
    θmin = Intθ_2(κ, μ)[1]
    θmax = Intθ_2(κ, μ)[2]
    δ = margin * (θmax - θmin)
    a, b = θmin + δ, θmax - δ
    brackets_signchanges(θ -> f₂θ(θ, κ, μ), a, b; n=n)
end


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
nk, nm = 100, 100
κs = collect(interior_linspace(κ_lo, κ_hi, nk))  # κ-grid
μs = collect(interior_linspace(μ_lo, μ_hi, nm))  # μ-grid

grid = [(κ, μ) for κ in κs, μ in μs] |> vec


function evaluate_grid(pairs; n::Int=400)
    N = length(pairs)
    results = Vector{Tuple{Float64,Float64,Float64}}(undef, N)
    log = NamedTuple[]

    for i in 1:N
        κ, μ = pairs[i]

        br = brackets_f₂θ(κ, μ; n=n)
        println("i=$i, (κ,μ)=($κ,$μ), brackets = ", br)

        if isempty(br)
            # no root: store NaN + :no_root and skip logging (or log if you like)
            results[i] = (κ, μ, NaN)

            continue
        end

        # numerical root
        θstar = find_zero(θ -> f₂θ(θ, κ, μ), br[1])


        # store compact result for plots
        results[i] = (κ, μ, θstar)


    end

    return results
end



function classify_f2_regions(results; tol::Real=0.0, verbose::Bool=false)
    out = Vector{Tuple{Float64,Float64,Float64,Float64,Symbol}}(undef, length(results))

    for (i, (κ, μ, θM)) in pairs(results)
        # Step 0: no maximum found
        if isnan(θM)
            out[i] = (κ, μ, θM, NaN, :no_max)
            verbose && println("i=$i κ=$κ μ=$μ → no maximum (θM=NaN)")
            continue
        end

        # Step 2: check θM ∈ [√κ q₂, π/2]
        left = sqrt(κ) * q₂(κ, μ)
        right = π / 2

        if !(left <= θM <= right)
            out[i] = (κ, μ, θM, NaN, :outside_range)
            verbose && println("i=$i κ=$κ μ=$μ θM=$θM outside [$left, $right] → zero roots")
            continue
        end

        # Step 3: evaluate f₂ at the maximum
        fval = f₂(θM, κ, μ)

        status =
            fval > tol ? :two_roots :
            fval < -tol ? :zero_roots :
            :unknown

        out[i] = (κ, μ, θM, fval, status)

        if verbose
            println("i=$i κ=$κ μ=$μ θM=$θM in [$left, $right], f₂(θM)=$fval → $status")
        end
    end

    return out
end


results = evaluate_grid(grid; n=400)
regions = classify_f2_regions(results; tol=1e-10, verbose=false)

using Plots

"""
    plot_f2_regions(regions; ms=5, title="f₂ root-count regions")

`regions` is a Vector of (κ, μ, θM, f2val, status).

Colors:
  :two_roots     → green
  :zero_roots    → red
  :outside_range → gray
  :no_max        → black
  :unknown       → orange
"""
function plot_f2_regions(regions;
    ms::Real=5,
    title::AbstractString="f₂ root-count regions")

    # group by status
    groups = Dict{Symbol,Vector{Tuple{Float64,Float64}}}()

    for (κ, μ, θM, fval, status) in regions
        push!(get!(groups, status, Tuple{Float64,Float64}[]), (κ, μ))
    end

    # color map
    color_of(status::Symbol) = status === :two_roots ? :green :
                               status === :zero_roots ? :red :
                               status === :outside_range ? :gray :
                               status === :no_max ? :black :
                               status === :unknown ? :orange :
                               :black

    plt = Plots.plot(; xlabel="κ", ylabel="μ",
        title=title, legend=true)

    for (status, pts) in groups
        xs = first.(pts)
        ys = last.(pts)
        col = color_of(status)
        Plots.scatter!(plt, xs, ys;
            ms=ms, color=col, label=String(status))
    end

    return plt
end


plt = plot_f2_regions(regions; ms=4, title="f₂: root-count regions")
display(plt)