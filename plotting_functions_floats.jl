using Printf
using Plots


#### Plotting stability functions 
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
    stabs = getindex.(results, 4)

    # parameter bounds (adjust if needed)
    μT = 0.0886765
    κ_lo, κ_hi = 0.0, (0.5 * π)^2
    μ_lo, μ_hi = 0.0, μT

    # Group by stability
    groups = Dict{Symbol,Vector{Tuple{Float64,Float64}}}()
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


